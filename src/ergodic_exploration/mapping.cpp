/**
 * @file mapping.cpp
 * @author Boston Cleek
 * @date 8 Nov 2020
 * @brief Occupancy grid mapping
 */

#include <cmath>

#include <ergodic_exploration/mapping.hpp>
#include <ergodic_exploration/numerics.hpp>

namespace ergodic_exploration
{
/**
 * @brief Convert log odds to a probability
 * @param l - log odds
 * @return probability
 */
double logOdds2Prob(double l)
{
  return 1.0 - (1.0 / (1.0 + std::exp(l)));
}

/**
 * @brief Convert probability to log odds
 * @param p - probability
 * @return log odds
 */
double prob2LogOdds(double p)
{
  return std::log(p / (1.0 - p));
}

OccupancyMapper::OccupancyMapper(const mat& Tbs)
  : Tbs_(Tbs)
  , prior_(0.5)
  , prob_occ_(0.90)
  , prob_free_(0.35)
  , log_odds_prior_(prob2LogOdds(prior_))
  , log_odds_occ_(prob2LogOdds(prob_occ_))
  , log_odds_free_(prob2LogOdds(prob_free_))
{
}

bool OccupancyMapper::updateMap(GridMap& grid,
                                const sensor_msgs::LaserScan::ConstPtr& scan,
                                const vec& pose)
{
  const auto psg = grid.world2Grid(pose(0), pose(1));
  if (!grid.gridBounds(psg.at(0), psg.at(1)))
  {
    std::cout << "WARNING: Cannot update map robot is not within bounds" << std::endl;
    return false;
  }

  // pm = Tmb * Tbs * ps
  const mat Tms = transform2d(pose(0), pose(1), pose(2)) * Tbs_;

  auto beam_angle = scan->angle_min;

  // Itertate over all measurements
  for (unsigned int i = 0; i < scan->ranges.size(); i++)
  {
    // Transform laser end point into frame of map
    const vec pt = Tms * polar2CartesianHomo(beam_angle, scan->ranges.at(i));
    const auto ptg = grid.world2Grid(pt(0), pt(1));

    // Only map using range measurements within map bounds
    if (grid.gridBounds(ptg.at(0), ptg.at(1)))
    {
      // Ray trace to find cell indices that are free space
      rayTrace(grid, psg.at(1), psg.at(0), ptg.at(1), ptg.at(0));

      // Update the end point
      updateCell(grid, ptg.at(1), ptg.at(0), log_odds_occ_);
    }

    else
    {
      std::cout << "WARNING: Skipping range measurement, not within bounds" << std::endl;
    }

    beam_angle += scan->angle_increment;
  }

  return true;
}

void OccupancyMapper::rayTrace(GridMap& grid, int x0, int y0, int x1, int y1) const
{
  auto dx = x1 - x0;
  auto dy = y1 - y0;

  // Case: vertical
  /////////////////////////////////////////////////////////////////////////////
  if (dx == 0)
  {
    // std::cout << "Plot Line Vertically" << std::endl;
    // move down
    if (dy < 0)
    {
      for (auto y = y0; y > y1; y--)
      {
        updateCell(grid, x0, y, log_odds_free_);
      }
    }

    // move up
    else
    {
      for (auto y = y0; y < y1; y++)
      {
        updateCell(grid, x0, y, log_odds_free_);
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////////

  // Case: horizontal
  /////////////////////////////////////////////////////////////////////////////
  else if (dy == 0)
  {
    // std::cout << "Plot Line Horizontally" << std::endl;
    // move left
    if (dx < 0)
    {
      for (auto x = x0; x > x1; x--)
      {
        updateCell(grid, x, y0, log_odds_free_);
      }
    }

    // move right
    else
    {
      for (auto x = x0; x < x1; x++)
      {
        updateCell(grid, x, y0, log_odds_free_);
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////////

  // Case: Draw line down and need to move farther in x than in y. (increment x every itertation)
  /////////////////////////////////////////////////////////////////////////////
  else if (std::abs(dy) < std::abs(dx))
  {
    // std::cout << "Plot Line Low" << std::endl;
    // Octant 5
    if (x0 > x1)
    {
      // add starting cell
      updateCell(grid, x0, y0, log_odds_free_);
      lineLow(grid, x1, y1, x0, y0);
    }

    // Octant 3
    else
    {
      // add starting cell
      updateCell(grid, x0, y0, log_odds_free_);
      lineLow(grid, x0, y0, x1, y1);
    }
  }
  /////////////////////////////////////////////////////////////////////////////

  // Case: Draw line down and need to move farther in y than in x. (increment y every itertation)
  /////////////////////////////////////////////////////////////////////////////
  else if (std::abs(dy) > std::abs(dx))
  {
    // std::cout << "Plot Line High" << std::endl;
    // Octant 7
    if (y0 > y1)
    {
      // add starting cell
      updateCell(grid, x0, y0, log_odds_free_);
      lineHigh(grid, x1, y1, x0, y0);
    }

    // Octant 1
    else
    {
      // add starting cell
      updateCell(grid, x0, y0, log_odds_free_);
      lineHigh(grid, x0, y0, x1, y1);
    }
  }
  /////////////////////////////////////////////////////////////////////////////

  // Case: Diagnol at +/- 45 degrees
  /////////////////////////////////////////////////////////////////////////////
  else if (std::abs(dy) == std::abs(dx))
  {
    lineDiag(grid, x0, y0, x1, y1);
  }
  /////////////////////////////////////////////////////////////////////////////

  else
  {
    throw std::invalid_argument("Invalid Bresenham's Line Algorithm State");
  }
}

void OccupancyMapper::lineLow(GridMap& grid, int x0, int y0, int x1, int y1) const
{
  auto dx = x1 - x0;
  auto dy = y1 - y0;
  auto yi = 1;

  if (dy < 0)
  {
    yi = -1;
    dy = -dy;
  }

  auto D = 2 * dy - dx;
  auto y = y0;

  // counter
  auto ctr = 0;

  for (auto x = x0; x < x1; x++)
  {
    // Start as been added already in function call
    // prevents adding the end point in the case dx < 0
    if (ctr != 0)
    {
      updateCell(grid, x, y, log_odds_free_);
    }

    if (D > 0)
    {
      y += yi;
      D -= 2 * dx;
    }

    D += 2 * dy;
    ctr++;
  }  // end loop
}

void OccupancyMapper::lineHigh(GridMap& grid, int x0, int y0, int x1, int y1) const
{
  auto dx = x1 - x0;
  auto dy = y1 - y0;
  auto xi = 1;

  if (dx < 0)
  {
    xi = -1;
    dx = -dx;
  }

  auto D = 2 * dx - dy;
  auto x = x0;

  // counter
  auto ctr = 0;

  for (auto y = y0; y < y1; y++)
  {
    // Start as been added already in function call
    // prevents adding the end point in the case dy < 0
    if (ctr != 0)
    {
      updateCell(grid, x, y, log_odds_free_);
    }

    if (D > 0)
    {
      x += xi;
      D -= 2 * dy;
    }

    D += 2 * dx;
    ctr++;
  }  // end loop
}

void OccupancyMapper::lineDiag(GridMap& grid, int x0, int y0, int x1, int y1) const
{
  const auto dx = x1 - x0;
  const auto dy = y1 - y0;

  const auto xi = (dx < 0) ? -1 : 1;
  const auto yi = (dy < 0) ? -1 : 1;

  auto x = x0;
  auto y = y0;

  while (x != x1 and y != y1)
  {
    updateCell(grid, x, y, log_odds_free_);

    x += xi;
    y += yi;
  }
}

void OccupancyMapper::updateCell(GridMap& grid, int x, int y, double log) const
{
  const auto idx = grid.grid2RowMajor(y, x);
  const auto log_odds = prob2LogOdds(grid.getCell(idx)) + log - log_odds_prior_;

  // std::cout << logOdds2Prob(log_odds)*100.0 << std::endl;
  grid.grid_data_.at(idx) = static_cast<int8_t>(logOdds2Prob(log_odds) * 100.0);
}

}  // namespace ergodic_exploration
