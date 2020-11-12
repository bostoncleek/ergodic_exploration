/**
 * @file mapping.hpp
 * @author Boston Cleek
 * @date 8 Nov 2020
 * @brief Occupancy grid mapping
 */

#pragma once

#include <armadillo>

#include <ergodic_exploration/grid.hpp>
#include <sensor_msgs/LaserScan.h>

namespace ergodic_exploration
{
using arma::mat;
using arma::vec;

/** @brief Occupancy grid mapping */
class OccupancyMapper
{
public:
  /**
   * @brief Constructor
   * @param Tbs - transformation from base link to frame of laser scanner
   */
  OccupancyMapper(const mat& Tbs);

  /**
   * @brief Update the occupancy grid map
   * @param grid[out] - current map
   * @param scan - laser scan
   * @param pose - robot pose in map frame
   * @return true if map is updated
   */
  bool updateMap(GridMap& grid, const sensor_msgs::LaserScan::ConstPtr& scan,
                 const vec& pose);

private:
  /**
   * @brief Bresenham's line algorithm on occupancy map
   * @param grid[out] - current map
   * @param x0 - line starting x position
   * @param y0 - line starting y position
   * @param x1 - line ending x position
   * @param y1 - line ending y position
   * @details Ray tracing is used to find the grid cells in free space.
   * Where x is the jth column and u is the ith row in the grid.
   */
  void rayTrace(GridMap& grid, int x0, int y0, int x1, int y1) const;

  /**
   * @brief Bresenham's algorithm for drawing lines downward
   * @param grid[out] - current map
   * @param x0 - line starting x position
   * @param y0 - line starting y position
   * @param x1 - line ending x position
   * @param y1 - line ending y position
   * @details x is the jth column and u is the ith row in the grid
   */
  void lineLow(GridMap& grid, int x0, int y0, int x1, int y1) const;

  /**
   * @brief Bresenham's algorithm for drawing lines upward
   * @param grid[out] - current map
   * @param x0 - line starting x position
   * @param y0 - line starting y position
   * @param x1 - line ending x position
   * @param y1 - line ending y position
   * @details x is the jth column and u is the ith row in the grid
   */
  void lineHigh(GridMap& grid, int x0, int y0, int x1, int y1) const;

  /**
   * @brief Bresenham's algorithm for drawing lines at +/- 45 deg
   * @param grid[out] - current map
   * @param x0 - line starting x position
   * @param y0 - line starting y position
   * @param x1 - line ending x position
   * @param y1 - line ending y position
   * @details x is the jth column and u is the ith row in the grid
   */
  void lineDiag(GridMap& grid, int x0, int y0, int x1, int y1) const;

  /**
   * @brief Update grid cell probability of occupancy
   * @param grid[out] - current map
   * @param x - grid cell x position
   * @param y - grid cell y position
   * @param log - either log_odds_occ_ or log_odds_free_
   * @details x is the jth column and u is the ith row in the grid
   */
  void updateCell(GridMap& grid, int x, int y, double log) const;

private:
  mat Tbs_;                              // transform base link to laser scanner
  double prior_, prob_occ_, prob_free_;  // probabilities for cell states
  double log_odds_prior_;                // log odds prior
  double log_odds_occ_;                  // log odds occupied
  double log_odds_free_;                 // log odds free
};

}  // namespace ergodic_exploration
