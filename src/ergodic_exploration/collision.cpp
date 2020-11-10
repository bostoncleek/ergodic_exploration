/**
 * @file collision.cpp
 * @author Boston Cleek
 * @date 27 Oct 2020
 * @brief Collision checking in 2D
 */

#include <iostream>
#include <stdexcept>

#include <ergodic_exploration/collision.hpp>
#include <ergodic_exploration/numerics.hpp>

namespace ergodic_exploration
{
Collision::Collision(double boundary_radius, double search_radius,
                     double obstacle_threshold, double occupied_threshold)
  : boundary_radius_(boundary_radius)
  , search_radius_(search_radius)
  , obstacle_threshold_(obstacle_threshold)
  , occupied_threshold_(occupied_threshold)
{
  if (search_radius_ < boundary_radius_)
  {
    throw std::invalid_argument(
        "Search radius must be at least the same size as the boundary radius");
  }

  if (occupied_threshold_ > 100.0 || occupied_threshold_ < 0.0)
  {
    throw std::invalid_argument("Occupied threshold must be between 0 and 100");
  }
}

CollisionMsg Collision::minDistance(double& dmin, const GridMap& grid,
                                    const vec& pose) const
{
  // pose in grid => {i,j}
  const auto psg = grid.world2Grid(pose(0), pose(1));

  CollisionConfig cfg(
      std::ceil(boundary_radius_ / grid.resolution()),
      std::ceil((boundary_radius_ + obstacle_threshold_) / grid.resolution()),
      std::ceil(search_radius_ / grid.resolution()), psg.at(1), psg.at(0));

  // std::cout << "r_bnd: " << cfg.r_bnd << std::endl;
  // std::cout << "r_col: " << cfg.r_col << std::endl;
  // std::cout << "r_max: " << cfg.r_max << std::endl;

  if (search(cfg, grid))
  {
    return CollisionMsg::crash;
  }

  else if (cfg.sqrd_obs != -1)
  {
    dmin = std::sqrt(static_cast<double>(cfg.sqrd_obs)) * grid.resolution() -
           boundary_radius_;

    return CollisionMsg::obstacle;
  }

  return CollisionMsg::none;
}

CollisionMsg Collision::minDirection(vec& disp, const GridMap& grid, const vec& pose) const
{
  // TODO: add parameter to enable collision detection
  // otherwise it will return the direction at the collision

  // TODO: consider weights in barrier function based on collisions

  // pose in grid => {i,j}
  const auto psg = grid.world2Grid(pose(0), pose(1));
  // std::cout << "cx: " << psg.at(1) << " cy: " << psg.at(0) << std::endl;

  CollisionConfig cfg(
      std::ceil(boundary_radius_ / grid.resolution()),
      std::ceil((boundary_radius_ + obstacle_threshold_) / grid.resolution()),
      std::ceil(search_radius_ / grid.resolution()), psg.at(1), psg.at(0));

  search(cfg, grid);

  // robot center to obstacle
  const auto d = std::sqrt(static_cast<double>(cfg.sqrd_obs)) * grid.resolution();
  // ratio of point on boundary to obstacle over robot center to obstacle
  const auto ratio = (d - boundary_radius_) / d;

  if (cfg.sqrd_obs != -1)
  {
    disp.resize(2);
    disp(0) = ratio * static_cast<double>(cfg.dx);
    disp(1) = ratio * static_cast<double>(cfg.dy);

    return CollisionMsg::obstacle;
  }

  return CollisionMsg::none;
}

bool Collision::search(CollisionConfig& cfg, const GridMap& grid) const
{
  // initialize start and stop radius
  const auto rmax = cfg.r_max;
  auto r = cfg.r_bnd;

  while (r <= rmax)
  {
    if (bresenhamCircle(cfg, grid, r))
    {
      return true;
    }
    r++;
  }

  return false;
}

bool Collision::bresenhamCircle(CollisionConfig& cfg, const GridMap& grid, int r) const
{
  int x = -r;
  int y = 0;
  int err = 2 - 2 * r;

  while (x < 0)
  {
    // Quadrant 1
    if (checkCell(cfg, grid, cfg.cx - x, cfg.cy + y))
    {
      return true;
    }

    // Quadrant 2
    if (checkCell(cfg, grid, cfg.cx - y, cfg.cy - x))
    {
      return true;
    }

    // Quadrant 3
    if (checkCell(cfg, grid, cfg.cx + x, cfg.cy - y))
    {
      return true;
    }

    // Quadrant 4
    if (checkCell(cfg, grid, cfg.cx + y, cfg.cy + x))
    {
      return true;
    }

    r = err;

    if (r <= y)
    {
      y++;
      err += 2 * y + 1;
    }

    if (r > x || err > y)
    {
      x++;
      err += 2 * x + 1;
    }
  }

  return false;
}

bool Collision::checkCell(CollisionConfig& cfg, const GridMap& grid, unsigned int cj,
                          unsigned int ci) const
{
  // within the bounds of the grid
  if (grid.gridBounds(ci, cj) && !(grid.getCell(ci, cj) < occupied_threshold_))
  {
    // squared distance circle center to obstacle
    const auto sqrd_obs = intDistSqaured(cfg.cx, cfg.cy, cj, ci);

    // update smallest squared distance to obstacle
    if (sqrd_obs < cfg.sqrd_obs || cfg.sqrd_obs == -1)
    {
      // std::cout << "ci: " << ci << " cj: " << cj << " p: " << grid.getCell(ci, cj)
      //           << std::endl;
      //
      // std::cout << "squared distance: " << sqrd_obs << std::endl;

      cfg.sqrd_obs = sqrd_obs;
      cfg.dx = cfg.cx - cj;
      cfg.dy = cfg.cy - ci;
    }

    // squared distance is >= 0 the inequality holds
    // check for collision
    if (sqrd_obs <= cfg.r_col * cfg.r_col)
    {
      return true;
    }
  }

  return false;
}
}  // namespace ergodic_exploration
