/*********************************************************************
 * BSD 3-Clause License
 *
 * Copyright (c) 2020 Northwestern University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 *  * Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/
/**
 * @file collision.cpp
 * @author Boston Cleek
 * @date 27 Oct 2020
 * @brief Collision checking in 2D
 */

#include <iostream>
#include <stdexcept>

#include <ergodic_exploration/collision.hpp>

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

tuple<CollisionMsg, double> Collision::minDistance(const GridMap& grid,
                                                   const vec& pose) const
{
  const auto psg = grid.world2Grid(pose(0), pose(1));

  CollisionConfig cfg(
      std::floor(boundary_radius_ / grid.resolution()),
      std::floor((boundary_radius_ + obstacle_threshold_) / grid.resolution()),
      std::floor(search_radius_ / grid.resolution()), psg.at(1), psg.at(0));

  // std::cout << "r_bnd: " << cfg.r_bnd << std::endl;
  // std::cout << "r_col: " << cfg.r_col << std::endl;
  // std::cout << "r_max: " << cfg.r_max << std::endl;

  if (search(cfg, grid))
  {
    return std::make_tuple(CollisionMsg::crash, 0.0);
  }

  else if (cfg.sqrd_obs != -1)
  {
    const double dmin = std::sqrt(static_cast<double>(cfg.sqrd_obs)) * grid.resolution() -
                        boundary_radius_;

    return std::make_tuple(CollisionMsg::obstacle, dmin);
  }

  return std::make_tuple(CollisionMsg::none, std::numeric_limits<double>::max());
}

tuple<CollisionMsg, vec> Collision::minDirection(const GridMap& grid,
                                                 const vec& pose) const
{
  const auto psg = grid.world2Grid(pose(0), pose(1));

  CollisionConfig cfg(
      std::floor(boundary_radius_ / grid.resolution()),
      std::floor((boundary_radius_ + obstacle_threshold_) / grid.resolution()),
      std::floor(search_radius_ / grid.resolution()), psg.at(1), psg.at(0));

  search(cfg, grid);

  if (search(cfg, grid))
  {
    const vec disp(2, arma::fill::zeros);
    return std::make_tuple(CollisionMsg::crash, disp);
  }

  if (cfg.sqrd_obs != -1)
  {
    // displacement robot center to obstacle
    const vec disp = { grid.resolution() * static_cast<double>(cfg.dx),
                       grid.resolution() * static_cast<double>(cfg.dy) };
    return std::make_tuple(CollisionMsg::obstacle, disp);
  }

  const vec disp = { std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max() };
  return std::make_tuple(CollisionMsg::none, disp);
}

bool Collision::collisionCheck(const GridMap& grid, const vec& pose) const
{
  const auto psg = grid.world2Grid(pose(0), pose(1));

  CollisionConfig cfg(
      std::floor(boundary_radius_ / grid.resolution()),
      std::floor((boundary_radius_ + obstacle_threshold_) / grid.resolution()),
      std::floor(search_radius_ / grid.resolution()), psg.at(1), psg.at(0));

  if (search(cfg, grid))
  {
    return true;
  }

  return false;
}

double Collision::totalPadding() const
{
  return boundary_radius_ + obstacle_threshold_;
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
    const auto sqrd_obs =
        static_cast<int>((cfg.cx - cj) * (cfg.cx - cj) + (cfg.cy - ci) * (cfg.cy - ci));

    // update smallest squared distance to obstacle
    if (sqrd_obs < cfg.sqrd_obs || cfg.sqrd_obs == -1)
    {
      cfg.sqrd_obs = sqrd_obs;
      cfg.dx = cj - cfg.cx;
      cfg.dy = ci - cfg.cy;
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
