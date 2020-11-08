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
                     double obstacle_threshold, int8_t occupied_threshold)
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

  if (occupied_threshold_ > 100 || occupied_threshold_ < 0)
  {
    throw std::invalid_argument("Occupied threshold must be between 0 and 100");
  }
}

const CollisionMap& Collision::collisionMap() const
{
  return collision_map_;
}

void Collision::bresenhamCircle(const GridMap& grid, int r, unsigned int cx,
                                unsigned int cy)
{
  int x = 0;
  int y = r;
  int d = 3 - 2 * r;
  cellCoordinates(grid, x, y, cx, cy);
  while (y >= x)
  {
    x++;
    if (d > 0)
    {
      y--;
      d += 4 * (x - y) + 10;
    }
    else
    {
      d += 4 * x + 6;
    }
    cellCoordinates(grid, x, y, cx, cy);
  }

  // int x = -r;
  // int y = 0;
  // int err = 2 - 2 * r;
  //
  // while (x < 0)
  // {
  //   // Quadrant 1
  //   checkCell(grid, cy + y, cx - x);
  //
  //   // Quadrant 2
  //   checkCell(grid, cy - x, cx - y);
  //
  //   // Quadrant 3
  //   checkCell(grid, cy - y, cx + x);
  //
  //   // Quadrant 4
  //   checkCell(grid, cy + x, cx + y);
  //
  //   r = err;
  //
  //   if (r <= y)
  //   {
  //     y++;
  //     err += 2 * y + 1;
  //   }
  //
  //   if (r > x || err > y)
  //   {
  //     x++;
  //     err += 2 * x + 1;
  //   }
  // }
}

void Collision::obstacleCells(const GridMap& grid, unsigned int cx, unsigned int cy)
{
  collision_map_.clear();

  const int rmax = std::ceil(search_radius_ / grid.resolution());
  int r = std::ceil(boundary_radius_ / grid.resolution());

  while (r <= rmax)
  {
    bresenhamCircle(grid, r, cx, cy);
    r++;
  }
}

void Collision::addObstacleCell(const GridMap& grid, unsigned int i, unsigned int j)
{
  // Signed and unsigned ints were mixed
  // Check if (i,j) are within the bounds of the grid
  if (grid.gridBounds(i, j))
  {
    // Safe to convert (i,j) to row major index
    const auto idx = grid.grid2RowMajor(i, j);
    if (!(grid.getCell(idx) < occupied_threshold_) && !collision_map_.contains(idx))
    {
      collision_map_.emplace(idx, std::make_pair(i, j));
    }
  }
}

void Collision::cellCoordinates(const GridMap& grid, unsigned int x, unsigned int y,
                                unsigned int cx, unsigned int cy)
{
  addObstacleCell(grid, cy + y, cx + x);
  addObstacleCell(grid, cy + y, cx - x);
  addObstacleCell(grid, cy - y, cx + x);
  addObstacleCell(grid, cy - y, cx - x);

  addObstacleCell(grid, cy + x, cx + y);
  addObstacleCell(grid, cy + x, cx - y);
  addObstacleCell(grid, cy - x, cx + y);
  addObstacleCell(grid, cy - x, cx - y);
}

}  // namespace ergodic_exploration
