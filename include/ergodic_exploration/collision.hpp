/**
 * @file collision.hpp
 * @author Boston Cleek
 * @date 27 Oct 2020
 * @brief Collision checking in 2D
 */

#pragma once

#include <unordered_map>
#include <utility>

#include <ergodic_exploration/grid.hpp>

namespace ergodic_exploration
{
/** @brief Map obstacle cell index to grid coordinates */
typedef std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int>>
    CollisionMap;

/** @brief 2D collision detection */
class Collision
{
public:
  /**
   * @brief Constructor
   * @param boundary_radius - collision boundary around robot represented as a circle
   * @param search_radius - radius outside of collision boundary to search for obstacles
   * @param obstacle_threshold - min distance to obstacle to be considered a collision
   * @param occupied_threshold - min probaility [0 100] for a cell to be considered an obstacle
   */
  Collision(double boundary_radius, double search_radius, double obstacle_threshold,
            int8_t occupied_threshold);

  /**
   * @brief Get collision map
   * @return collision map
   */
  const CollisionMap& collisionMap() const;

  /**
   * @brief Find occupied cells in grid within the search radius
   * @param [in] <name> <parameter_description>
   * @return <return_description>
   * @details <details>
   */
  void obstacleCells(const GridMap& grid, unsigned int cx, unsigned int cy);

private:
  /**
   * @brief Bresenham's circle algorithm
   * @param grid - grid map
   * @param r - radius of circle
   * @param cx - x-coordinate (jth column) center of circle
   * @param cy - y-coordinate (ith row) center of circle
   * @details Perfroms Bresenham's circle algorithm to determine the obstalce cells
   */
  void bresenhamCircle(const GridMap& grid, int r, unsigned int cx, unsigned int cy);

  /**
   * @brief Adds cell to collision map
   * @param grid - grid map
   * @param i - row in the grid (y-coordinate)
   * @param j - column in the grid (x-coordinate)
   * @details Cell is added to collision map if within the boundary of the grid and meets
   * the occupied threshold
   */
  void addObstacleCell(const GridMap& grid, unsigned int i, unsigned int j);

  /**
   * @brief Compose the coordinates of cells on circle
   * @param grid - grid map
   * @param x - x-axis offset from center
   * @param y - y-axis offset from center
   * @param cx - x-coordinate (jth column) center of circle
   * @param cy - y-coordinate (ith row) center of circle
   */
  void cellCoordinates(const GridMap& grid, unsigned int x, unsigned int y,
                       unsigned int cx, unsigned int cy);

private:
  double boundary_radius_, search_radius_, obstacle_threshold_;
  int8_t occupied_threshold_;
  // TODO: pass this in as a parameter to obstacleCells() ??
  CollisionMap collision_map_;
};

}  // namespace ergodic_exploration
