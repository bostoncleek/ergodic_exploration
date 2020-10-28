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
    collisionMap;

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
  const collisionMap& getCollisionMap() const;

  // TODO: move to private after testing
  /**
   * @brief Bresenham's circle algorithm
   * @param [in] <name> <parameter_description>
   * @return <return_description>
   * @details Perfroms Bresenham's circle algorithm to determine the obstalce cells
   */
  void bresenhamCircle(const GridMap& grid, int r, unsigned int cx, unsigned int cy);

  /**
   * @brief Find occupied cells in grid within the search radius
   * @param [in] <name> <parameter_description>
   * @return <return_description>
   * @details <details>
   */
  void obstacleCells(const GridMap& grid, unsigned int cx, unsigned int cy);

private:
  /**
  * @brief Adds cell to collision map
  * @param [in] <name> <parameter_description>
  * @return <return_description>
  * @details Cell is added if within the boundary of the grid and meets the occupied threshold
  */
  void checkCell(const GridMap& grid, unsigned int i, unsigned int j);

  void circleCoordinates(const GridMap& grid, int x, int y, unsigned int cx, unsigned int cy);



private:
  double boundary_radius_, search_radius_, obstacle_threshold_;
  int8_t occupied_threshold_;
  // TODO: pass this in as a parameter to obstacleCells() ??
  collisionMap collision_map_;
};

}  // namespace ergodic_exploration
