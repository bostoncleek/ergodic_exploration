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
 * @file collision.hpp
 * @author Boston Cleek
 * @date 27 Oct 2020
 * @brief Collision checking in 2D
 */
#ifndef COLLISION_HPP
#define COLLISION_HPP

#include <utility>
#include <unordered_map>
#include <tuple>
#include <limits>
#include <armadillo>

#include <ergodic_exploration/grid.hpp>

namespace ergodic_exploration
{
using arma::vec;
using std::tuple;

/** @brief State of the collision detector */
enum class CollisionMsg
{
  crash,
  obstacle,
  none
};

/** @brief Collision detection parameters */
struct CollisionConfig
{
  /**
   * @brief Constructor
   * @param rb - collision boundary
   * @param rc - collision threshold
   * @param rm - max search distance
   * @param cx - x-coordinate (jth column) center of circle
   * @param cy - y-coordinate (ith row) center of circle
   */
  CollisionConfig(int rb, int rc, int rm, int cx, int cy)
    : r_bnd(rb), r_col(rc), r_max(rm), cx(cx), cy(cy), dx(0), dy(0), sqrd_obs(-1)
  {
  }

  const int r_bnd, r_col, r_max;  // bounding , collision, and max search radii
  const int cx, cy;               // circle center
  int dx, dy;                     // displacement from circle center to obstacle
  int sqrd_obs;                   // squared distance circle center to obstacle
};

/** @brief 2D collision detection */
class Collision
{
public:
  /**
   * @brief Constructor
   * @param boundary_radius - collision boundary around robot represented as a circle
   * @param search_radius - radius outside of collision boundary to search for obstacles
   * @param obstacle_threshold - distance from boundary radius to obstacle to be considered a collision
   * @param occupied_threshold - min probaility [0 1] for a cell to be considered an obstacle
   */
  Collision(double boundary_radius, double search_radius, double obstacle_threshold,
            double occupied_threshold);

  /**
   * @brief Compose distance to closest obstacle
   * @param grid - grid map
   * @param pose - robot state [x, y, theta]
   * @return state of collision detector and distance to closest obstacle
   * @details if collision distance to obstacle is zero and if there are no
   * obstacles distance is max numerical value for a double
   */
  tuple<CollisionMsg, double> minDistance(const GridMap& grid, const vec& pose) const;

  /**
   * @brief Compose displacement vector from closest obstacle to robot
   * @param grid - grid map
   * @param pose - robot state [x, y, theta]
   * @return state of collision detector and vector from closest obstacle
   * @details if collision distance to vector is all zeros and if there are no
   * obstacles the vector contains the max numerical value for a double
   */
  tuple<CollisionMsg, vec> minDirection(const GridMap& grid, const vec& pose) const;

  /**
   * @brief Check if the given state is a collision
   * @param grid - grid map
   * @param pose - robot state [x, y, theta]
   * @return true if collision
   */
  bool collisionCheck(const GridMap& grid, const vec& pose) const;

  /** @brief Return the distance from robot center to collision boundary */
  double totalPadding() const;

private:
  /**
   * @brief Find occupied cells between collision boundary and max search radius
   * @param cfg[out] - collision configuration
   * @param grid - grid map
   * @return true if there is a collision
   */
  bool search(CollisionConfig& cfg, const GridMap& grid) const;

  /**
   * @brief Bresenham's circle algorithm
   * @param cfg[out] - collision configuration
   * @param grid - grid map
   * @param r - radius of circle
   * @return true if there is a collision
   * @details Perfroms Bresenham's circle algorithm to detect obstacle cells
   */
  bool bresenhamCircle(CollisionConfig& cfg, const GridMap& grid, int r) const;

  /**
   * @brief Check if cell is an obstacle
   * @param cfg[out] - collision configuration
   * @param grid - grid map
   * @param cj - x-coordinate (jth column) of cell on perimeter of circle
   * @param ci - y-coordinate (ith row) of cell on perimeter of circle
   * @return true if there is a collision
   */
  bool checkCell(CollisionConfig& cfg, const GridMap& grid, unsigned int cj,
                 unsigned int ci) const;

private:
  double boundary_radius_;     // circular radius around robot
  double search_radius_;       // search radius for obstacles
  double obstacle_threshold_;  // collision distance threshold
  double occupied_threshold_;  // probaility cell is occupied
};
}  // namespace ergodic_exploration
#endif
