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
 * @file mapping.hpp
 * @author Boston Cleek
 * @date 8 Nov 2020
 * @brief Occupancy grid mapping
 */
#ifndef MAPPING_HPP
#define MAPPING_HPP

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
#endif
