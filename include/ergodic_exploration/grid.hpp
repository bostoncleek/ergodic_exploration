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
 * @file grid.hpp
 * @author Boston Cleek
 * @date 23 Oct 2020
 * @brief 2D grid represented in row major order
 */
#ifndef GRID_HPP
#define GRID_HPP

#include <iosfwd>
#include <cmath>
#include <vector>
// #include <memory>

#include <geometry_msgs/Pose.h>
#include <nav_msgs/OccupancyGrid.h>

namespace ergodic_exploration
{
/** @brief Occupancy grid data */
typedef std::vector<int8_t> GridData;

/**
 * @brief Length of a grid axis
 * @param lower - lower limit
 * @param upper - upper limit
 * @param resolution - resolution of grid
 * @return length
 */
inline unsigned int axis_length(double lower, double upper, double resolution)
{
  return static_cast<unsigned int>(std::round((upper - lower) / resolution));
}

/**
 * @brief Upper axis limit
 * @param lower - position of lower axis limit
 * @param resolution - grid resolution
 * @param size - length of axis
 * @return axis upper limit
 */
inline double axis_upper(double lower, double resolution, unsigned int size)
{
  return static_cast<double>(resolution * size) + lower;
}

/** @brief Constructs an 2D grid */
class GridMap
{
private:
  friend class OccupancyMapper;

public:
  /**
   * @brief Constructor
   * @param xmin - x low bound of grid
   * @param xmax - x upper bound of grid
   * @param ymin - y low bound of grid
   * @param ymax - y upper bound of grid
   * @param resolution - grid resolution
   * @param grid_ptr - shared pointer to current grid data
   */
  GridMap(double xmin, double xmax, double ymin, double ymax, double resolution,
          const GridData& grid_data);

  /**
   * @brief Constructor
   * @param grid_msg - ROS occupancy grid message
   */
  GridMap(const nav_msgs::OccupancyGrid::ConstPtr& grid_msg);

  /**
   * @brief Constructor a empty grid
   */
  GridMap();

  /**
   * @brief Update the grid
   * @param grid_msg - ROS occupancy grid message
   */
  void update(const nav_msgs::OccupancyGrid::ConstPtr& grid_msg);

  /**
   * @brief Test if coordinates are within the grid
   * @param i - row in the grid (y-coordinate)
   * @param j - column in the grid (x-coordinate)
   * @return true if coordinates are within the grid
   */
  bool gridBounds(unsigned int i, unsigned int j) const;

  /**
   * @brief Test if coordinates are within the grid
   * @param idx - row major index
   * @return true if index is within the grid
   */
  bool gridBounds(unsigned int idx) const;

  /**
   * @brief Converts grid indices to row major order index
   * @param i - row in the grid (y-coordinate)
   * @param j - column in the grid (x-coordinate)
   * @return grid index in row major order
   */
  unsigned int grid2RowMajor(unsigned int i, unsigned int j) const;

  /**
   * @brief Convert row major index to row and column indices
   * @param idx - row major index
   * @return row and column indices
   */
  std::vector<unsigned int> rowMajor2Grid(unsigned int idx) const;

  /**
   * @brief Convert grid coordinates to world position
   * @param i - row in the grid (y-coordinate)
   * @param j - column in the grid (x-coordinate)
   * @return the x and y coordinates in the world
   */
  std::vector<double> grid2World(unsigned int i, unsigned int j) const;

  /**
   * @brief Convert row major index to world position
   * @param idx - row major index
   * @return the x and y coordinates in the world
   */
  std::vector<double> grid2World(unsigned int idx) const;

  /**
   * @brief Convert world coordinates to row and column index in grid
   * @param x - position in world (jth column)
   * @param y - position in world (ith row)
   * @return ith row and jth column in grid
   */
  std::vector<unsigned int> world2Grid(double x, double y) const;

  /**
   * @brief Convert world coordinates to row major index
   * @param x - position in world (jth column)
   * @param y - position in world (ith row)
   * @return row major index
   */
  unsigned int world2RowMajor(double x, double y) const;

  /**
   * @brief Get the value of the cell at a specified position
   * @param x - position in world (jth column)
   * @param y - position in world (ith row)
   * @return value of the cell
   */
  double getCell(double x, double y) const;

  /**
   * @brief Get the value of the cell at a specified index
   * @param i - row in the grid (y-coordinate)
   * @param j - column in the grid (x-coordinate)
   * @return value of the cell
   */
  double getCell(unsigned int i, unsigned int j) const;

  /**
   * @brief Get the value of the cell at a specified index
   * @param idx - row major index
   * @return value of the cell
   */
  double getCell(unsigned int idx) const;

  /** @brief Print grid properties */
  void print() const;

  /** @brief Return the grid data  */
  const GridData& gridData() const
  {
    return grid_data_;
  }

  /** @brief Return grid resolution */
  double resolution() const
  {
    return resolution_;
  }

  /** @brief Return grid xmin */
  double xmin() const
  {
    return xmin_;
  }

  /** @brief Return grid ymin */
  double ymin() const
  {
    return ymin_;
  }

  /** @brief Return grid xmax */
  double xmax() const
  {
    return xmax_;
  }

  /** @brief Return grid ymax */
  double ymax() const
  {
    return ymax_;
  }

  /** @brief Return grid x-axis size */
  unsigned int xsize() const
  {
    return xsize_;
  }

  /** @brief Return y-axis size */
  unsigned int ysize() const
  {
    return ysize_;
  }

  /** @brief Return grid size */
  unsigned int size() const
  {
    return xsize_ * ysize_;
  }

private:
  unsigned int xsize_, ysize_;
  double resolution_, xmin_, ymin_, xmax_, ymax_;
  GridData grid_data_;
};

}  // namespace ergodic_exploration
#endif
