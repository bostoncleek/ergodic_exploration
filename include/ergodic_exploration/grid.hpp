/**
 * @file grid.hpp
 * @author Boston Cleek
 * @date 23 Oct 2020
 * @brief 2D grid represented in row major order
 */

#pragma once

#include <iosfwd>
#include <cmath>
#include <vector>
// #include <memory>

#include <geometry_msgs/Pose.h>
#include <nav_msgs/OccupancyGrid.h>

namespace ergodic_exploration
{
typedef std::vector<int8_t> gridData;

/**
 * @brief Length of a grid axis
 * @param lower - lower limit
 * @param upper - upper limit
 * @param resolution - resolution of grid
 * @return length
 */
constexpr unsigned int axisLength(double lower, double upper, double resolution)
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
constexpr double axisUpper(double lower, double resolution, unsigned int size)
{
  return static_cast<double>(resolution * size) + lower;
}

/** @brief Constructs an 2D grid */
class GridMap
{
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
          const gridData& grid_data);

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
   * @brief Get the grid data
   * @return grid data
   */
  const gridData& getGridData() const;

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
  int8_t getCell(double x, double y) const;

  /**
   * @brief Get the value of the cell at a specified index
   * @param i - row in the grid (y-coordinate)
   * @param j - column in the grid (x-coordinate)
   * @return value of the cell
   */
  int8_t getCell(unsigned int i, unsigned int j) const;

  /**
   * @brief Get the value of the cell at a specified index
   * @param idx - row major index
   * @return value of the cell
   */
  int8_t getCell(unsigned int idx) const;

private:
  unsigned int xsize_, ysize_;
  double resolution_, xmin_, ymin_, xmax_, ymax_;
  gridData grid_data_;
};

}  // namespace ergodic_exploration
