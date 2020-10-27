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
#include <memory>

namespace ergodic_exploration
{
/**
 * @brief Length of a grid axis
 * @param lower - lower limit
 * @param upper - upper limit
 * @param resolution - resolution of grid
 * @return lenght
 */
unsigned int axisLength(double lower, double upper, double resolution)
{
  return static_cast<unsigned int>(std::round((upper - lower) / resolution));
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
          const std::shared_ptr<std::vector<int8_t>>& grid_ptr);

  /**
   * @brief Get the grid data
   * @return grid data
   */
  const std::shared_ptr<std::vector<int8_t>>& getGridData() const;

  /**
   * @brief Converts grid indices to row major order index
   * @param i - row in the grid (y-coordinate)
   * @param j - column in the grid (x-coordinate)
   * @return grid index in row major order
   */
  unsigned int grid2RowMajor(unsigned int i, unsigned int j) const;

  /**
   * @brief Convert row major index to row and column indices
   * @param idx - row major column index
   * @return row and column indices
   */
  std::vector<unsigned int> rowMajor2Grid(unsigned int idx) const;

private:
  double xmin_, xmax_, ymin_, ymax_, resolution_;
  unsigned int xsize_, ysize_;

protected:
  std::shared_ptr<std::vector<int8_t>> grid_ptr_;
};

}  // namespace ergodic_exploration
