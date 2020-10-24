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

/** @brief Grid cell */
struct Cell
{
  Cell(): i(-1), j(-1), marked(false) {}
  unsigned int i, j;
  bool marked;
};

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
   */
  GridMap(double xmin, double xmax, double ymin, double ymax, double resolution);

private:
  /**
   * @brief Converts grid indices to row major order index
   * @param i - row in the grid (y-coordinate)
   * @param j - column in the grid (x-coordinate)
   * @return grid index in row major order
   */
  unsigned int grid2RowMajor(int i, int j) const;

private:
  double xmin_, xmax_, ymin_, ymax_, resolution_;
  unsigned int xsize_, ysize_;

protected:
  std::vector<Cell> grid_;

};

}  // namespace ergodic_exploration
