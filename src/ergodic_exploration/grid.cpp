/**
 * @file grid.cpp
 * @author Boston Cleek
 * @date 23 Oct 2020
 * @brief 2D grid represented in row major order
 */

#include <iostream>

#include <ergodic_exploration/grid.hpp>

namespace ergodic_exploration
{
GridMap::GridMap(double xmin, double xmax, double ymin, double ymax, double resolution,
                 const std::shared_ptr<std::vector<int8_t>>& grid_ptr)
  : xmin_(xmin)
  , xmax_(xmax)
  , ymin_(ymin)
  , ymax_(ymax)
  , resolution_(resolution)
  , xsize_(axisLength(xmin, xmax, resolution))
  , ysize_(axisLength(ymin, ymax, resolution))
  , grid_ptr_(grid_ptr)
{
}

const std::shared_ptr<std::vector<int8_t>>& GridMap::getGridData() const
{
  return grid_ptr_;
}

unsigned int GridMap::grid2RowMajor(unsigned int i, unsigned int j) const
{
  // ith_row * Ncols + jth_col
  return i * xsize_ + j;
}

std::vector<unsigned int> GridMap::rowMajor2Grid(unsigned int idx) const
{
  const auto i = static_cast<unsigned int>(idx / xsize_);
  const auto j = idx - i * xsize_;
  return { i, j };
}

}  // namespace ergodic_exploration
