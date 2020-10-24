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
GridMap::GridMap(double xmin, double xmax, double ymin, double ymax, double resolution)
  : xmin_(xmin)
  , xmax_(xmax)
  , ymin_(ymin)
  , ymax_(ymax)
  , resolution_(resolution)
  , xsize_(axisLength(xmin, xmax, resolution))
  , ysize_(axisLength(ymin, ymax, resolution))
  , grid_(xsize_ * ysize_)
{
}

unsigned int GridMap::grid2RowMajor(int i, int j) const
{
  // ith_row * Ncols + jth_col
  return i * xsize_ + j;
}

}  // namespace ergodic_exploration
