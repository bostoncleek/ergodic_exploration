/**
 * @file grid.cpp
 * @author Boston Cleek
 * @date 23 Oct 2020
 * @brief 2D grid represented in row major order
 */

#include <iostream>
#include <stdexcept>

#include <ergodic_exploration/grid.hpp>

namespace ergodic_exploration
{
GridMap::GridMap(double xmin, double xmax, double ymin, double ymax, double resolution,
                 const gridData& grid_data)
  : xsize_(axisLength(xmin, xmax, resolution))
  , ysize_(axisLength(ymin, ymax, resolution))
  , resolution_(resolution)
  , xmin_(xmin)
  , ymin_(ymin)
  , xmax_(xmax)
  , ymax_(ymax)
  , grid_data_(grid_data)
{
}

GridMap::GridMap(const nav_msgs::OccupancyGrid::ConstPtr& grid_msg)
  : xsize_(grid_msg->info.width)
  , ysize_(grid_msg->info.height)
  , resolution_(grid_msg->info.resolution)
  , xmin_(grid_msg->info.origin.position.x)
  , ymin_(grid_msg->info.origin.position.y)
  , xmax_(axisUpper(xmin_, resolution_, xsize_))
  , ymax_(axisUpper(ymin_, resolution_, ysize_))
  , grid_data_(grid_msg->data)
{
}

GridMap::GridMap()
  : xsize_(0), ysize_(0), resolution_(0), xmin_(0.0), ymin_(0.0), xmax_(0.0), ymax_(0.0)
{
}

const gridData& GridMap::getGridData() const
{
  return grid_data_;
}

bool GridMap::gridBounds(unsigned int i, unsigned int j) const
{
  // unsigned int is always greater than 0
  return ((i <= ysize_ - 1) and (j <= xsize_ - 1)) ? true : false;
}

bool GridMap::gridBounds(unsigned int idx) const
{
  const auto indices = rowMajor2Grid(idx);
  return ((indices.at(0) <= ysize_ - 1) and (indices.at(1) <= xsize_ - 1)) ? true : false;
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

std::vector<double> GridMap::grid2World(unsigned int i, unsigned int j) const
{
  const auto x = static_cast<double>(j * resolution_) + resolution_ / 2.0 + xmin_;
  const auto y = static_cast<double>(i * resolution_) + resolution_ / 2.0 + ymin_;
  return { x, y };
}

std::vector<double> GridMap::grid2World(unsigned int idx) const
{
  const auto indices = rowMajor2Grid(idx);
  const auto x =
      static_cast<double>(indices.at(1) * resolution_) + resolution_ / 2.0 + xmin_;
  const auto y =
      static_cast<double>(indices.at(0) * resolution_) + resolution_ / 2.0 + ymin_;
  return { x, y };
}

std::vector<unsigned int> GridMap::world2Grid(double x, double y) const
{
  unsigned int j = static_cast<unsigned int>(std::floor((x - xmin_) / resolution_));
  unsigned int i = static_cast<unsigned int>(std::floor((y - ymin_) / resolution_));

  if (j == xsize_)
  {
    j--;
  }

  if (i == ysize_)
  {
    i--;
  }

  return { i, j };
}

unsigned int GridMap::world2RowMajor(double x, double y) const
{
  const auto indices = world2Grid(x, y);
  return grid2RowMajor(indices.at(0), indices.at(1));
}

int8_t GridMap::getCell(double x, double y) const
{
  return getCell(world2RowMajor(x, y));
}

int8_t GridMap::getCell(unsigned int i, unsigned int j) const
{
  return getCell(grid2RowMajor(i, j));
}

int8_t GridMap::getCell(unsigned int idx) const
{
  if (!gridBounds(idx))
  {
    throw std::invalid_argument("Grid index out of range");
  }
  return grid_data_.at(idx);
}

}  // namespace ergodic_exploration
