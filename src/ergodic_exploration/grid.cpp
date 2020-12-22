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
                 const GridData& grid_data)
  : xsize_(axis_length(xmin, xmax, resolution))
  , ysize_(axis_length(ymin, ymax, resolution))
  , resolution_(resolution)
  , xmin_(xmin)
  , ymin_(ymin)
  , xmax_(xmax)
  , ymax_(ymax)
  , grid_data_(grid_data)
{
  if (xsize_ * ysize_ != grid_data_.size())
  {
    throw std::invalid_argument("Grid data size does not match the grid size");
  }
}

GridMap::GridMap(const nav_msgs::OccupancyGrid::ConstPtr& grid_msg)
  : xsize_(grid_msg->info.width)
  , ysize_(grid_msg->info.height)
  , resolution_(grid_msg->info.resolution)
  , xmin_(grid_msg->info.origin.position.x)
  , ymin_(grid_msg->info.origin.position.y)
  , xmax_(axis_upper(xmin_, resolution_, xsize_))
  , ymax_(axis_upper(ymin_, resolution_, ysize_))
  , grid_data_(grid_msg->data)
{
  if (xsize_ * ysize_ != grid_data_.size())
  {
    throw std::invalid_argument("Grid data size does not match the grid size");
  }
}

GridMap::GridMap()
  : xsize_(0), ysize_(0), resolution_(0), xmin_(0.0), ymin_(0.0), xmax_(0.0), ymax_(0.0)
{
}

void GridMap::update(const nav_msgs::OccupancyGrid::ConstPtr& grid_msg)
{
  xsize_ = grid_msg->info.width;
  ysize_ = grid_msg->info.height;
  resolution_ = grid_msg->info.resolution;
  xmin_ = grid_msg->info.origin.position.x;
  ymin_ = grid_msg->info.origin.position.y;
  xmax_ = axis_upper(xmin_, resolution_, xsize_);
  ymax_ = axis_upper(ymin_, resolution_, ysize_);
  grid_data_ = grid_msg->data;
}

bool GridMap::gridBounds(unsigned int i, unsigned int j) const
{
  // unsigned int is always greater than 0
  return ((i <= ysize_ - 1) and (j <= xsize_ - 1)) ? true : false;
}

bool GridMap::gridBounds(unsigned int idx) const
{
  return (idx <= (xsize_ * ysize_ - 1)) ? true : false;
  // std::vector<unsigned int> indices = rowMajor2Grid(idx);
  // return gridBounds(indices.at(0), indices.at(1));
}

unsigned int GridMap::grid2RowMajor(unsigned int i, unsigned int j) const
{
  if (!gridBounds(i, j))
  {
    std::cout << "WARNING (grid2RowMajor) i and j NOT within bounds" << std::endl;
  }
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

double GridMap::getCell(double x, double y) const
{
  return getCell(world2RowMajor(x, y));
}

double GridMap::getCell(unsigned int i, unsigned int j) const
{
  return getCell(grid2RowMajor(i, j));
}

double GridMap::getCell(unsigned int idx) const
{
  if (!gridBounds(idx))
  {
    throw std::invalid_argument("Grid index out of range");
  }
  return static_cast<double>(grid_data_.at(idx)) / 100.0;
}

void GridMap::print() const
{
  std::cout << "Grid \n"
            << "x-axis: [" << xmin_ << ", " << xmax_ << "] \n"
            << "y-axis: [" << ymin_ << ", " << ymax_ << "] \n"
            << "xsize: " << xsize_ << "\n"
            << "ysize: " << ysize_ << "\n"
            << "size: " << grid_data_.size() << std::endl;
}

}  // namespace ergodic_exploration
