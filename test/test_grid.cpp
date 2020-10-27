/**
 * @file test_grid.cpp
 * @author Boston Cleek
 * @date 26 Oct 2020
 * @brief Test grid functionality
 */

#include <gtest/gtest.h>
#include <ergodic_exploration/grid.hpp>

TEST(GridTest, Grid2RowMajor)
{
  const auto xmin = 0.0;
  const auto xmax = 2.0;
  const auto ymin = 0.0;
  const auto ymax = 3.0;
  const auto resolution = 1.0;
  const std::shared_ptr<std::vector<int8_t>> grid_ptr;

  ergodic_exploration::GridMap grid_map(xmin, xmax, ymin, ymax, resolution, grid_ptr);

  // Convert grid index to row major order
  const auto i = 2;
  const auto j = 1;

  const auto idx = grid_map.grid2RowMajor(i, j);
  ASSERT_EQ(idx, 5);
}

TEST(GridTest, RowMajor2Grid)
{
  const auto xmin = 0.0;
  const auto xmax = 2.0;
  const auto ymin = 0.0;
  const auto ymax = 3.0;
  const auto resolution = 1.0;
  const std::shared_ptr<std::vector<int8_t>> grid_ptr;

  ergodic_exploration::GridMap grid_map(xmin, xmax, ymin, ymax, resolution, grid_ptr);

  // Convert row major order index to row and column ndex
  const auto idx = 5;

  const std::vector<unsigned int> indices = grid_map.rowMajor2Grid(idx);
  ASSERT_EQ(indices.at(0), 2);
  ASSERT_EQ(indices.at(1), 1);
}
