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
  const auto xsize = ergodic_exploration::axis_length(xmin, xmax, resolution);
  const auto ysize = ergodic_exploration::axis_length(ymin, ymax, resolution);
  const ergodic_exploration::GridData grid_data(xsize * ysize, 0);

  const ergodic_exploration::GridMap grid_map(xmin, xmax, ymin, ymax, resolution,
                                              grid_data);

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
  const auto xsize = ergodic_exploration::axis_length(xmin, xmax, resolution);
  const auto ysize = ergodic_exploration::axis_length(ymin, ymax, resolution);
  const ergodic_exploration::GridData grid_data(xsize * ysize, 0);

  const ergodic_exploration::GridMap grid_map(xmin, xmax, ymin, ymax, resolution,
                                              grid_data);

  // Convert row major order index to row and column ndex
  const auto idx = 5;

  const auto indices = grid_map.rowMajor2Grid(idx);
  ASSERT_EQ(indices.at(0), 2);
  ASSERT_EQ(indices.at(1), 1);
}

TEST(GridTest, GridBoundsRowMajor)
{
  const auto xmin = 0.0;
  const auto xmax = 2.0;
  const auto ymin = 0.0;
  const auto ymax = 3.0;
  const auto resolution = 1.0;
  const auto xsize = ergodic_exploration::axis_length(xmin, xmax, resolution);
  const auto ysize = ergodic_exploration::axis_length(ymin, ymax, resolution);
  const ergodic_exploration::GridData grid_data(xsize * ysize, 0);

  const ergodic_exploration::GridMap grid_map(xmin, xmax, ymin, ymax, resolution,
                                              grid_data);

  const auto idx1 = 5;
  ASSERT_TRUE(grid_map.gridBounds(idx1));

  const auto idx2 = 6;
  ASSERT_FALSE(grid_map.gridBounds(idx2));
}

TEST(GridTest, GridBoundsCoordinates)
{
  const auto xmin = 0.0;
  const auto xmax = 2.0;
  const auto ymin = 0.0;
  const auto ymax = 3.0;
  const auto resolution = 1.0;
  const auto xsize = ergodic_exploration::axis_length(xmin, xmax, resolution);
  const auto ysize = ergodic_exploration::axis_length(ymin, ymax, resolution);
  const ergodic_exploration::GridData grid_data(xsize * ysize, 0);

  const ergodic_exploration::GridMap grid_map(xmin, xmax, ymin, ymax, resolution,
                                              grid_data);

  const auto i1 = 2;
  const auto j1 = 0;
  ASSERT_TRUE(grid_map.gridBounds(i1, j1));

  const auto i2 = 0;
  const auto j2 = 2;
  ASSERT_FALSE(grid_map.gridBounds(i2, j2));
}

TEST(GridTest, Grid2World)
{
  const auto xmin = -0.5;
  const auto xmax = 0.5;
  const auto ymin = 0.0;
  const auto ymax = 1.5;
  const auto resolution = 0.5;
  const auto xsize = ergodic_exploration::axis_length(xmin, xmax, resolution);
  const auto ysize = ergodic_exploration::axis_length(ymin, ymax, resolution);
  const ergodic_exploration::GridData grid_data(xsize * ysize, 0);

  const ergodic_exploration::GridMap grid_map(xmin, xmax, ymin, ymax, resolution,
                                              grid_data);

  const auto i1 = 0;
  const auto j1 = 1;

  const auto coords1 = grid_map.grid2World(i1, j1);
  ASSERT_DOUBLE_EQ(coords1.at(0), 0.25);
  ASSERT_DOUBLE_EQ(coords1.at(1), 0.25);

  const auto idx2 = 4;
  const auto coords2 = grid_map.grid2World(idx2);
  ASSERT_DOUBLE_EQ(coords2.at(0), -0.25);
  ASSERT_DOUBLE_EQ(coords2.at(1), 1.25);
}

TEST(GridTest, World2Grid)
{
  const auto xmin = -0.5;
  const auto xmax = 0.5;
  const auto ymin = 0.0;
  const auto ymax = 1.5;
  const auto resolution = 0.5;
  const auto xsize = ergodic_exploration::axis_length(xmin, xmax, resolution);
  const auto ysize = ergodic_exploration::axis_length(ymin, ymax, resolution);
  const ergodic_exploration::GridData grid_data(xsize * ysize, 0);

  const ergodic_exploration::GridMap grid_map(xmin, xmax, ymin, ymax, resolution,
                                              grid_data);

  const auto x1 = -0.25;
  const auto y1 = 1.25;

  const auto indices1 = grid_map.world2Grid(x1, y1);
  ASSERT_EQ(indices1.at(0), 2);
  ASSERT_EQ(indices1.at(1), 0);

  const auto x2 = 0.25;
  const auto y2 = 0.75;
  const auto idx = grid_map.world2RowMajor(x2, y2);
  ASSERT_EQ(idx, 3);
}

TEST(GridTest, GetCellValue)
{
  const auto xmin = -0.5;
  const auto xmax = 0.5;
  const auto ymin = 0.0;
  const auto ymax = 1.5;
  const auto resolution = 0.5;
  const auto xsize = ergodic_exploration::axis_length(xmin, xmax, resolution);
  const auto ysize = ergodic_exploration::axis_length(ymin, ymax, resolution);
  ergodic_exploration::GridData grid_data(xsize * ysize, 0);

  grid_data.at(3) = 100;
  grid_data.at(4) = 90;

  const ergodic_exploration::GridMap grid_map(xmin, xmax, ymin, ymax, resolution,
                                              grid_data);

  const auto x = -0.25;
  const auto y = 1.25;

  const unsigned int i = 1;
  const unsigned int j = 1;

  ASSERT_EQ(grid_map.getCell(4), 0.9);
  ASSERT_EQ(grid_map.getCell(x, y), 0.9);
  ASSERT_EQ(grid_map.getCell(i, j), 1.0);
}
