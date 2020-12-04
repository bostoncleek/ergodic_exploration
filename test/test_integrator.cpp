/**
 * @file test_integrator.cpp
 * @author Boston Cleek
 * @date 28 Oct 2020
 * @brief Test numerical integration
 */

#include <gtest/gtest.h>
#include <ergodic_exploration/models/cart.hpp>
#include <ergodic_exploration/integrator.hpp>

TEST(RungeKuttaTest, Simulate)
{
  ergodic_exploration::models::Cart cart(0.1, 2.0);

  const auto horizon = 0.4;
  const auto dt = 0.1;
  const auto N = static_cast<unsigned int>(horizon / dt);
  arma::mat ut(2, N);
  ut.row(0).fill(1.0);
  ut.row(1).fill(1.0);

  arma::vec x0 = { 0.0, 0.0, 0.0 };
  ergodic_exploration::RungeKutta rk4(0.1);
  arma::mat xt;
  rk4.solve(xt, cart, x0, ut, horizon);

  // Expected to drive in straight line along the x-axis
  ASSERT_DOUBLE_EQ(xt(0, 0), 0.01);
  ASSERT_DOUBLE_EQ(xt(0, 1), 0.02);
  ASSERT_DOUBLE_EQ(xt(0, 2), 0.03);
  ASSERT_DOUBLE_EQ(xt(0, 3), 0.04);

  ASSERT_DOUBLE_EQ(xt(1, 0), 0.0);
  ASSERT_DOUBLE_EQ(xt(1, 1), 0.0);
  ASSERT_DOUBLE_EQ(xt(1, 2), 0.0);
  ASSERT_DOUBLE_EQ(xt(1, 3), 0.0);

  ASSERT_DOUBLE_EQ(xt(2, 0), 0.0);
  ASSERT_DOUBLE_EQ(xt(2, 1), 0.0);
  ASSERT_DOUBLE_EQ(xt(2, 2), 0.0);
  ASSERT_DOUBLE_EQ(xt(2, 3), 0.0);
}
