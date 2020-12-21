/**
 * @file test_omni.cpp
 * @author Boston Cleek
 * @date 30 Oct 2020
 * @brief Test omni robot models
 */

#include <gtest/gtest.h>
#include <ergodic_exploration/models/omni.hpp>

TEST(OmniTest, MecanumKinematics)
{
  const auto wheel_radius = 0.1;
  const auto wheel_base_x = 0.5;
  const auto wheel_base_y = 0.5;

  const ergodic_exploration::models::Mecanum mecanum(wheel_radius, wheel_base_x,
                                                     wheel_base_y);

  const arma::vec x = { 1.0, 2.0, 0.707 };
  const arma::vec u = { 0.5, 0.4, 0.6, 0.3 };
  const arma::vec xdot = mecanum(x, u);

  ASSERT_NEAR(xdot(0), 0.040709, 1e-6);
  ASSERT_NEAR(xdot(1), 0.021626, 1e-6);
  ASSERT_NEAR(xdot(2), 0.005, 1e-6);
}

TEST(OmniTest, MecanumJacobianState)
{
  const auto wheel_radius = 0.1;
  const auto wheel_base_x = 0.5;
  const auto wheel_base_y = 0.5;

  const ergodic_exploration::models::Mecanum mecanum(wheel_radius, wheel_base_x,
                                                     wheel_base_y);

  const arma::vec x = { 1.0, 2.0, 0.707 };
  const arma::vec u = { 0.5, 0.4, 0.6, 0.3 };
  const arma::mat A = mecanum.fdx(x, u);

  ASSERT_NEAR(A(0, 2), -0.021626, 1e-6);
  ASSERT_NEAR(A(1, 2), 0.040709, 1e-6);
}

TEST(OmniTest, MecanumJacobianControl)
{
  const auto wheel_radius = 0.1;
  const auto wheel_base_x = 0.5;
  const auto wheel_base_y = 0.5;

  const ergodic_exploration::models::Mecanum mecanum(wheel_radius, wheel_base_x,
                                                     wheel_base_y);

  const arma::vec x = { 1.0, 2.0, 0.707 };
  const arma::vec u = { 0.5, 0.4, 0.6, 0.3 };
  const arma::mat B = mecanum.fdu(x);

  ASSERT_NEAR(B(0, 0), 0.035246, 1e-6);
  ASSERT_NEAR(B(0, 1), 0.002768, 1e-6);
  ASSERT_NEAR(B(0, 2), 0.035246, 1e-6);
  ASSERT_NEAR(B(0, 3), 0.002768, 1e-6);

  ASSERT_NEAR(B(1, 0), -0.002768, 1e-6);
  ASSERT_NEAR(B(1, 1), 0.035246, 1e-6);
  ASSERT_NEAR(B(1, 2), -0.002768, 1e-6);
  ASSERT_NEAR(B(1, 3), 0.035246, 1e-6);

  ASSERT_NEAR(B(2, 0), -0.025, 1e-6);
  ASSERT_NEAR(B(2, 1), 0.025, 1e-6);
  ASSERT_NEAR(B(2, 2), 0.025, 1e-6);
  ASSERT_NEAR(B(2, 3), -0.025, 1e-6);
}
