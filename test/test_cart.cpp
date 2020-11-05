/**
 * @file test_cart.cpp
 * @author Boston Cleek
 * @date 30 Oct 2020
 * @brief Test cart model
 */

#include <gtest/gtest.h>
#include <ergodic_exploration/cart.hpp>
#include <ergodic_exploration/types.hpp>

TEST(CartTest, CartKinematics)
{
  const auto wheel_radius = 0.033;
  const auto wheel_base = 0.08;
  ergodic_exploration::Cart cart(wheel_radius, wheel_base);

  const arma::vec x = { 1.0, 2.0, 0.707 };
  const arma::vec u = { 0.1, 0.2 };
  const arma::vec xdot = cart(x, u);

  ASSERT_NEAR(xdot(0), 0.003763, 1e-6);
  ASSERT_NEAR(xdot(1), 0.003215, 1e-6);
  ASSERT_NEAR(xdot(2), 0.020625, 1e-6);
}

TEST(CartTest, CartJacobianState)
{
  const auto wheel_radius = 0.033;
  const auto wheel_base = 0.08;
  ergodic_exploration::Cart cart(wheel_radius, wheel_base);

  const arma::vec x = { 1.0, 2.0, 0.707 };
  const arma::vec u = { 0.1, 0.2 };
  const arma::mat A = cart.fdx(x, u);

  ASSERT_NEAR(A(0, 2), -0.003215, 1e-6);
  ASSERT_NEAR(A(1, 2), 0.003763, 1e-6);
}

TEST(CartTest, CartJacobianControl)
{
  const auto wheel_radius = 0.033;
  const auto wheel_base = 0.08;
  ergodic_exploration::Cart cart(wheel_radius, wheel_base);

  const arma::vec x = { 1.0, 2.0, 0.707 };
  const arma::mat B = cart.fdu(x);

  ASSERT_NEAR(B(0, 0), 0.012545, 1e-6);
  ASSERT_NEAR(B(0, 1), 0.012545, 1e-6);

  ASSERT_NEAR(B(1, 0), 0.010717, 1e-6);
  ASSERT_NEAR(B(1, 1), 0.010717, 1e-6);

  ASSERT_NEAR(B(2, 0), -0.20625, 1e-6);
  ASSERT_NEAR(B(2, 1), 0.20625, 1e-6);
}

TEST(CartTest, CartWheels2TwistStraight)
{
  const auto wheel_radius = 0.033;
  const auto wheel_base = 0.08;
  ergodic_exploration::Cart cart(wheel_radius, wheel_base);

  const arma::vec u = { 1.0, 1.0 };
  ergodic_exploration::Twist2D vb = cart.wheels2Twist(u);

  ASSERT_NEAR(vb.vx, 0.033, 1e-6);
  ASSERT_NEAR(vb.vy, 0.0, 1e-6);
  ASSERT_NEAR(vb.w, 0.0, 1e-6);
}

TEST(CartTest, CartWheels2LeftTurn)
{
  const auto wheel_radius = 0.033;
  const auto wheel_base = 0.08;
  ergodic_exploration::Cart cart(wheel_radius, wheel_base);

  const arma::vec u = { -1.0, 1.0 };
  ergodic_exploration::Twist2D vb = cart.wheels2Twist(u);

  ASSERT_NEAR(vb.vx, 0.0, 1e-6);
  ASSERT_NEAR(vb.vy, 0.0, 1e-6);
  ASSERT_NEAR(vb.w, 0.4125, 1e-6);
}

TEST(CartTest, CartWheels2RightTurn)
{
  const auto wheel_radius = 0.033;
  const auto wheel_base = 0.08;
  ergodic_exploration::Cart cart(wheel_radius, wheel_base);

  const arma::vec u = { 1.0, -1.0 };
  ergodic_exploration::Twist2D vb = cart.wheels2Twist(u);

  ASSERT_NEAR(vb.vx, 0.0, 1e-6);
  ASSERT_NEAR(vb.vy, 0.0, 1e-6);
  ASSERT_NEAR(vb.w, -0.4125, 1e-6);
}

TEST(CartTest, SimpleCartKinematics)
{
  ergodic_exploration::SimpleCart cart;

  const arma::vec x = { 1.0, 2.0, 0.707 };
  const arma::vec u = { 0.5, 0.01 };
  const arma::vec xdot = cart(x, u);

  ASSERT_NEAR(xdot(0), 0.380156, 1e-6);
  ASSERT_NEAR(xdot(1), 0.324777, 1e-6);
  ASSERT_NEAR(xdot(2), 0.01, 1e-6);
}

TEST(CartTest, SimpleCartJacobianState)
{
  ergodic_exploration::SimpleCart cart;

  const arma::vec x = { 1.0, 2.0, 0.707 };
  const arma::vec u = { 0.5, 0.01 };
  const arma::mat A = cart.fdx(x, u);

  ASSERT_NEAR(A(0, 2), -0.324777, 1e-6);
  ASSERT_NEAR(A(1, 2), 0.380156, 1e-6);
}

TEST(CartTest, SimpleCartJacobianControl)
{
  ergodic_exploration::SimpleCart cart;

  const arma::vec x = { 1.0, 2.0, 0.707 };
  const arma::vec u = { 0.5, 0.01 };
  const arma::mat B = cart.fdu(x);

  ASSERT_NEAR(B(0, 0), 0.760313, 1e-6);
  ASSERT_NEAR(B(1, 0), 0.649555, 1e-6);
}
