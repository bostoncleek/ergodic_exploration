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
 * @file test_cart.cpp
 * @author Boston Cleek
 * @date 30 Oct 2020
 * @brief Test cart model
 */

#include <gtest/gtest.h>
#include <ergodic_exploration/models/cart.hpp>

TEST(CartTest, CartKinematics)
{
  const auto wheel_radius = 0.033;
  const auto wheel_base = 0.08;
  const ergodic_exploration::models::Cart cart(wheel_radius, wheel_base);

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
  const ergodic_exploration::models::Cart cart(wheel_radius, wheel_base);

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
  const ergodic_exploration::models::Cart cart(wheel_radius, wheel_base);

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
  const ergodic_exploration::models::Cart cart(wheel_radius, wheel_base);

  const arma::vec u = { 1.0, 1.0 };
  const arma::vec vb = cart.wheels2Twist(u);

  ASSERT_NEAR(vb(0), 0.033, 1e-6);
  ASSERT_NEAR(vb(1), 0.0, 1e-6);
  ASSERT_NEAR(vb(2), 0.0, 1e-6);
}

TEST(CartTest, CartWheels2LeftTurn)
{
  const auto wheel_radius = 0.033;
  const auto wheel_base = 0.08;
  const ergodic_exploration::models::Cart cart(wheel_radius, wheel_base);

  const arma::vec u = { -1.0, 1.0 };
  const arma::vec vb = cart.wheels2Twist(u);

  ASSERT_NEAR(vb(0), 0.0, 1e-6);
  ASSERT_NEAR(vb(1), 0.0, 1e-6);
  ASSERT_NEAR(vb(2), 0.4125, 1e-6);
}

TEST(CartTest, CartWheels2RightTurn)
{
  const auto wheel_radius = 0.033;
  const auto wheel_base = 0.08;
  const ergodic_exploration::models::Cart cart(wheel_radius, wheel_base);

  const arma::vec u = { 1.0, -1.0 };
  const arma::vec vb = cart.wheels2Twist(u);

  ASSERT_NEAR(vb(0), 0.0, 1e-6);
  ASSERT_NEAR(vb(1), 0.0, 1e-6);
  ASSERT_NEAR(vb(2), -0.4125, 1e-6);
}

TEST(CartTest, SimpleCartKinematics)
{
  const ergodic_exploration::models::SimpleCart cart;

  const arma::vec x = { 1.0, 2.0, 0.707 };
  const arma::vec u = { 0.5, 0.0, 0.01 };
  const arma::vec xdot = cart(x, u);

  ASSERT_NEAR(xdot(0), 0.380156, 1e-6);
  ASSERT_NEAR(xdot(1), 0.324777, 1e-6);
  ASSERT_NEAR(xdot(2), 0.01, 1e-6);
}

TEST(CartTest, SimpleCartJacobianState)
{
  const ergodic_exploration::models::SimpleCart cart;

  const arma::vec x = { 1.0, 2.0, 0.707 };
  const arma::vec u = { 0.5, 0.0, 0.01 };
  const arma::mat A = cart.fdx(x, u);

  ASSERT_NEAR(A(0, 2), -0.324777, 1e-6);
  ASSERT_NEAR(A(1, 2), 0.380156, 1e-6);
}

TEST(CartTest, SimpleCartJacobianControl)
{
  const ergodic_exploration::models::SimpleCart cart;

  const arma::vec x = { 1.0, 2.0, 0.707 };
  const arma::vec u = { 0.5, 0.01 };
  const arma::mat B = cart.fdu(x);

  ASSERT_NEAR(B(0, 0), 0.760313, 1e-6);
  ASSERT_NEAR(B(1, 0), 0.649555, 1e-6);
}
