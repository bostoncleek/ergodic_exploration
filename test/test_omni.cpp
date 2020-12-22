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
