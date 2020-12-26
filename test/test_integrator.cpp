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
 * @file test_integrator.cpp
 * @author Boston Cleek
 * @date 28 Oct 2020
 * @brief Test numerical integration
 */

#include <gtest/gtest.h>
#include <ergodic_exploration/models/cart.hpp>
#include <ergodic_exploration/models/omni.hpp>
#include <ergodic_exploration/integrator.hpp>

TEST(IntegratorTest, RungeKutta)
{
  const ergodic_exploration::models::Cart cart(0.1, 2.0);
  const auto horizon = 0.4;
  const auto dt = 0.1;
  const auto N = static_cast<unsigned int>(horizon / dt);
  arma::mat ut(2, N);
  ut.row(0).fill(1.0);
  ut.row(1).fill(1.0);

  const arma::vec x0 = { 0.0, 0.0, 0.0 };
  const ergodic_exploration::RungeKutta rk4(0.1);
  const arma::mat xt = rk4.solve(cart, x0, ut, horizon);

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

TEST(IntegratorTest, OmniPoseChange)
{
  const auto dt = 0.1;
  const arma::vec vb = { 1.0, 0.5, 0.5 };
  const arma::vec x0(3, arma::fill::zeros);

  const ergodic_exploration::models::Omni omni;
  const ergodic_exploration::RungeKutta rk4(dt);

  const arma::vec pose = ergodic_exploration::integrate_twist(x0, vb, dt);
  // pose.print("integrate twist");
  //
  // const arma::vec pose2 = rk4.step(omni, x0, vb);
  // pose2.print("rk4");

  ASSERT_NEAR(pose(0), 0.0987, 1e-4);
  ASSERT_NEAR(pose(1), 0.0525, 1e-4);
  ASSERT_NEAR(pose(2), 0.0500, 1e-4);
}

TEST(IntegratorTest, OmniTrajectory)
{
  const auto horizon = 0.5;
  const auto dt = 0.1;
  const auto steps = static_cast<unsigned int>(horizon / dt);
  const arma::vec vb = { 1.0, 0.5, 0.5 };
  const arma::vec x0 = { 1.0, 0.0, 0.0 };

  const ergodic_exploration::models::Omni omni;
  const ergodic_exploration::RungeKutta rk4(dt);

  arma::mat traj(3, steps, arma::fill::zeros);
  arma::vec x = x0;

  for (unsigned int i = 0; i < steps; i++)
  {
    x = ergodic_exploration::integrate_twist(x, vb, dt);
    traj.col(i) = x;
  }

  // arma::mat ut(3, steps);
  // ut.row(0).fill(vb(0));
  // ut.row(1).fill(vb(1));
  // ut.row(2).fill(vb(2));
  //
  // arma::mat traj2 = rk4.solve(omni, x0, ut, horizon);
  //
  // traj.print("integrate twist");
  // traj2.print("rk4");

  ASSERT_NEAR(traj(0, 0), 1.0987, 1e-4);
  ASSERT_NEAR(traj(0, 1), 1.1947, 1e-4);
  ASSERT_NEAR(traj(0, 2), 1.2876, 1e-4);
  ASSERT_NEAR(traj(0, 3), 1.3774, 1e-4);
  ASSERT_NEAR(traj(0, 4), 1.4637, 1e-4);

  ASSERT_NEAR(traj(1, 0), 0.0525, 1e-4);
  ASSERT_NEAR(traj(1, 1), 0.1098, 1e-4);
  ASSERT_NEAR(traj(1, 2), 0.1719, 1e-4);
  ASSERT_NEAR(traj(1, 3), 0.2385, 1e-4);
  ASSERT_NEAR(traj(1, 4), 0.3096, 1e-4);

  ASSERT_NEAR(traj(2, 0), 0.0500, 1e-4);
  ASSERT_NEAR(traj(2, 1), 0.1000, 1e-4);
  ASSERT_NEAR(traj(2, 2), 0.1500, 1e-4);
  ASSERT_NEAR(traj(2, 3), 0.2000, 1e-4);
  ASSERT_NEAR(traj(2, 4), 0.2500, 1e-4);
}
