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
 * @file dynamic_window.cpp
 * @author Boston Cleek
 * @date 15 Nov 2020
 * @brief Dynamic window control
 */

#include <ergodic_exploration/dynamic_window.hpp>
#include <ergodic_exploration/numerics.hpp>

namespace ergodic_exploration
{
using arma::span;

DynamicWindow::DynamicWindow(const Collision& collision, double dt, double horizon,
                             double acc_dt, double acc_lim_x, double acc_lim_y,
                             double acc_lim_th, double max_vel_x, double min_vel_x,
                             double max_vel_y, double min_vel_y, double max_rot_vel,
                             double min_rot_vel, unsigned int vx_samples,
                             unsigned int vy_samples, unsigned int vth_samples)
  : collision_(collision)
  , dt_(dt)
  , horizon_(horizon)
  , acc_dt_(acc_dt)
  , acc_lim_x_(acc_lim_x)
  , acc_lim_y_(acc_lim_y)
  , acc_lim_th_(acc_lim_th)
  , max_vel_x_(max_vel_x)
  , min_vel_x_(min_vel_x)
  , max_vel_y_(max_vel_y)
  , min_vel_y_(min_vel_y)
  , max_rot_vel_(max_rot_vel)
  , min_rot_vel_(min_rot_vel)
  , vx_samples_(vx_samples)
  , vy_samples_(vy_samples)
  , vth_samples_(vth_samples)
  , steps_(static_cast<unsigned int>(std::abs(horizon / dt)))
{
  if (vx_samples_ == 0)
  {
    std::cout << "vx samples set to 0 but need at least 1... setting this to 1"
              << std::endl;
    vx_samples_ = 1;
  }

  if (vy_samples_ == 0)
  {
    std::cout << "vy samples set to 0 but need at least 1... setting this to 1"
              << std::endl;
    vy_samples_ = 1;
  }

  if (vth_samples_ == 0)
  {
    std::cout << "vth samples set to 0 but need at least 1... setting this to 1"
              << std::endl;
    vth_samples_ = 1;
  }
}

tuple<bool, vec> DynamicWindow::control(const GridMap& grid, const vec& x0, const vec& vb,
                                        const vec& vref) const
{
  // window limits and discretization
  const std::tuple<vec, vec> window_config = window(vb);

  // Search velocity space
  auto min_cost = std::numeric_limits<double>::max();

  vec u_opt(3, arma::fill::zeros);  // optimal twist
  vec u(3, arma::fill::zeros);      // sample twist

  auto vx = std::get<0>(window_config)(0);
  for (unsigned int i = 0; i < vx_samples_; i++)
  {
    auto vy = std::get<0>(window_config)(1);
    for (unsigned int j = 0; j < vy_samples_; j++)
    {
      auto w = std::get<0>(window_config)(2);
      for (unsigned int k = 0; k < vth_samples_; k++)
      {
        u(0) = vx;
        u(1) = vy;
        u(2) = w;

        const auto cost = objective(grid, x0, vref, u);
        // minimize cost
        if (cost < min_cost)
        {
          min_cost = cost;
          u_opt = u;
        }

        w += std::get<1>(window_config)(2);
      }  // end w loop
      vy += std::get<1>(window_config)(1);
    }  // end vy loop
    vx += std::get<1>(window_config)(0);
  }  // end vx loop

  if (almost_equal(min_cost, std::numeric_limits<double>::max()))
  {
    std::cout << "DWA Failed! Not even 1 solution found" << std::endl;
    return std::make_tuple(false, u_opt);
  }

  return std::make_tuple(true, u_opt);
}

tuple<bool, vec> DynamicWindow::control(const GridMap& grid, const vec& x0, const vec& vb,
                                        const mat& xt_ref, double dt_ref) const
{
  // window limits and discretization
  const tuple<vec, vec> window_config = window(vb);

  // Time parameterization of the reference trajectory
  const double tf = static_cast<double>(xt_ref.n_cols) * dt_ref;

  // Search velocity space
  auto min_cost = std::numeric_limits<double>::max();

  vec u_opt(3, arma::fill::zeros);  // optimal twist
  vec u(3, arma::fill::zeros);      // sample twist

  auto vx = std::get<0>(window_config)(0);
  for (unsigned int i = 0; i < vx_samples_; i++)
  {
    auto vy = std::get<0>(window_config)(1);
    for (unsigned int j = 0; j < vy_samples_; j++)
    {
      auto w = std::get<0>(window_config)(2);
      for (unsigned int k = 0; k < vth_samples_; k++)
      {
        u(0) = vx;
        u(1) = vy;
        u(2) = w;

        const auto cost = objective(grid, x0, u, xt_ref, tf);
        // minimize cost
        if (cost < min_cost)
        {
          min_cost = cost;
          u_opt = u;
        }

        w += std::get<1>(window_config)(2);
      }  // end w loop
      vy += std::get<1>(window_config)(1);
    }  // end vy loop
    vx += std::get<1>(window_config)(0);
  }  // end vx loop

  if (almost_equal(min_cost, std::numeric_limits<double>::max()))
  {
    std::cout << "DWA Failed! Not even 1 solution found" << std::endl;
    return std::make_tuple(false, u_opt);
  }

  return std::make_tuple(true, u_opt);
}

tuple<vec, vec> DynamicWindow::window(const vec& vb) const
{
  // Window bounds and discretization
  vec vel_upper(3, arma::fill::zeros);  // twist upper limits
  vec vel_lower(3, arma::fill::zeros);  // twist lower limits
  vec delta_vb(3, arma::fill::zeros);   // twist discretization

  vel_lower(0) = std::max(vb(0) - acc_lim_x_ * acc_dt_, min_vel_x_);
  vel_upper(0) = std::min(vb(0) + acc_lim_x_ * acc_dt_, max_vel_x_);

  vel_lower(1) = std::max(vb(1) - acc_lim_y_ * acc_dt_, min_vel_y_);
  vel_upper(1) = std::min(vb(1) + acc_lim_y_ * acc_dt_, max_vel_y_);

  vel_lower(2) = std::max(vb(2) - acc_lim_th_ * acc_dt_, min_rot_vel_);
  vel_upper(2) = std::min(vb(2) + acc_lim_th_ * acc_dt_, max_rot_vel_);

  // std::cout << "Window " << std::endl;
  // std::cout << "vx: [ " << vel_lower(0) << ", " << vel_upper(0) << " ] \n" ;
  // std::cout << "vy: [ " << vel_lower(1) << ", " << vel_upper(1) << " ] \n";
  // std::cout << "w: [ " << vel_lower(2) << ", " << vel_upper(2) << " ]" << std::endl;

  if (vx_samples_ > 1)
  {
    delta_vb(0) = (vel_upper(0) - vel_lower(0)) / static_cast<double>(vx_samples_ - 1);
  }

  if (vy_samples_ > 1)
  {
    delta_vb(1) = (vel_upper(1) - vel_lower(1)) / static_cast<double>(vy_samples_ - 1);
  }

  if (vth_samples_ > 1)
  {
    delta_vb(2) = (vel_upper(2) - vel_lower(2)) / static_cast<double>(vth_samples_ - 1);
  }

  // std::cout << "Window discretization " << std::endl;
  // std::cout << "dvx: " << delta_vb(0) << std::endl;
  // std::cout << "dvy: " << delta_vb(1) << std::endl;
  // std::cout << "dw: " << delta_vb(2) << std::endl;

  return std::make_tuple(vel_lower, delta_vb);
}

double DynamicWindow::objective(const GridMap& grid, const vec& x0, const vec& vref,
                                const vec& u) const
{
  // follow constant twist
  vec pose = x0;
  for (unsigned int i = 0; i < steps_; i++)
  {
    pose = integrate_twist(pose, u, dt_);
    pose(2) = normalize_angle_PI(pose(2));

    if (collision_.collisionCheck(grid, pose))
    {
      return std::numeric_limits<double>::max();
    }
  }

  // control error;
  const vec cntrl_error = vref - u;
  return dot(cntrl_error, cntrl_error);
}

double DynamicWindow::objective(const GridMap& grid, const vec& x0, const vec& u,
                                const mat& xt_ref, double tf) const
{
  vec pose = x0;
  auto t = 0.0;
  auto cost = 0.0;

  // follow constant twist
  for (unsigned int i = 0; i < steps_; i++)
  {
    pose = integrate_twist(pose, u, dt_);
    pose(2) = normalize_angle_PI(pose(2));

    if (collision_.collisionCheck(grid, pose))
    {
      return std::numeric_limits<double>::max();
    }

    // index into reference trajectory
    const auto j = static_cast<unsigned int>(std::round((xt_ref.n_cols - 1) * t / tf));

    cost += arma::norm(xt_ref(span(0, 1), span(j, j)) - pose.rows(0, 1));
    cost += std::abs(normalize_angle_PI(normalize_angle_PI(xt_ref(2, j)) - pose(2)));

    t += dt_;
  }

  return cost;
}
}  // namespace ergodic_exploration
