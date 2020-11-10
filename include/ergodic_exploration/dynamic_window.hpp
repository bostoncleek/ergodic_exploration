/**
 * @file dynamic_window.hpp
 * @author Boston Cleek
 * @date 5 Nov 2020
 * @brief Dynamic window control
 */

#pragma once

#include <cmath>
#include <limits>

// #include <geometry_msgs/PoseStamped.h>

#include <ergodic_exploration/grid.hpp>
#include <ergodic_exploration/collision.hpp>

namespace ergodic_exploration
{
template <class ModelT>
class DynamicWindow
{
public:
  DynamicWindow(const ModelT& model, double dt, double horizon, double acc_lim_x,
                double acc_lim_y, double acc_lim_th, double max_vel_trans,
                double min_vel_trans, double max_vel_x, double min_vel_x,
                double max_vel_y, double min_vel_y, double max_rot_vel,
                double min_rot_vel, unsigned int vx_samples, unsigned int vy_samples,
                unsigned int vth_samples);

  vec control(Collision collision, const GridMap& grid, const vec& xs, const vec& xg,
              const vec& vb);

  double trajectory(Collision collision, const GridMap& grid, const vec& xs,
                    const vec& xg, const vec u);

  double loss(Collision collision, const GridMap& grid, const vec& x, const vec& xg,
              const vec u);

private:
  ModelT model_;                               // robot's dynamic model
  double dt_;                                  // time step in integration
  double horizon_;                             // control horizon
  double acc_lim_x_, acc_lim_y_, acc_lim_th_;  // acceleration limits
  double max_vel_trans_, min_vel_trans_;       // translational velocity limits
  double max_vel_x_, min_vel_x_;               // velocity limits in x-direction
  double max_vel_y_, min_vel_y_;               // velocity limits in y-direction
  double max_rot_vel_, min_rot_vel_;           // rotational velocity limits
  unsigned int vx_samples_;   // number of velocity samples in x-direction
  unsigned int vy_samples_;   // number of velocity samples in y-direction
  unsigned int vth_samples_;  // number of angular velocity samples
  RungeKutta rk4_;            // integrator
};

template <class ModelT>
DynamicWindow<ModelT>::DynamicWindow(const ModelT& model, double dt, double horizon,
                                     double acc_lim_x, double acc_lim_y,
                                     double acc_lim_th, double max_vel_trans,
                                     double min_vel_trans, double max_vel_x,
                                     double min_vel_x, double max_vel_y, double min_vel_y,
                                     double max_rot_vel, double min_rot_vel,
                                     unsigned int vx_samples, unsigned int vy_samples,
                                     unsigned int vth_samples)
  : model_(model)
  , dt_(dt)
  , horizon_(horizon)
  , acc_lim_x_(acc_lim_x)
  , acc_lim_y_(acc_lim_y)
  , acc_lim_th_(acc_lim_th)
  , max_vel_trans_(max_vel_trans)
  , min_vel_trans_(min_vel_trans)
  , max_vel_x_(max_vel_x)
  , min_vel_x_(min_vel_x)
  , max_vel_y_(max_vel_y)
  , min_vel_y_(min_vel_y)
  , max_rot_vel_(max_rot_vel)
  , min_rot_vel_(min_rot_vel)
  , vx_samples_(vx_samples)
  , vy_samples_(vy_samples)
  , vth_samples_(vth_samples)
  , rk4_(dt)
{
}

template <class ModelT>
vec DynamicWindow<ModelT>::control(Collision collision, const GridMap& grid,
                                   const vec& xs, const vec& xg, const vec& vb)
{
  // TODO: Add frequency parameter to determine the delta t the max accel is applied
  // Compose dynamic window
  const auto vdx_low = std::max(vb(0) - acc_lim_x_ * dt_, min_vel_x_);
  const auto vdx_high = std::min(vb(0) + acc_lim_x_ * dt_, max_vel_x_);

  const auto vdy_low = std::max(vb(1) - acc_lim_y_ * dt_, min_vel_y_);
  const auto vdy_high = std::min(vb(1) + acc_lim_y_ * dt_, max_vel_y_);

  const auto wd_low = std::max(vb(2) - acc_lim_th_ * dt_, min_rot_vel_);
  const auto wd_high = std::min(vb(2) + acc_lim_th_ * dt_, max_rot_vel_);

  // TODO: make sure samples are at least 1
  // Search space
  const auto dx = (vdx_high - vdx_low) / static_cast<double>(vx_samples_);
  const auto dy = (vdy_high - vdy_low) / static_cast<double>(vy_samples_);
  const auto dw = (wd_high - wd_low) / static_cast<double>(vth_samples_);

  // TODO: randomly sample and/or add thread pool
  // Search velocity space
  vec u_opt(model_.action_space);
  auto max_cost = std::numeric_limits<double>::min();
  auto vx = vdx_low;
  auto vy = vdy_low;
  auto w = wd_low;
  for (unsigned int i = 0; i < vx_samples_; i++)
  {
    for (unsigned int j = 0; j < vy_samples_; j++)
    {
      for (unsigned int k = 0; k < vy_samples_; k++)
      {
        // Froward simulate trajectory and compose cost
        const vec u = { vx, vy, w };
        const auto cost = trajectory(collision, grid, xs, xg, u);

        // trajectory with highst cost
        if (cost > max_cost)
        {
          max_cost = cost;
          u_opt = u;
        }

        w += dw;
      }  // end w loop
      vy += dy;
    }  // end vy loop
    vx += dx;
  }  // end vx loop

  return u_opt;
}

template <class ModelT>
double DynamicWindow<ModelT>::trajectory(Collision collision, const GridMap& grid,
                                         const vec& xs, const vec& xg, const vec u)
{
  vec x = xs;
  const auto steps = static_cast<unsigned int>(horizon_ / std::abs(dt_));

  auto cost = 0.0;
  for (unsigned int i = 0; i < steps; i++)
  {
    x = rk4_.step(model_, x, u);
    cost += loss(collision, grid, x, xg, u);
  }

  return cost;
}

template <class ModelT>
double DynamicWindow<ModelT>::loss(Collision collision, const GridMap& grid, const vec& x,
                                   const vec& xg, const vec u)
{
}

}  // namespace ergodic_exploration
