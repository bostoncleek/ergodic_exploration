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
DynamicWindow::DynamicWindow(double dt, double horizon, double frequency,
                             double acc_lim_x, double acc_lim_y, double acc_lim_th,
                             double max_vel_x, double min_vel_x, double max_vel_y,
                             double min_vel_y, double max_rot_vel, double min_rot_vel,
                             unsigned int vx_samples, unsigned int vy_samples,
                             unsigned int vth_samples)
  : dt_(dt)
  , horizon_(horizon)
  , frequency_(frequency)
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
{
}

vec DynamicWindow::control(const Collision& collision, const GridMap& grid, const vec& x,
                           const vec& vb, const vec& vref)
{
  const auto cntrl_dt = 1.0 / frequency_;
  // const auto cntrl_dt = 0.5;
  // const auto cntrl_dt = 1.0/ 10.0;

  const auto vdx_low = std::max(vb(0) - acc_lim_x_ * cntrl_dt, min_vel_x_);
  const auto vdx_high = std::min(vb(0) + acc_lim_x_ * cntrl_dt, max_vel_x_);

  const auto vdy_low = std::max(vb(1) - acc_lim_y_ * cntrl_dt, min_vel_y_);
  const auto vdy_high = std::min(vb(1) + acc_lim_y_ * cntrl_dt, max_vel_y_);

  const auto wd_low = std::max(vb(2) - acc_lim_th_ * cntrl_dt, min_rot_vel_);
  const auto wd_high = std::min(vb(2) + acc_lim_th_ * cntrl_dt, max_rot_vel_);

  // std::cout << "Window " << std::endl;
  // std::cout << "vx: [ " << vdx_low << ", " << vdx_high << " ]" << std::endl;
  // std::cout << "vy: [ " << vdy_low << ", " << vdy_high << " ]" << std::endl;
  // std::cout << "w: [ " << wd_low << ", " << wd_high << " ]" << std::endl;

  // TODO: make sure no division by 0
  const auto dvx = (vdx_high - vdx_low) / static_cast<double>(vx_samples_ - 1);
  const auto dvy = (vdy_high - vdy_low) / static_cast<double>(vy_samples_ - 1);
  const auto dw = (wd_high - wd_low) / static_cast<double>(vth_samples_ - 1);

  // // number of steps in each rollout
  // const auto steps = static_cast<unsigned int>(horizon_ / dt_);

  // std::cout << "Window discretization " << std::endl;
  // std::cout << "dvx: " << dvx << std::endl;
  // std::cout << "dvy: " << dvy << std::endl;
  // std::cout << "dw: " << dw << std::endl;

  // Search velocity space
  vec u_opt(3, arma::fill::zeros);
  auto min_cost = std::numeric_limits<double>::max();

  bool soln_found = false;

  vec u(3);
  auto vx = vdx_low;
  for (unsigned int i = 0; i < vx_samples_; i++)
  {
    auto vy = vdy_low;
    for (unsigned int j = 0; j < vy_samples_; j++)
    {
      auto w = wd_low;
      for (unsigned int k = 0; k < vth_samples_; k++)
      {
        u(0) = vx;
        u(1) = vy;
        u(2) = w;

        auto cost = 0.0;
        if (!objective(cost, collision, grid, x, vref, u))
        {
          // control error;
          const vec error = vref - u;
          cost += dot(error, error);

          const auto v_trans =
              std::sqrt(max_vel_x_ * max_vel_x_ + max_vel_y_ * max_vel_y_);
          const auto v_sample = std::sqrt(u(0) * u(0) + u(1) * u(1));

          cost += v_trans - v_sample;

          // cost += std::abs(u(2) - max_rot_vel_);
          // cost += std::abs(u(2) - max_rot_vel);

          // minimize cost
          if (cost < min_cost)
          {
            min_cost = cost;
            u_opt = u;
            soln_found = true;
          }
        }

        w += dw;
      }  // end w loop
      vy += dvy;
    }  // end vy loop
    vx += dvx;
  }  // end vx loop

  // simple recovery behavior rotate in place
  if (!soln_found)
  {
    std::cout << "DWA Failed! Not even 1 solution found" << std::endl;
  }

  // std::cout << "min cost " << min_cost << std::endl;

  return u_opt;
}

bool DynamicWindow::objective(double& loss, const Collision& collision,
                              const GridMap& grid, const vec& x, const vec& vref,
                              const vec& u)
{
  // number of steps in each rollout
  const auto steps = static_cast<unsigned int>(horizon_ / dt_);

  vec pose = x;
  for (unsigned int t = 0; t < steps; t++)
  {
    pose = integrate_twist(pose, u, dt_);

    auto dmin = std::numeric_limits<double>::max();
    if (collision.minDistance(dmin, grid, pose) == CollisionMsg::crash)
    {
      return true;
    }

    loss += (1.0 / dmin);
  }

  // // Clearance
  // // assume obstacles are very far away
  // auto dmin = std::numeric_limits<double>::max();
  // if (collision.minDistance(dmin, grid, x) == CollisionMsg::crash)
  // {
  //   // collison receives maximum penalty
  //   // return std::numeric_limits<double>::max();
  //   return true;
  // }

  // // control error;
  // const vec error = vref - u;
  //
  // // TODO: add weight parameters
  // // return 100.0 * (1.0 / dmin) + 1.0 * dot(error, error);
  // loss = (1.0 / dmin);
  return false;
}
}  // namespace ergodic_exploration
