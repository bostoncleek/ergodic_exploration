/**
 * @file dynamic_window.hpp
 * @author Boston Cleek
 * @date 5 Nov 2020
 * @brief Dynamic window control
 */

#pragma once

#include <cmath>
#include <limits>

#include <ergodic_exploration/grid.hpp>
#include <ergodic_exploration/collision.hpp>

namespace ergodic_exploration
{
template <class ModelT>
class DynamicWindow
{
public:
  DynamicWindow(const ModelT& model, double dt, double horizon, double frequency,
                double acc_lim_x, double acc_lim_y, double acc_lim_th,
                double max_vel_trans, double min_vel_trans, double max_vel_x,
                double min_vel_x, double max_vel_y, double min_vel_y, double max_rot_vel,
                double min_rot_vel, unsigned int vx_samples, unsigned int vy_samples,
                unsigned int vth_samples);

  /**
   * @brief Update control
   * @param collision - collision detector
   * @param grid - grid map
   * @param x0 - current state [x, y, theta]
   * @param vb - current twist [vx, vy, w]
   * @param vref - desired twist to follow [vx, vy, w]
   * @return twist
   */
  vec control(const Collision& collision, const GridMap& grid, const vec& x0,
              const vec& vb, const vec& vref);

private:
  bool trajectory(double& cost, const Collision& collision, const GridMap& grid,
                  const vec& x0, const vec& vref, const vec& u);

  bool collisionLoss(double& loss, const Collision& collision, const GridMap& grid,
                     const vec& x);

  double controlLoss(const vec& vref, const vec& u);

private:
  ModelT model_;                               // robot's dynamic model
  double dt_;                                  // time step in integration
  double horizon_;                             // control horizon
  double frequency_;                           // control loop frequency
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
                                     double frequency, double acc_lim_x, double acc_lim_y,
                                     double acc_lim_th, double max_vel_trans,
                                     double min_vel_trans, double max_vel_x,
                                     double min_vel_x, double max_vel_y, double min_vel_y,
                                     double max_rot_vel, double min_rot_vel,
                                     unsigned int vx_samples, unsigned int vy_samples,
                                     unsigned int vth_samples)
  : model_(model)
  , dt_(dt)
  , horizon_(horizon)
  , frequency_(frequency)
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
vec DynamicWindow<ModelT>::control(const Collision& collision, const GridMap& grid,
                                   const vec& x0, const vec& vb, const vec& vref)
{
  // const auto cntrl_dt = 1.0 / frequency_;
  const auto cntrl_dt = 0.2;

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

  const vec vx_vec = arma::linspace(vdx_low, vdx_high, vx_samples_ + 1);
  const vec vy_vec = arma::linspace(vdy_low, vdy_high, vy_samples_ + 1);
  const vec w_vec = arma::linspace(wd_low, wd_high, vth_samples_ + 1);

  // vx_vec.print("vxs:");
  // vy_vec.print("vys:");
  // w_vec.print("ws:");

  // TODO: randomly sample and/or add thread pool
  // Search velocity space
  vec u_opt(model_.action_space, arma::fill::zeros);
  auto min_cost = std::numeric_limits<double>::max();

  bool soln_found = false;

  for (const auto vx : vx_vec)
  {
    for (const auto vy : vy_vec)
    {
      for (const auto w : w_vec)
      {
        // Froward simulate trajectory and compose cost
        const vec u = { vx, vy, w };

        auto cost = 0.0;
        if (!trajectory(cost, collision, grid, x0, vref, u))
        {
          // minimize cost
          if (cost < min_cost)
          {
            min_cost = cost;
            u_opt = u;
            soln_found = true;
          }
        }

      }  // end w loop
    }    // end vy loop
  }      // end vx loop

  if (!soln_found)
  {
    std::cout << "DWA Failed! Not even 1 solution found" << std::endl;
  }

  // std::cout << "min cost " << min_cost << std::endl;

  return u_opt;
}

template <class ModelT>
bool DynamicWindow<ModelT>::trajectory(double& cost, const Collision& collision,
                                       const GridMap& grid, const vec& x0,
                                       const vec& vref, const vec& u)
{
  vec x = x0;
  const auto steps = static_cast<unsigned int>(horizon_ / std::abs(dt_));

  // Check for collisions and distance to obstacles
  auto obs_loss = 0.0;
  for (unsigned int i = 0; i < steps; i++)
  {
    x = rk4_.step(model_, x, u);

    auto loss = 0.0;
    if (collisionLoss(loss, collision, grid, x))
    {
      return true;
    }
    obs_loss += loss;
  }

  // How closely follows reference control
  auto cntrl_loss = controlLoss(vref, u);

  cost = obs_loss + cntrl_loss;

  std::cout << "obstacle loss: " << obs_loss << std::endl;
  std::cout << "control loss: " << cntrl_loss << std::endl;
  // std::cout << "cost: " << cost << std::endl;

  return false;
}

template <class ModelT>
bool DynamicWindow<ModelT>::collisionLoss(double& loss, const Collision& collision,
                                          const GridMap& grid, const vec& x)
{
  // Clearance
  auto dmin = std::numeric_limits<double>::max();  // assume obstacles are very far away
  if (collision.minDistance(dmin, grid, x) == CollisionMsg::crash)
  {
    // std::cout << "collision" << std::endl;
    return true;
  }

  loss = 1.0 / dmin;

  return false;
}

template <class ModelT>
double DynamicWindow<ModelT>::controlLoss(const vec& vref, const vec& u)
{
  // TODO: add weights as param
  mat W(3, 3, arma::fill::zeros);
  W(0, 0) = 1e2;
  W(1, 1) = 1e2;
  W(2, 2) = 1e2;

  const vec error = vref - u;
  return dot(error.t() * W, error);
}

}  // namespace ergodic_exploration
