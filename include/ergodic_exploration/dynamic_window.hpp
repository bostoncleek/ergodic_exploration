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

  vec control(const Collision& collision, const GridMap& grid, const vec& xs,
              const vec& xg, const vec& vb);

  bool trajectory(double& cost, const Collision& collision, const GridMap& grid,
                  const vec& xs, const vec& xg, const vec& u);

  bool objective(double& loss, const Collision& collision, const GridMap& grid,
                 const vec& x, const vec& xg, const vec& u);

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
vec DynamicWindow<ModelT>::control(const Collision& collision, const GridMap& grid,
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
  // const auto dx = (vdx_high - vdx_low) / static_cast<double>(vx_samples_);
  // const auto dy = (vdy_high - vdy_low) / static_cast<double>(vy_samples_);
  // const auto dw = (wd_high - wd_low) / static_cast<double>(vth_samples_);

  // std::cout << "Window " << std::endl;
  // std::cout << "vx: [ " << vdx_low << ", " << vdx_high << " ]" << std::endl;
  // std::cout << "vy: [ " << vdy_low << ", " << vdy_high << " ]" << std::endl;
  // std::cout << "w: [ " << wd_low << ", " << wd_high << " ]" << std::endl;
  // std::cout << "Window Discretization " << std::endl;
  // std::cout << "dx: " << dx  << std::endl;
  // std::cout << "dy: " << dy  << std::endl;
  // std::cout << "dw: " << dw  << std::endl;

  vec vx_vec = arma::linspace(vdx_low, vdx_high, vx_samples_ + 1);
  vec vy_vec = arma::linspace(vdy_low, vdy_high, vy_samples_ + 1);
  vec w_vec = arma::linspace(wd_low, wd_high, vth_samples_ + 1);

  // vx_vec.print("vxs:");
  // vy_vec.print("vys:");
  // w_vec.print("ws:");

  // TODO: randomly sample and/or add thread pool
  // Search velocity space
  vec u_opt(model_.action_space, arma::fill::zeros);
  auto min_cost = std::numeric_limits<double>::max();

  for (const auto vx : vx_vec)
  {
    for (const auto vy : vy_vec)
    {
      for (const auto w : w_vec)
      {
        // Froward simulate trajectory and compose cost
        const vec u = { vx, vy, w };
        // u.print("Control sample:");

        auto cost = 0.0;
        if(!trajectory(cost, collision, grid, xs, xg, u))
        {
          // u.print("u:");
          // std::cout << "cost " << cost << std::endl;

          // minimize cost
          if (cost < min_cost)
          {
            min_cost = cost;
            u_opt = u;
          }
        }

      }  // end w loop
    }  // end vy loop
  }  // end vx loop

  // std::cout << "min cost " << min_cost << std::endl;

  return u_opt;
}

template <class ModelT>
bool DynamicWindow<ModelT>::trajectory(double& cost, const Collision& collision,
                                       const GridMap& grid, const vec& xs, const vec& xg,
                                       const vec& u)
{
  vec x = xs;
  const auto steps = static_cast<unsigned int>(horizon_ / std::abs(dt_));

  for (unsigned int i = 0; i < steps; i++)
  {
    x = rk4_.step(model_, x, u);

    auto loss = 0.0;
    if (objective(loss, collision, grid, x, xg, u))
    {
      return true;
    }

    cost += loss;
  }

  return false;
}

template <class ModelT>
bool DynamicWindow<ModelT>::objective(double& loss, const Collision& collision,
                                      const GridMap& grid, const vec& x, const vec& xg,
                                      const vec& u)
{
  // TODO: Add weight parameters
  // Heading
  const auto h =
      std::abs(normalize_angle_PI(normalize_angle_PI(xg(2)) - normalize_angle_PI(x(2))));

  // std::cout << "heading term " << h << std::endl;

  // Distance to goal
  const auto dgoal = distance(x(0), x(1), xg(0), xg(1));

  // std::cout << "distance to goal: " << dgoal << std::endl;


  // Clearance
  auto dmin = std::numeric_limits<double>::max(); // assume obstacles are very far away
  if (collision.minDistance(dmin, grid, x) == CollisionMsg::crash)
  {
    // std::cout << "collision" << std::endl;
    return true;
  }

  // // TODO: add threshold parameter for distance to goal
  // // TODO: make sure u is clamped
  // // TODO: add rotational velocity to cost
  // Velocity loss
  // auto vel_loss = 0.0;
  // if (dgoal < 0.2)
  // {
  //   vel_loss = norm(u.rows(0,1), 2) / max_vel_trans_;
  // }
  // else
  // {
  //   vel_loss = 1.0 - norm(u.rows(0,1), 2) / max_vel_trans_;
  // }

  // std::cout << "obstacle term " << (1.0 / dmin) << std::endl;

  // normalize each component
  vec loss_vec = {h, dgoal, (1.0 / dmin)};
  loss_vec /= sum(loss_vec);

  // weigh each component
  vec w = {1.5, 2.0, 10.0};

  loss = dot(loss_vec, w);


  // heading_loss + obstacle_loss + velocity_loss
  // loss = h + dgoal + (1.0 / dmin) + vel_loss;
  // loss = h + dgoal + vel_loss;
  // loss = h + dgoal + (1.0 / dmin);
  //loss = h;
  // loss = vel_loss;

  // std::cout << "updated loss: " << loss << std::endl;

  return false;
}

}  // namespace ergodic_exploration
