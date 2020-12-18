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

vec DynamicWindow::control(const GridMap& grid, const vec& x0, const vec& vb,
                           const vec& vref) const
{
  // Window bounds and discretization
  vec vel_lower(3, arma::fill::zeros);  // twist lower limits
  vec delta_vb(3, arma::fill::zeros);   // twist discretization
  vec u_opt(3, arma::fill::zeros);      // optimal twist
  vec u(3);                             // sample twist

  window(vel_lower, delta_vb, vb);

  // Search velocity space
  auto min_cost = std::numeric_limits<double>::max();

  bool soln_found = false;
  auto vx = vel_lower(0);
  for (unsigned int i = 0; i < vx_samples_; i++)
  {
    auto vy = vel_lower(1);
    for (unsigned int j = 0; j < vy_samples_; j++)
    {
      auto w = vel_lower(2);
      for (unsigned int k = 0; k < vth_samples_; k++)
      {
        u(0) = vx;
        u(1) = vy;
        u(2) = w;

        auto cost = 0.0;
        if (!objective(cost, grid, x0, vref, u))
        {
          // minimize cost
          if (cost < min_cost)
          {
            min_cost = cost;
            u_opt = u;
            soln_found = true;
          }
        }

        w += delta_vb(2);
      }  // end w loop
      vy += delta_vb(1);
    }  // end vy loop
    vx += delta_vb(0);
  }  // end vx loop

  // simple recovery behavior rotate in place
  if (!soln_found)
  {
    std::cout << "DWA Failed! Not even 1 solution found" << std::endl;
  }

  return u_opt;
}

bool DynamicWindow::control(vec& u_opt, const GridMap& grid, const vec& x0, const vec& vb,
                            const mat& xt_ref, double dt_ref) const
{
  // Window bounds and discretization
  vec vel_lower(3, arma::fill::zeros);  // twist lower limits
  vec delta_vb(3, arma::fill::zeros);   // twist discretization

  window(vel_lower, delta_vb, vb);

  // Time parameterization of the reference trajectory
  const double tf = static_cast<double>(xt_ref.n_cols) * dt_ref;

  // Search velocity space
  auto min_cost = std::numeric_limits<double>::max();

  bool soln_found = false;

  vec u(3);
  auto vx = vel_lower(0);
  for (unsigned int i = 0; i < vx_samples_; i++)
  {
    auto vy = vel_lower(1);
    for (unsigned int j = 0; j < vy_samples_; j++)
    {
      auto w = vel_lower(2);
      for (unsigned int k = 0; k < vth_samples_; k++)
      {
        u(0) = vx;
        u(1) = vy;
        u(2) = w;

        auto cost = 0.0;
        if (!objective(cost, grid, x0, u, xt_ref, tf))
        {
          // minimize cost
          if (cost < min_cost)
          {
            min_cost = cost;
            u_opt = u;
            soln_found = true;
          }
        }

        w += delta_vb(2);
      }  // end w loop
      vy += delta_vb(1);
    }  // end vy loop
    vx += delta_vb(0);
  }  // end vx loop

  // simple recovery behavior rotate in place
  if (!soln_found)
  {
    std::cout << "DWA Failed! Not even 1 solution found" << std::endl;
    u_opt.zeros();
  }

  return soln_found;
}

void DynamicWindow::window(vec& vel_lower, vec& delta_vb, const vec& vb) const
{
  vec vel_upper(3, arma::fill::zeros);  // twist upper limits

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
}

bool DynamicWindow::objective(double& cost, const GridMap& grid, const vec& x0,
                              const vec& vref, const vec& u) const
{
  vec pose = x0;
  // follow constant twist
  const vec delta = integrate_twist(pose, u, dt_);
  for (unsigned int i = 0; i < steps_; i++)
  {
    pose += delta;
    pose(2) = normalize_angle_PI(pose(2));

    auto dmin = std::numeric_limits<double>::max();
    if (collision_.minDistance(dmin, grid, pose) == CollisionMsg::crash)
    {
      return true;
    }
  }

  // control error;
  const vec cntrl_error = vref - u;
  cost += dot(cntrl_error, cntrl_error);

  return false;
}

bool DynamicWindow::objective(double& cost, const GridMap& grid, const vec& x0,
                              const vec& u, const mat& xt_ref, double tf) const
{
  vec pose = x0;
  double t = 0.0;

  // follow constant twist
  const vec delta = integrate_twist(pose, u, dt_);
  for (unsigned int i = 0; i < steps_; i++)
  {
    pose += delta;
    pose(2) = normalize_angle_PI(pose(2));

    auto dmin = std::numeric_limits<double>::max();
    if (collision_.minDistance(dmin, grid, pose) == CollisionMsg::crash)
    {
      return true;
    }

    // index into reference trajectory
    const auto j = static_cast<unsigned int>(std::round((xt_ref.n_cols - 1) * t / tf));

    cost += arma::norm(xt_ref(span(0, 1), span(j, j)) - pose.rows(0, 1));
    cost += std::abs(normalize_angle_PI(normalize_angle_PI(xt_ref(2, j)) - pose(2)));

    t += dt_;
  }

  return false;
}

}  // namespace ergodic_exploration
