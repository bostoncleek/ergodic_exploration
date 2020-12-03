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
  const auto vdx_low = std::max(vb(0) - acc_lim_x_ * acc_dt_, min_vel_x_);
  const auto vdx_high = std::min(vb(0) + acc_lim_x_ * acc_dt_, max_vel_x_);

  const auto vdy_low = std::max(vb(1) - acc_lim_y_ * acc_dt_, min_vel_y_);
  const auto vdy_high = std::min(vb(1) + acc_lim_y_ * acc_dt_, max_vel_y_);

  const auto wd_low = std::max(vb(2) - acc_lim_th_ * acc_dt_, min_rot_vel_);
  const auto wd_high = std::min(vb(2) + acc_lim_th_ * acc_dt_, max_rot_vel_);

  // std::cout << "Window " << std::endl;
  // std::cout << "vx: [ " << vdx_low << ", " << vdx_high << " ]" << std::endl;
  // std::cout << "vy: [ " << vdy_low << ", " << vdy_high << " ]" << std::endl;
  // std::cout << "w: [ " << wd_low << ", " << wd_high << " ]" << std::endl;

  auto dvx = 0.0;
  auto dvy = 0.0;
  auto dw = 0.0;

  if (vx_samples_ > 1)
  {
    dvx = (vdx_high - vdx_low) / static_cast<double>(vx_samples_ - 1);
  }

  if (vy_samples_ > 1)
  {
    dvy = (vdy_high - vdy_low) / static_cast<double>(vy_samples_ - 1);
  }

  if (vth_samples_ > 1)
  {
    dw = (wd_high - wd_low) / static_cast<double>(vth_samples_ - 1);
  }

  // TODO: make sure no division by 0
  // const auto dvx = (vdx_high - vdx_low) / static_cast<double>(vx_samples_ - 1);
  // const auto dvy = (vdy_high - vdy_low) / static_cast<double>(vy_samples_ - 1);
  // const auto dw = (wd_high - wd_low) / static_cast<double>(vth_samples_ - 1);

  // std::cout << "Window discretization " << std::endl;
  // std::cout << "dvx: " << dvx << std::endl;
  // std::cout << "dvy: " << dvy << std::endl;
  // std::cout << "dw: " << dw << std::endl;

  // Search velocity space
  auto min_cost = std::numeric_limits<double>::max();

  bool soln_found = false;

  vec u_opt(3, arma::fill::zeros);
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

        // u.print("u tests");

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

vec DynamicWindow::control(const GridMap& grid, const vec& x0, const vec& vb,
                           const mat& xt_ref, double delta) const
{
  const auto vdx_low = std::max(vb(0) - acc_lim_x_ * acc_dt_, min_vel_x_);
  const auto vdx_high = std::min(vb(0) + acc_lim_x_ * acc_dt_, max_vel_x_);

  const auto vdy_low = std::max(vb(1) - acc_lim_y_ * acc_dt_, min_vel_y_);
  const auto vdy_high = std::min(vb(1) + acc_lim_y_ * acc_dt_, max_vel_y_);

  const auto wd_low = std::max(vb(2) - acc_lim_th_ * acc_dt_, min_rot_vel_);
  const auto wd_high = std::min(vb(2) + acc_lim_th_ * acc_dt_, max_rot_vel_);

  // std::cout << "Window " << std::endl;
  // std::cout << "vx: [ " << vdx_low << ", " << vdx_high << " ]" << std::endl;
  // std::cout << "vy: [ " << vdy_low << ", " << vdy_high << " ]" << std::endl;
  // std::cout << "w: [ " << wd_low << ", " << wd_high << " ]" << std::endl;

  auto dvx = 0.0;
  auto dvy = 0.0;
  auto dw = 0.0;

  if (vx_samples_ > 1)
  {
    dvx = (vdx_high - vdx_low) / static_cast<double>(vx_samples_ - 1);
  }

  if (vy_samples_ > 1)
  {
    dvy = (vdy_high - vdy_low) / static_cast<double>(vy_samples_ - 1);
  }

  if (vth_samples_ > 1)
  {
    dw = (wd_high - wd_low) / static_cast<double>(vth_samples_ - 1);
  }

  // TODO: make sure no division by 0
  // const auto dvx = (vdx_high - vdx_low) / static_cast<double>(vx_samples_ - 1);
  // const auto dvy = (vdy_high - vdy_low) / static_cast<double>(vy_samples_ - 1);
  // const auto dw = (wd_high - wd_low) / static_cast<double>(vth_samples_ - 1);

  // std::cout << "Window discretization " << std::endl;
  // std::cout << "dvx: " << dvx << std::endl;
  // std::cout << "dvy: " << dvy << std::endl;
  // std::cout << "dw: " << dw << std::endl;

  // Time parameterization of the reference trajectory
  const double tf = static_cast<double>(xt_ref.n_cols) * delta;

  // if (horizon_ > tf)
  // {
  //   std::cout << "WARNING: Ergodic control horizon is less \
  //                 than the dynamic window horizon causing the search \
  //                 to extrapolate which will lead to uncertainty."
  //             << std::endl;
  // }

  // const vec tvec_ref = arma::linspace(0.0, tf, xt_ref.n_cols);
  //
  // // Time parameterization of the dwa trajectory
  // const vec tvec_dwa = arma::linspace(0.0, horizon_, steps_);
  //
  // // Fit quartic polynomials to the reference trajectory
  // const vec poly_x = polyfit(tvec_ref, xt_ref.row(0), 4);
  // const vec poly_y = polyfit(tvec_ref, xt_ref.row(1), 4);
  // const vec poly_th = polyfit(tvec_ref, xt_ref.row(2), 4);
  //
  // // Interpolate
  // mat xt_param(3, steps_);
  // xt_param.row(0) = polyval(poly_x, tvec_dwa).t();
  // xt_param.row(1) = polyval(poly_y, tvec_dwa).t();
  // xt_param.row(2) = polyval(poly_th, tvec_dwa).t();

  // Search velocity space
  auto min_cost = std::numeric_limits<double>::max();

  bool soln_found = false;

  vec u_opt(3, arma::fill::zeros);
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

        // u.print("u tests");

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

    // cost += (1.0 / dmin);
  }

  // control error;
  const vec cntrl_error = vref - u;
  cost += dot(cntrl_error, cntrl_error);

  // const auto v_trans =
  //     std::sqrt(max_vel_x_ * max_vel_x_ + max_vel_y_ * max_vel_y_);
  // const auto v_sample = std::sqrt(u(0) * u(0) + u(1) * u(1));
  //
  // loss += v_trans - v_sample;

  // loss += std::abs(u(2) - max_rot_vel_);
  // loss += std::abs(u(2) - max_rot_vel);

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

    // index into
    const auto j = static_cast<unsigned int>(std::round((xt_ref.n_cols - 1) * t / tf));

    // std::cout << "t: " << t << " j: " << j << std::endl;


    // if (j > xt_ref.n_cols - 1)
    // {
    //   std::cout << "ERROR" << std::endl;
    // }

    cost += arma::norm(xt_ref(span(0, 1), span(j, j)) - pose.rows(0, 1));
    cost += std::abs(normalize_angle_PI(normalize_angle_PI(xt_ref(2, j)) - pose(2)));

    // cost += (1.0 / dmin);

    t += dt_;
  }

  return false;
}

void DynamicWindow::path(nav_msgs::Path& path, const vec& x, const vec& u,
                         std::string frame) const
{
  path.header.frame_id = frame;
  path.poses.resize(steps_);

  vec pose = x;
  const vec delta = integrate_twist(x, u, dt_);

  for (unsigned int i = 0; i < steps_; i++)
  {
    pose += delta;

    path.poses.at(i).pose.position.x = pose(0);
    path.poses.at(i).pose.position.y = pose(1);

    tf2::Quaternion quat;
    quat.setRPY(0.0, 0.0, normalize_angle_PI(pose(2)));

    path.poses.at(i).pose.orientation.x = quat.x();
    path.poses.at(i).pose.orientation.y = quat.y();
    path.poses.at(i).pose.orientation.z = quat.z();
    path.poses.at(i).pose.orientation.w = quat.w();
  }
}
}  // namespace ergodic_exploration
