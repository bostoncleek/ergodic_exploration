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
 * @file ergodic_control.hpp
 * @aut_hor Boston Cleek
 * @date 5 Nov 2020
 * @brief Ergodic control strategy for exploration
 */
#ifndef ERGODIC_CONTROL_HPP
#define ERGODIC_CONTROL_HPP

#include <cmath>
#include <stdexcept>
#include <algorithm>

#include <nav_msgs/Path.h>
#include <tf2/LinearMath/Quaternion.h>

#include <ergodic_exploration/basis.hpp>
#include <ergodic_exploration/buffer.hpp>
#include <ergodic_exploration/integrator.hpp>
#include <ergodic_exploration/target.hpp>

namespace ergodic_exploration
{
using arma::span;

/**
 * @brief Time derivatve of the co-state variable
 * @param rho - co-state variable
 * @param gdx - gradient of the ergodic metric
 * @param dbar - derivatve of barrier function
 * @param fdx - jacobian of the dynamics w.r.t state A = D1(f(x,u))
 * @return  d/dt[rho]
 */
inline vec rhodot(const vec& rho, const vec& gdx, const vec& dbar, const mat& fdx)
{
  return -gdx - dbar - fdx.t() * rho;
  // return -gdx - fdx.t() * rho;
}

/** @brief Receding horizon ergodic trajectory optimization */
template <class ModelT>
class ErgodicControl
{
public:
  /**
   * @brief Constructor
   * @param model - robot's dynamic model
   * @param collision - collision detector
   * @param dt - time step in integration
   * @param horizon - control horizon
   * @param resolution - discretization of target grid
   * @param num_basis - number of cosine basis functions per dimension
   * @param buffer_size - number of past states in memory
   * @param batch_size - number of states sampled from memory
   * @param Rinv - inverse of a positive definite matrix that penalizes controls
   * @param umin - body twist lower limits
   * @param umax - body twist upper limits
   */
  ErgodicControl(const ModelT& model, const Collision& collision, double dt,
                 double horizon, double resolution, double exploration_weight,
                 unsigned int num_basis, unsigned int buffer_size,
                 unsigned int batch_size, const mat& Rinv, const vec& umin,
                 const vec& umax);

  /**
   * @brief Update the control signal
   * @param grid - grid map
   * @param x - current state [x, y, theta]
   * @return first twist in the updated control signal [vx, vy, w]
   */
  vec control(const GridMap& grid, const vec& x);

  /**
   * @brief Return optimized trajectory
   */
  mat optTraj() const;

  /**
   * @brief Return optimized trajectory
   */
  nav_msgs::Path path(const std::string& map_frame_id) const;

  /**
   * @brief Add the robot's state to memory
   * @param x - current state [x, y, theta]
   */
  void addStateMemory(const vec& x);

  /** @brief return time step */
  double timeStep() const;

  /**
   * @brief Set the target distribution
   * @param target - target distribution
   */
  void setTarget(const Target& target);

  /**
   * @brief Compose target distribution grid and spatial coefficients
   * @param grid - grid map
   * @details Updates the spatial coefficients if the map is larger than before
   */
  void configTarget(const GridMap& grid);

private:
  /**
   * @brief Compose the gradient of the ergodic metric
   * @param ck - trajectory fourier coefficients
   * @param xt - forward simulated trajectory in fourier domain
   * @return ergodic measure gradients
   * @details Updates a matrix containing ergodic metric gradient for each state in xt
   */
  mat gradErgodicMetric(const vec& ck, const mat& xt);

  /**
   * @brief Update the control signal
   * @param xt - forward simulated trajectory in fourier domain
   * @param rhot - co-state variable solution
   * @details rhot is assumed to already be sorted from [t0 tf] with t0 at index 0
   */
  void updateControl(const mat& xt, const mat& rhot);

  /**
   * @brief Gradient of the barrier function
   * @param xt - forward simulated trajectory in fourier domain
   * @return barrier funcion gradients
   * @details Updates a matrix conatining the barrier function derivatve
   * for each state in xt. The barrier function tries to keep the robot's
   * predicted trajectory within the fourier domain.
   */
  mat gradBarrier(const mat& xt);

private:
  ModelT model_;         // robot's dynamic model
  Collision collision_;  // collision detection
  double dt_;            // time step in integration
  double horizon_;       // control horizon
  double resolution_;    // target grid resolution in meters
  double expl_weight_;   // ergodic exploration weight
  unsigned int steps_;   // number of steps used in integration
  mat Rinv_;             // inverse of the control weight matrix
  mat ut_;               // control signal
  vec phik_;             // target distribution fourier coefficients
  vec rhoT_;             // co-state terminal condition
  vec map_pos_;          // position of map
  vec umin_;             // lower limit on controls
  vec umax_;             // upper limit on controls
  vec pose_;             // current state
  Basis basis_;          // fourier basis
  ReplayBuffer buffer_;  // store past states in frame of occupancy map
  RungeKutta rk4_;       // Runge-Kutta integrator
  CoStateFunc rhodot_;   // co-state function
  Target target_;        // target distribution
};

template <class ModelT>
ErgodicControl<ModelT>::ErgodicControl(const ModelT& model, const Collision& collision,
                                       double dt, double horizon, double resolution,
                                       double exploration_weight, unsigned int num_basis,
                                       unsigned int buffer_size, unsigned int batch_size,
                                       const mat& Rinv, const vec& umin, const vec& umax)
  : model_(model)
  , collision_(collision)
  , dt_(dt)
  , horizon_(horizon)
  , resolution_(resolution)
  , expl_weight_(exploration_weight)
  , steps_(static_cast<unsigned int>(std::abs(horizon / dt)))
  , Rinv_(Rinv)
  , ut_(3, steps_, arma::fill::zeros)  // body twist Vb = [vx, vy, w]
  , phik_(num_basis * num_basis)
  , rhoT_(model.state_space, arma::fill::zeros)
  , map_pos_(2, arma::fill::zeros)
  , umin_(umin)
  , umax_(umax)
  , pose_(model.state_space, arma::fill::zeros)
  , basis_(0.0, 0.0, num_basis)  // fourier domain init to 0
  , buffer_(buffer_size, batch_size)
  , rk4_(dt)
{
  if (steps_ == 1)
  {
    throw std::invalid_argument("Need at least two steps in forward simulation. \
                                 Increase the horizon or decrease the time step.");
  }

  // arma::arma_rng::set_seed_random();

  rhodot_ = std::bind(rhodot, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, std::placeholders::_4);
}

template <class ModelT>
vec ErgodicControl<ModelT>::control(const GridMap& grid, const vec& x)
{
  pose_ = x;

  // Update target grid if needed
  configTarget(grid);

  // Shift columns to the left by 1 and set last column to zeros
  ut_.cols(0, ut_.n_cols - 2) = ut_.cols(1, ut_.n_cols - 1);
  ut_.col(ut_.n_cols - 1).fill(0.0);

  // Forward simulation
  const mat traj = rk4_.solve(model_, pose_, ut_, horizon_);

  // Sample past states
  mat xt_total = buffer_.sampleMemory(traj);

  // Transform from map frame to fourier frame
  xt_total.row(0) -= map_pos_(0);
  xt_total.row(1) -= map_pos_(1);

  //////////////////////////////////////////////////////////////////////////////
  // DEBUG
  // for (unsigned int i = 0; i < xt_total.n_cols; i++)
  // {
  //   if (xt_total(2, i) < -PI|| xt_total(2, i) > PI)
  //   {
  //       std::cout << "WARNING: heading is not normalized" << std::endl;
  //   }
  //
  //   // if (xt_total(0, i) < 0.0 || xt_total(1, i) < 0.0)
  //   // {
  //   //   std::cout << "WARNING: Trajectory is not within fourier domain" << std::endl;
  //   //   // xt_total.rows(0, 1).print();
  //   // }
  // }
  //////////////////////////////////////////////////////////////////////////////

  // Extract optimized trajectory in fourier domain
  const mat xt = xt_total.cols(xt_total.n_cols - steps_, xt_total.n_cols - 1);

  // Update the trajectory fourier coefficients
  const vec ck = basis_.trajCoeff(xt_total);

  // Gradient of the ergodic measure w.r.t the state
  const mat edx = gradErgodicMetric(ck, xt);

  // Gradient of the barrier function w.r.t the state
  const mat bdx = gradBarrier(xt);
  // bdx_.print("bdx:");

  // Backwards pass
  const mat rhot = rk4_.solve(rhodot_, model_, rhoT_, xt, ut_, edx, bdx, horizon_);

  //////////////////////////////////////////////////////////////////////////////
  // DEBUG
  // const auto max_rho = max(max(rhot, 1));
  // const auto min_rho = min(min(rhot, 1));
  //
  // std::cout << "max rho: " << max(max(rhot, 1)) << std::endl;
  // std::cout << "min rho: " << min(min(rhot, 1)) << std::endl;

  // if (max_rho > 100.0)
  // {
  //   std::cout << "OVERFLOW max rho: " << max_rho << std::endl;
  // }
  //
  // if (min_rho < -100.0)
  // {
  //   std::cout << "OVERFLOW min rho: " << min_rho << std::endl;
  // }
  //
  // // if (any(abs(rhot.row(0)) < 1.0e-12) || any(abs(rhot.row(1)) < 1.0e-12) ||
  // // any(abs(rhot.row(2)) < 1.0e-12))
  // // {
  // //   std::cout << "UNDERFLOW: gradient is zero" << std::endl;
  // // }
  //////////////////////////////////////////////////////////////////////////////

  // Update control signal
  updateControl(xt, rhot);

  // Store current state in frame of map
  // buffer_.append(x);

  return ut_.col(0);
}

template <class ModelT>
mat ErgodicControl<ModelT>::optTraj() const
{
  return rk4_.solve(model_, pose_, ut_, horizon_);
}

template <class ModelT>
nav_msgs::Path ErgodicControl<ModelT>::path(const std::string& map_frame_id) const
{
  nav_msgs::Path path;
  path.header.frame_id = map_frame_id;
  path.poses.resize(steps_);

  const mat opt_traj = rk4_.solve(model_, pose_, ut_, horizon_);
  for (unsigned int i = 0; i < opt_traj.n_cols; i++)
  {
    path.poses.at(i).pose.position.x = opt_traj(0, i);
    path.poses.at(i).pose.position.y = opt_traj(1, i);

    tf2::Quaternion quat;
    quat.setRPY(0.0, 0.0, normalize_angle_PI(opt_traj(2, i)));

    path.poses.at(i).pose.orientation.x = quat.x();
    path.poses.at(i).pose.orientation.y = quat.y();
    path.poses.at(i).pose.orientation.z = quat.z();
    path.poses.at(i).pose.orientation.w = quat.w();
  }

  return path;
}

template <class ModelT>
void ErgodicControl<ModelT>::addStateMemory(const vec& x)
{
  buffer_.append(x);
}

template <class ModelT>
double ErgodicControl<ModelT>::timeStep() const
{
  return dt_;
}

template <class ModelT>
void ErgodicControl<ModelT>::setTarget(const Target& target)
{
  target_ = target;
}

template <class ModelT>
void ErgodicControl<ModelT>::configTarget(const GridMap& grid)
{
  // translation from map to fourier domain
  map_pos_(0) = grid.xmin();
  map_pos_(1) = grid.ymin();

  // size of map
  const auto mx = grid.xmax() - grid.xmin();
  const auto my = grid.ymax() - grid.ymin();

  // Map is the same size
  if (almost_equal(mx, basis_.lx_) && almost_equal(my, basis_.ly_))
  {
    return;
  }

  std::cout << "Updating target distribution grid... ";

  // update phi_grid if map has grown
  basis_.lx_ = mx;
  basis_.ly_ = my;

  // Add 1 to include the boundary
  // Use resolution that is specified for the target grid
  const auto nx = axis_length(0.0, basis_.lx_, resolution_) + 1;
  const auto ny = axis_length(0.0, basis_.ly_, resolution_) + 1;

  // Construct the grid
  mat phi_grid(2, nx * ny);

  unsigned int col = 0;
  auto y = 0.0;
  for (unsigned int i = 0; i < ny; i++)
  {
    auto x = 0.0;
    for (unsigned int j = 0; j < nx; j++)
    {
      phi_grid(0, col) = x;
      phi_grid(1, col) = y;
      // std::cout << x << " " << y << std::endl;

      col++;
      x += resolution_;
    }
    y += resolution_;
  }

  // Evaluate each grid cell
  const vec phi_vals = target_.fill(map_pos_, phi_grid);

  phik_ = basis_.spatialCoeff(phi_vals, phi_grid);

  std::cout << "done" << std::endl;
}

template <class ModelT>
mat ErgodicControl<ModelT>::gradErgodicMetric(const vec& ck, const mat& xt)
{
  // element wise multiplication
  const vec fourier_diff = basis_.lamdak_ % (ck - phik_);

  mat edx(model_.state_space, steps_);
  for (unsigned int i = 0; i < steps_; i++)
  {
    edx(span(0, 1), span(i, i)) = basis_.gradFourierBasis(xt.col(i)) * fourier_diff;

    // set heading to zero
    edx(2, i) = 0.0;
  }

  edx.rows(0, 1) *= expl_weight_;
  // edx_ *= 1.0 / steps_ * expl_weight_;
  return edx;
}

template <class ModelT>
void ErgodicControl<ModelT>::updateControl(const mat& xt, const mat& rhot)
{
  for (unsigned int i = 0; i < xt.n_cols; i++)
  {
    // rhot is already sorted from t0 to tf
    ut_.col(i) = -Rinv_ * model_.fdu(xt.col(i)).t() * rhot.col(i);
    // ut_.col(i).print("u:");

    ut_(0, i) = std::clamp(ut_(0, i), umin_(0), umax_(0));
    ut_(1, i) = std::clamp(ut_(1, i), umin_(1), umax_(1));
    ut_(2, i) = std::clamp(ut_(2, i), umin_(2), umax_(2));
  }
}

template <class ModelT>
mat ErgodicControl<ModelT>::gradBarrier(const mat& xt)
{
  // TODO: load from param sever
  const auto weight = 25.0;
  const auto eps = 0.05;

  // set heading column to zeros
  mat bdx(model_.state_space, steps_, arma::fill::zeros);
  for (unsigned int i = 0; i < steps_; i++)
  {
    bdx(0, i) += 2.0 * (xt(0, i) > basis_.lx_ - eps) * (xt(0, i) - (basis_.lx_ - eps));
    bdx(1, i) += 2.0 * (xt(1, i) > basis_.ly_ - eps) * (xt(1, i) - (basis_.ly_ - eps));

    bdx(0, i) += 2.0 * (xt(0, i) < eps) * (xt(0, i) - eps);
    bdx(1, i) += 2.0 * (xt(1, i) < eps) * (xt(1, i) - eps);
  }

  bdx.rows(0, 1) *= weight;

  return bdx;
}
}  // namespace ergodic_exploration
#endif
