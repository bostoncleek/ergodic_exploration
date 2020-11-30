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
 * @brief Determine if control will cause a collision
 * @param collision - collision detector
 * @param grid - grid map
 * @param x0 - initial state
 * @param u - twist [vx, vy, w]
 * @param dt - time step in integration
 * @param horizon - length of integration
 * @return true if the control is collision free
 * @details The control is assumed to be constant and a twist is
 * integrated for a fixed amout of time
 */
inline bool validate_control(const Collision& collision, const GridMap& grid,
                             const vec& x0, const vec& u, double dt, double horizon)
{
  vec x = x0;
  const vec delta = integrate_twist(x, u, dt);
  const auto steps = static_cast<unsigned int>(std::abs(horizon / dt));

  for (unsigned int i = 0; i < steps; i++)
  {
    x += delta;
    if (collision.collisionCheck(grid, x))
    {
      return false;
    }
  }

  return true;
}

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
   * @param R - positive definite matrix that penalizes controls
   * @param umin - body twist lower limits
   * @param umax - body twist upper limits
   */
  ErgodicControl(const ModelT& model, const Collision& collision, double dt,
                 double horizon, double resolution, double exploration_weight,
                 unsigned int num_basis, unsigned int buffer_size,
                 unsigned int batch_size, const mat& R, const vec& umin, const vec& umax);
  /**
   * @brief Update the control signal
   * @param grid - grid map
   * @param x - current state [x, y, theta]
   * @return first twist in the updated control signal [vx, vy, w]
   */
  vec control(const GridMap& grid, const vec& x);

  /**
   * @brief Previous optimized trajectory
   * @param path - trajectory
   * @param frame - trajectory frame
   */
  void path(nav_msgs::Path& path, std::string frame) const;

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
   * @param xt - forward simulated trajectory in fourier domain
   * @details Updates a matrix containing ergodic metric gradient for each state in xt
   */
  void gradErgodicMetric(const mat& xt);

  /**
   * @brief Update the control signal
   * @param xt - forward simulated trajectory in fourier domain
   * @param rhot - co-state variable solution
   * @details rhot is assumed to already be sorted from [t0 tf] with t0 at index 0
   */
  void updateControl(const mat& xt, const mat& rhot);

  /**
   * @brief Gradient of the barrier function to obstacles of the form (x-b)*(x-b)
   * @param xt - forward simulated trajectory in fourier domain
   * @details Updates a matrix conatining the barrier function derivatve
   * for each state in xt
   */
  void barrier(const mat& xt);

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
  mat edx_;              // ergodic measure derivatives
  mat bdx_;              // log loss barrier derivatives
  mat traj_;             // current trajectory
  mat phi_grid_;         // target grid
  vec phi_vals_;         // target values
  vec phik_;             // target distribution fourier coefficients
  vec ck_;               // trajectory fourier coefficients
  vec rhoT_;             // co-state terminal condition
  vec map_pos_;          // position of map
  vec umin_;             // lower limit on controls
  vec umax_;             // upper limit on controls
  Basis basis_;          // fourier basis
  ReplayBuffer buffer_;  // store past states in frame of occupancy map
  CoStateFunc rhodot_;   // co-state function
  Target target_;        // target distribution
};

template <class ModelT>
ErgodicControl<ModelT>::ErgodicControl(const ModelT& model, const Collision& collision,
                                       double dt, double horizon, double resolution,
                                       double exploration_weight, unsigned int num_basis,
                                       unsigned int buffer_size, unsigned int batch_size,
                                       const mat& R, const vec& umin, const vec& umax)
  : model_(model)
  , collision_(collision)
  , dt_(dt)
  , horizon_(horizon)
  , resolution_(resolution)
  , expl_weight_(exploration_weight)
  , steps_(static_cast<unsigned int>(std::abs(horizon / dt)))
  , Rinv_(inv(R))
  , ut_(3, steps_, arma::fill::zeros)                   // body twist Vb = [vx, vy, w]
  , edx_(model.state_space, steps_, arma::fill::zeros)  // set heading column to zeros
  , bdx_(model.state_space, steps_, arma::fill::zeros)  // set heading column to zeros
  , traj_(model.state_space, steps_)
  , phik_(num_basis * num_basis)
  , ck_(num_basis * num_basis)
  , rhoT_(model.state_space, arma::fill::zeros)
  , map_pos_(2, arma::fill::zeros)
  , umin_(umin)
  , umax_(umax)
  , basis_(0.0, 0.0, num_basis)  // fourier domain init to 0
  , buffer_(buffer_size, batch_size)
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
  // Update target grid if needed
  // configTarget(grid);

  // Shift columns to the left by 1 and set last column to zeros
  ut_.cols(0, ut_.n_cols - 2) = ut_.cols(1, ut_.n_cols - 1);
  ut_.col(ut_.n_cols - 1).fill(0.0);

  // Forward simulation
  RungeKutta rk4(dt_);
  rk4.solve(traj_, model_, x, ut_, horizon_);

  // Sample past states
  mat xt_total;
  buffer_.sampleMemory(xt_total, traj_);

  // Transform from map frame to fourier frame
  xt_total.row(0) -= map_pos_(0);
  xt_total.row(1) -= map_pos_(1);

  //////////////////////////////////////////////////////////////////////////////
  // DEBUG
  // for (unsigned int i = 0; i < xt_total.n_cols; i++)
  // {
  //   if (xt_total(0, i) < 0.0 || xt_total(1, i) < 0.0)
  //   {
  //     std::cout << "WARNING: Trajectory is not within fourier domain" << std::endl;
  //     // xt_total.rows(0, 1).print();
  //   }
  // }
  //////////////////////////////////////////////////////////////////////////////

  // Extract optimized trajectory in fourier domain
  mat xt = xt_total.cols(xt_total.n_cols - steps_, xt_total.n_cols - 1);

  // Update the trajectory fourier coefficients
  basis_.trajCoeff(ck_, xt_total);

  // Gradient of he ergodic measure w.r.t the state
  gradErgodicMetric(xt);
  edx_.rows(0, 1) *= expl_weight_;

  // Barrier function
  barrier(xt);
  // bdx_.print("bdx_:");

  // Backwards pass
  mat rhot;
  rk4.solve(rhot, rhodot_, model_, rhoT_, xt, ut_, edx_, bdx_, horizon_);

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
  buffer_.append(x);

  return ut_.col(0);
}

template <class ModelT>
void ErgodicControl<ModelT>::path(nav_msgs::Path& path, std::string frame) const
{
  path.header.frame_id = frame;
  path.poses.resize(traj_.n_cols);
  for (unsigned int i = 0; i < traj_.n_cols; i++)
  {
    path.poses.at(i).pose.position.x = traj_(0, i);
    path.poses.at(i).pose.position.y = traj_(1, i);

    tf2::Quaternion quat;
    quat.setRPY(0.0, 0.0, traj_(2, i));

    path.poses.at(i).pose.orientation.x = quat.x();
    path.poses.at(i).pose.orientation.y = quat.y();
    path.poses.at(i).pose.orientation.z = quat.z();
    path.poses.at(i).pose.orientation.w = quat.w();
  }
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

  std::cout << "Updating target distribution grid... ";

  // size of map
  basis_.lx_ = grid.xmax() - grid.xmin();
  basis_.ly_ = grid.ymax() - grid.ymin();

  // Add 1 to include the boundary
  // Use resolution that is specified for the target grid
  const auto nx = axis_length(0.0, basis_.lx_, resolution_) + 1;
  const auto ny = axis_length(0.0, basis_.ly_, resolution_) + 1;

  phi_grid_.resize(2, nx * ny);
  phi_vals_.resize(nx * ny);

  const GridData mi_data = grid.gridData();

  // Construct the grid
  unsigned int col = 0;
  auto y = 0.0;
  for (unsigned int i = 0; i < ny; i++)
  {
    auto x = 0.0;
    for (unsigned int j = 0; j < nx; j++)
    {
      phi_grid_(0, col) = x;
      phi_grid_(1, col) = y;
      // std::cout << x << " " << y << std::endl;

      const std::vector<unsigned int> gidx =
          grid.world2Grid(x + grid.xmin(), y + grid.ymin());

      if (grid.gridBounds(gidx.at(0), gidx.at(1)))
      {
        // std::cout << "HERE !!!!!!!!!!!!" << std::endl;

        const auto idx = grid.grid2RowMajor(gidx.at(0), gidx.at(1));

        phi_vals_(col) = static_cast<double>(mi_data.at(idx));
        // std::cout << phi_vals_(col) << std::endl;
      }

      col++;
      x += resolution_;
    }
    y += resolution_;
  }
  // phi_grid_.t().print("phi_grid");

  // spatial fourier coefficients
  phi_vals_ /= sum(phi_vals_);
  // std::cout << "sum phi vals: " << sum(phi_vals_) << std::endl;

  basis_.spatialCoeff(phik_, phi_vals_, phi_grid_);
  // phik_.print("phik_");

  std::cout << "done" << std::endl;
}

template <class ModelT>
void ErgodicControl<ModelT>::gradErgodicMetric(const mat& xt)
{
  const vec fourier_diff = basis_.lamdak_ % (ck_ - phik_);

  mat dfk;
  for (unsigned int i = 0; i < steps_; i++)
  {
    basis_.gradFourierBasis(dfk, xt.col(i));
    edx_(span(0, 1), span(i, i)) = dfk * fourier_diff;
  }

  // edx_ *= 1.0 / steps_;
}

template <class ModelT>
void ErgodicControl<ModelT>::updateControl(const mat& xt, const mat& rhot)
{
  // TODO: apply filter to control signal (savgol filter)
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
void ErgodicControl<ModelT>::barrier(const mat& xt)
{
  // TODO: load from param sever
  const auto weight = 25.0;
  const auto eps = 0.05;

  for (unsigned int i = 0; i < steps_; i++)
  {
    vec dbar(2, arma::fill::zeros);

    dbar(0) += 2.0 * (xt(0, i) > basis_.lx_ - eps) * (xt(0, i) - (basis_.lx_ - eps));
    dbar(1) += 2.0 * (xt(1, i) > basis_.ly_ - eps) * (xt(1, i) - (basis_.ly_ - eps));

    dbar(0) += 2.0 * (xt(0, i) < eps) * (xt(0, i) - eps);
    dbar(1) += 2.0 * (xt(1, i) < eps) * (xt(1, i) - eps);

    // dbar.print("dbar");

    bdx_(span(0, 1), span(i, i)) = weight * dbar;
  }
}

}  // namespace ergodic_exploration
#endif
