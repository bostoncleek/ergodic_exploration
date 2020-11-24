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
  // return -gdx - dbar - fdx.t() * rho;
  return -gdx - fdx.t() * rho;
}

/** @brief Receding horizon ergodic trajectory optimization */
template <class ModelT>
class ErgodicControl
{
public:
  /**
   * @brief Constructor
   * @param model - robot's dynamic model
   * @param dt - time step in integration
   * @param horizon - control horizon
   * @param num_samples - samples drawn from map
   * @param buffer_size - past states in memory
   * @param batch_size - states sampled from memory
   * @param R - positive definite matrix that penalizes controls
   */
  ErgodicControl(const ModelT& model, const Collision& collision, double dt,
                 double horizon, double resolution, double exploration_weight,
                 unsigned int num_basis, unsigned int buffer_size,
                 unsigned int batch_size, const mat& R, const vec& umin, const vec& umax);
  /**
   * @brief Update the control signal
   * @param grid - grid map
   * @param x - current state [x, y, theta]
   * @return first twist in the updated control signal
   */
  vec control(const GridMap& grid, const vec& x);

  void path(nav_msgs::Path& path, std::string frame) const;

  void setTarget(const Target& target);

private:
  void configTarget(const GridMap& grid);

  /**
   * @brief Compose time averaged trajectory fourier coefficients
   * @param xt - trajectory
   * @details xt contains the predicted trajectory + sampled states from memory
   */
  // void trajCoeff(const mat& xt);

  /**
   * @brief Compose the gradient of the ergodic metric
   * @param xt - predicted trajectory
   * @details Updates a matrix containing ergodic metric gradient for each state in xt
   */
  void gradErgodicMetric(const mat& xt);

  /**
   * @brief Update the control signal
   * @param xt - predicted trajectory
   * @param rhot - co-state variable solution
   * @details rhot is assumed to already be sorted from [t0 tf] with t0 at index 0
   */
  void updateControl(const mat& xt, const mat& rhot);

  /**
   * @brief Log loss barrier function to obstacles of the form -log(b-x)
   * @param grid - grid map
   * @param xt - predicted trajectory
   * @details Updates a matrix conatining the log loss barrier function derivatve
   * for each state in xt
   */
  void barrier(const GridMap& grid);

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
  vec umin_;             // lower limif on controls
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
  , ut_(model.action_space, steps_, arma::fill::zeros)
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
  configTarget(grid);

  // Shift columns to the left by 1 and set last column to zeros
  ut_.cols(0, ut_.n_cols - 2) = ut_.cols(1, ut_.n_cols - 1);
  ut_.col(ut_.n_cols - 1).fill(0.0);

  // Forward simulation
  RungeKutta rk4(dt_);
  rk4.solve(traj_, model_, x, ut_, horizon_);
  // RungeKutta45 model_rk45(0.001, 0.01, 1e-3, 1e6);
  // model_rk45.solve(xt_, model_, x, ut_, dt_, horizon_);
  // traj_.print("traj_:");

  // Sample past states
  mat xt_total;
  buffer_.sampleMemory(xt_total, traj_);
  // std::cout << xt_total.n_cols << std::endl;

  // Transform from map frame to fourier frame
  xt_total.row(0) -= map_pos_(0);
  xt_total.row(1) -= map_pos_(1);

  for (unsigned int i = 0; i < xt_total.n_cols; i++)
  {
    if (xt_total(0, i) < 0.0 || xt_total(1, i) < 0.0)
    {
      std::cout << "WARNING: Trajectory is not within fourier domain" << std::endl;
      xt_total.rows(0, 1).print();
    }
  }

  // Extract optimized trajectory in fourier domain
  mat xt = xt_total.cols(xt_total.n_cols - steps_, xt_total.n_cols - 1);
  // std::cout << steps_ << " " << xt.n_cols << std::endl;

  // Update the trajectory fourier coefficients
  basis_.trajCoeff(ck_, xt_total);

  // Gradient of he ergodic measure w.r.t the state
  gradErgodicMetric(xt);
  edx_.rows(0, 1) *= expl_weight_;
  // edx_.print("edx:");

  // barrier(grid);
  // bdx_.print("bdx_:");

  // Backwards pass
  mat rhot;
  rk4.solve(rhot, rhodot_, model_, rhoT_, xt, ut_, edx_, bdx_, horizon_);
  // rhot.print("rho: rk4");

  // RungeKutta45 costate_rk45(0.001, 0.01, 1e-2, 1e6);
  // costate_rk45.solve(rhot, rhodot_, model_, rhoT_, xt_, ut_, edx_, bdx_, dt_,
  // horizon_);
  // rhot.print("rho: rk45");

  const auto max_rho = max(max(rhot, 1));
  const auto min_rho = min(min(rhot, 1));

  // std::cout << "max rho: " << max(max(rhot, 1)) << std::endl;
  // std::cout << "min rho: " << min(min(rhot, 1)) << std::endl;

  if (max_rho > 100.0)
  {
    std::cout << "OVERFLOW max rho: " << max_rho << std::endl;
  }

  if (min_rho < -100.0)
  {
    std::cout << "OVERFLOW min rho: " << min_rho << std::endl;
  }

  // std::cout << "max rho: " << max(max(rhot, 1)) << std::endl;
  // std::cout << "min rho: " << min(min(rhot, 1)) << std::endl;

  // max(rhot, 1).print("max rho");
  // min(rhot, 1).print("min rho");

  updateControl(xt, rhot);
  // ut_.print("Control:");

  // Store current state in frame of map
  buffer_.append(x);

  // xt.row(0) += map_pos_(0);
  // xt.row(1) += map_pos_(1);
  // xt.print("xt");

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

  // size of map
  const auto mx = grid.xmax() - grid.xmin();
  const auto my = grid.ymax() - grid.ymin();

  // Map is the same size
  // if (mx > basis.lx_ || my > basis.ly_)
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

  phi_grid_.resize(2, nx * ny);
  phi_vals_.resize(nx * ny);

  // Construction the grid
  auto col = 0;
  auto y = 0.0;
  for (unsigned int i = 0; i < ny; i++)
  {
    auto x = 0.0;
    for (unsigned int j = 0; j < nx; j++)
    {
      phi_grid_(0, col) = x;
      phi_grid_(1, col) = y;

      // std::cout << x << " " << y << std::endl;
      col++;
      x += resolution_;
    }
    y += resolution_;
  }
  // phi_grid_.t().print("phi_grid");

  // Evaluate each grid cell
  target_.fill(phi_vals_, map_pos_, phi_grid_);
  // phi_vals_.print("phi_vals");

  // spatial fourier coefficients
  basis_.spatialCoeff(phik_, phi_vals_, phi_grid_);
  // phik_.print("phik_");

  std::cout << "done" << std::endl;
}

template <class ModelT>
void ErgodicControl<ModelT>::gradErgodicMetric(const mat& xt)
{
  // ck_.print("ck_");
  // phik_.print("phik_");

  const vec fourier_diff = basis_.lamdak_ % (ck_ - phik_);
  // fourier_diff.print("fourier_diff");

  // TODO: add weight to influence ergodicity
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
  for (unsigned int i = 0; i < xt.n_cols; i++)
  {
    // rhot is already sorted from t0 to tf
    ut_.col(i) = -Rinv_ * model_.fdu(xt.col(i)).t() * rhot.col(i);

    // TODO: apply filter to control signal (savgol filter)
    // ut_.col(i).print("u:");

    for (unsigned int j = 0; j < ut_.n_rows; j++)
    {
      ut_(j, i) = std::clamp(ut_(j, i), umin_(j), umax_(j));
    }

    // ut_(0, i) = std::clamp(ut_(0, i), -0.1, 0.1);
    // ut_(1, i) = std::clamp(ut_(1, i), -0.1, 0.1);
    // ut_(2, i) = std::clamp(ut_(2, i), -0.5, 0.5);

    // ut_(0, i) = std::clamp(ut_(0, i), -1.0, 1.0);
    // ut_(1, i) = std::clamp(ut_(1, i), -1.0, 1.0);
    // ut_(2, i) = std::clamp(ut_(2, i), -2.0, 2.0);

    // if (any(ut_.col(i) > 1.0) || any(ut_.col(i) < -1.0))
    // {
    //   // euclidean norm
    //   ut_.col(i) /= norm(ut_.col(i), 2);
    // }
  }
}

template <class ModelT>
void ErgodicControl<ModelT>::barrier(const GridMap& grid)
{
  const auto weight = 1e12;
  const auto eps = collision_.totalPadding();

  for (unsigned int i = 0; i < steps_; i++)
  {
    vec dbar;
    if (collision_.minDirection(dbar, grid, traj_.col(i)) == CollisionMsg::obstacle)
    {
      // x
      if (std::abs(dbar(0)) < eps)
      {
        dbar(0) *= weight;
      }
      else
      {
        dbar(0) *= 1.0 / dbar(0);
        // std::cout << "dbar x: " << dbar(0) << std::endl;
      }

      // y
      if (std::abs(dbar(1)) < eps)
      {
        dbar(1) *= weight;
      }
      else
      {
        dbar(1) = 1.0 / dbar(1);
        // std::cout << "dbar y: " << dbar(1) << std::endl;
      }
    }

    bdx_(span(0, 1), span(i, i)) = dbar;
  }
}

}  // namespace ergodic_exploration
#endif
