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

#include <ergodic_exploration/buffer.hpp>
#include <ergodic_exploration/grid.hpp>
#include <ergodic_exploration/numerics.hpp>
#include <ergodic_exploration/integrator.hpp>

namespace ergodic_exploration
{
using arma::span;

// TODO: check base 2 or e ???
/**
 * @brief Entropy of a single grid cell
 * @param prob_occu - probability grid cell is occupied represented as a decimal
 * @return entropy
 */
inline double entropy(double p)
{
  // if ( p > 0.0 && p < 1.0)
  // {
  //   std::cout << "p: " << p << std::endl;
  // }

  // Assign zero information gain
  if (almost_equal(0.0, p) || almost_equal(1.0, p) /*|| p < 0.0*/)
  {
    return 0.0;
  }

  else if (p < 0.0)
  {
    return 0.7;
  }

  return -p * std::log(p) - (1.0 - p) * std::log(1.0 - p);
}

/**
 * @brief Time derivatve of the co-state variable
 * @param rho - co-state variable
 * @param kldx - ergodic measure derivative
 * @param dbar - derivatve of barrier function
 * @param fdx - jacobian of the dynamics w.r.t state A = D1(f(x,u))
 * @return  d/dt[rho]
 */
inline vec rhodot(const vec& rho, const vec& kldx, const vec& dbar, const mat& fdx)
{
  return kldx - dbar - fdx.t() * rho;
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
   * @param Sigma - uncertainty in robot's (x,y) position
   */
  ErgodicControl(const ModelT& model, double dt, double horizon, unsigned int num_samples,
                 unsigned int buffer_size, unsigned int batch_size, const mat& R,
                 const mat& Sigma);

  /**
   * @brief Update the control signal
   * @param collision - collision detector
   * @param grid - grid map
   * @param x - current state [x, y, theta]
   * @return first twist in the updated control signal
   */
  vec control(const Collision& collision, const GridMap& grid, const vec& x);

private:
  /**
   * @brief Initialize proposal distribution
   * @param grid - occupancy map
   */
  void initProposal(const GridMap& grid);

  /**
   * @brief Update proposal distribution
   */
  void updateProposal();

  /**
   * @brief Sample proposal distribution
   * @param grid - occupancy map
   */
  void sampleProposal(const GridMap& grid);

  /**
   * @brief Compose time averaged trajectory statistics
   * @param xt - trajectory
   * @details xt contains the predicted trajectory + sampled states from memory
   */
  void trajStat(const mat& xt);

  /**
   * @brief Compose the derivative of the ergodic measure
   * @param xt - predicted trajectory
   * @details Updates a matrix containing ergodic measure derivatve for each state in xt
   */
  void ergMeasDeriv(const mat& xt);

  /**
   * @brief Log loss barrier function to obstacles of the form -log(b-x)
   * @param collision - collision detector
   * @param grid - grid map
   * @param xt - predicted trajectory
   * @details Updates a matrix conatining the log loss barrier function derivatve
   * for each state in xt
   */
  void barrier(const Collision& collision, const GridMap& grid, const mat& xt);

  /**
   * @brief Update the control signal
   * @param xt - predicted trajectory
   * @param rhot - co-state variable solution
   * @details rhot is assumed to already be sorted from [t0 tf] with t0 at index 0
   */
  void updateControl(const mat& xt, const mat& rhot);

private:
  ModelT model_;              // robot's dynamic model
  double dt_;                 // time step in integration
  double horizon_;            // control horizon
  unsigned int num_samples_;  // number of samples drawn from map
  unsigned int steps_;        // number of steps used in integration
  bool proposal_;             // proposal distribution is initialized
  mat Rinv_;                  // inverse of the control weight matrix
  mat Sigmainv_;              // inverse of the robot's (x,y) position uncertainty
  mat ut_;                    // control signal
  mat s_;                     // pair of (x,y) points sampled from map
  mat edx_;                   // ergodic measure derivatives
  mat bdx_;                   // log loss barrier derivatives
  mat cov_;                   // target distribution covariance
  vec mu_;                    // target distribution mean
  vec rhoT_;                  // co-state terminal condition
  vec w_;                     // importance sample weights
  vec p_;                     // occupancy map entropy vector
  vec q_;                     // trajectory statistics vector
  RungeKutta rk4_;            // integrator
  ReplayBuffer buffer_;       // store past states
  CoStateFunc rhodot_;        // co-state function
};

template <class ModelT>
ErgodicControl<ModelT>::ErgodicControl(const ModelT& model, double dt, double horizon,
                                       unsigned int num_samples, unsigned int buffer_size,
                                       unsigned int batch_size, const mat& R,
                                       const mat& Sigma)
  : model_(model)
  , dt_(dt)
  , horizon_(horizon)
  , num_samples_(num_samples)
  , steps_(static_cast<unsigned int>(horizon / dt))
  , proposal_(false)
  , Rinv_(inv(R))
  , Sigmainv_(inv(Sigma))
  , ut_(model.action_space, steps_, arma::fill::zeros)
  , s_(2, num_samples)  // exploration space is set to (x,y)
  , edx_(model.state_space, steps_)
  , bdx_(model.state_space, steps_, arma::fill::zeros)  // set heading column to zeros
  , cov_(2, 2)
  , mu_(2)
  , rhoT_(model.state_space, arma::fill::zeros)
  , w_(num_samples)
  , p_(num_samples)
  , q_(num_samples)
  , rk4_(dt)
  , buffer_(buffer_size, batch_size)
{
  if (steps_ == 1)
  {
    throw std::invalid_argument("Need at least two steps in forward simulation. \
                                 Increase the horizon or decrease the time step.");
  }

  arma::arma_rng::set_seed_random();
  // arma::arma_rng::set_seed(0);

  rhodot_ = std::bind(rhodot, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, std::placeholders::_4);
}

template <class ModelT>
vec ErgodicControl<ModelT>::control(const Collision& collision, const GridMap& grid,
                                    const vec& x)
{
  // initialize the proposal distribution
  // if (!proposal_)
  // {
  //   initProposal(grid);
  // }

  // Shift columns to the left by 1 and set last column to zeros
  ut_.cols(0, ut_.n_cols - 2) = ut_.cols(1, ut_.n_cols - 1);
  ut_.col(ut_.n_cols - 1).fill(0.0);

  // Forward simulation
  mat xt;
  rk4_.solve(xt, model_, x, ut_, horizon_);
  // xt.print("xt:");

  // Sample (x,y) space
  sampleProposal(grid);
  // p_.print("p:");

  // Sample past states
  mat xt_total;
  buffer_.sampleMemory(xt_total, xt);

  // Time average trajectory statistics using memory and predicted trajectory
  trajStat(xt_total);
  // q_.print("q:");

  // Derivative of the ergodic measure w.r.t the state
  ergMeasDeriv(xt);
  // edx_.print("edx:");

  // Derivative of the barrier function
  barrier(collision, grid, xt);

  // Backwards pass
  mat rhot;
  rk4_.solve(rhot, rhodot_, model_, rhoT_, xt, ut_, edx_, bdx_, horizon_);
  // rhot.print("rho:");

  // Update controls
  updateControl(xt, rhot);
  // ut_.print("Control:");

  // add current state to memory
  buffer_.append(x);

  return ut_.col(0);
}

template <class ModelT>
void ErgodicControl<ModelT>::initProposal(const GridMap& grid)
{
  proposal_ = true;

  // Uniformly generate sample [0 1]
  s_.randu(2, num_samples_);

  // Scale based on grid domain
  s_.row(0) *= (grid.xmax() - grid.xmin()) + grid.xmin();
  s_.row(1) *= (grid.ymax() - grid.ymin()) + grid.ymin();

  for (unsigned int i = 0; i < num_samples_; i++)
  {
    w_(i) = grid.getCell(s_(0, i), s_(1, i));
  }

  // Normalize the sampled values
  const auto total = sum(w_);
  if (!almost_equal(0.0, total))
  {
    w_ /= total;
  }
  // std::cout << "Sum: " << sum << std::endl;

  updateProposal();
}

template <class ModelT>
void ErgodicControl<ModelT>::updateProposal()
{
  vec mu_prev = mu_;

  // Update mean
  mu_ = s_ * w_;
  // mu_.print("mean");

  // // Update covariance
  // mat diff1 = s_;
  // diff1.each_col() -= mu_;
  //
  // mat diff2 = diff1;
  //
  // diff1.each_row() %= p_.t();
  //
  // cov_ = diff1 * diff2.t();

  // TODO: does the cov update use mu at i not i+1?
  vec diff(2);
  for (unsigned int i = 0; i < num_samples_; i++)
  {
    diff = s_(i) - mu_prev;
    cov_ += w_(i) * diff * diff.t();
  }

  // cov_.print("covariance");

  // std::cout << "EES: " << sum(square(w_)) << std::endl;
}

template <class ModelT>
void ErgodicControl<ModelT>::sampleProposal(const GridMap& grid)
{
  // vec sample(2);
  // unsigned int i = 0;
  // while (i < num_samples_)
  // {
  //   if (!mvnrnd(sample, mu_, cov_))
  //   {
  //     std::cout << "FAILURE! Unable to sample proposal distribution. " << std::endl;
  //   }
  //
  //   if ((sample(0) > grid.xmax()) || (sample(0) < grid.xmin()) ||
  //       (sample(1) > grid.ymax()) || (sample(1) < grid.ymin()))
  //   {
  //     // std::cout << "FAILURED! Sample outside if grid bounds. " << std::endl;
  //     // s_.col(i).print("sample:");
  //     continue;
  //   }
  //
  //   s_.col(i) = sample;
  //   p_(i) = entropy(grid.getCell(s_(0, i), s_(1, i)));
  //   i++;
  // }
  //
  // // Normalize the sampled values
  // const auto total = sum(p_);
  // if (!almost_equal(0.0, total))
  // {
  //   p_ /= total;
  // }
  //
  // // std::cout << "EES: " << sum(square(p_)) << std::endl;
  //
  // updateProposal();

  // Uniformly generate sample [0 1]
  s_.randu(2, num_samples_);

  // Scale based on grid domain
  s_.row(0) = s_.row(0) * (grid.xmax() - grid.xmin()) + grid.xmin();
  s_.row(1) = s_.row(1) * (grid.ymax() - grid.ymin()) + grid.ymin();

  for (unsigned int i = 0; i < num_samples_; i++)
  {
    p_(i) = entropy(grid.getCell(s_(0, i), s_(1, i)));
  }

  // Normalize the sampled values
  const auto sum = accu(p_);
  if (!almost_equal(0.0, sum))
  {
    p_ /= sum;
  }
}

template <class ModelT>
void ErgodicControl<ModelT>::trajStat(const mat& xt)
{
  for (unsigned int i = 0; i < num_samples_; i++)
  {
    auto qval = 0.0;
    for (unsigned int j = 0; j < xt.n_cols; j++)
    {
      const vec diff = s_.col(i) - xt(span(0, 1), span(j, j));
      qval += std::exp(-0.5 * dot(diff.t() * Sigmainv_, diff));
      // std::cout << std::exp(-0.5 * dot(diff.t() * Sigmainv_, diff)) << std::endl;
    }

    if (almost_equal(0.0, qval))
    {
      std::cout << "WARNING: trajectory statistics are zero" << std::endl;
      qval = 1e-8;
    }

    q_(i) = qval;
  }

  // Normalize the trajectory statistics
  const auto total = sum(q_);
  if (!almost_equal(0.0, total))
  {
    q_ /= total;
    // std::cout << "sum q: " << sum << std::endl;
  }
}

template <class ModelT>
void ErgodicControl<ModelT>::ergMeasDeriv(const mat& xt)
{
  for (unsigned int i = 0; i < steps_; i++)
  {
    vec kldx(3, arma::fill::zeros);
    for (unsigned int j = 0; j < num_samples_; j++)
    {
      const vec diff = s_.col(j) - xt(span(0, 1), span(i, i));
      const auto qval = std::exp(-0.5 * dot(diff.t() * Sigmainv_, diff));
      kldx.rows(0, 1) += (p_(j) / q_(j)) * qval * (Sigmainv_ * diff);
    }

    edx_.col(i) = kldx;
  }
}

template <class ModelT>
void ErgodicControl<ModelT>::barrier(const Collision& collision, const GridMap& grid,
                                     const mat& xt)
{
  // TODO: add eps param as obstacle padding
  const auto weight = 1e12;
  const auto eps = collision.totalPadding();

  for (unsigned int i = 0; i < steps_; i++)
  {
    vec dbar;
    if (collision.minDirection(dbar, grid, xt.col(i)) == CollisionMsg::obstacle)
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

template <class ModelT>
void ErgodicControl<ModelT>::updateControl(const mat& xt, const mat& rhot)
{
  for (unsigned int i = 0; i < xt.n_cols; i++)
  {
    // rhot is already sorted from t0 to tf
    ut_.col(i) = -Rinv_ * model_.fdu(xt.col(i)).t() * rhot.col(i);

    // TODO: add control limits as a parameter
    // TODO: apply filter to control signal (savgol filter)

    // ut_.col(i).print();

    // see armadillo clamp()
    ut_(0, i) = std::clamp(ut_(0, i), -0.1, 0.1);
    ut_(1, i) = std::clamp(ut_(1, i), -0.1, 0.1);
    ut_(2, i) = std::clamp(ut_(2, i), -0.5, 0.5);

    // ut_(0, i) = std::clamp(ut_(0, i), -1.0, 1.0);
    // ut_(1, i) = std::clamp(ut_(1, i), -1.0, 1.0);
    // ut_(2, i) = std::clamp(ut_(2, i), -2.0, 2.0);

    // ut_(0, i) = std::clamp(ut_(0, i), -1.0, 1.0);
    // ut_(1, i) = std::clamp(ut_(1, i), -1.0, 1.0);
    // ut_(2, i) = std::clamp(ut_(2, i), -0.5, 0.5);

    // ut_(0,i) = std::clamp(ut_(0,i), -1.0, 1.0);
    // ut_(1,i) = std::clamp(ut_(1,i), -1.0, 1.0);
    // ut_(2,i) = std::clamp(ut_(2,i), -1.0, 1.0);
    // ut_(3,i) = std::clamp(ut_(3,i), -1.0, 1.0);

    // if (any(ut_.col(i) > 1.0) || any(ut_.col(i) < -1.0))
    // {
    //   // euclidean norm
    //   ut_.col(i) /= norm(ut_.col(i), 2);
    // }
  }
}
}  // namespace ergodic_exploration
#endif
