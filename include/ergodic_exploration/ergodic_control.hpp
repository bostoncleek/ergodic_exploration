/**
 * @file ergodic_control.hpp
 * @aut_hor Boston Cleek
 * @date 5 Nov 2020
 * @brief Ergodic control strategy for exploration
 */

#pragma once

#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <algorithm>

#include <ergodic_exploration/grid.hpp>
#include <ergodic_exploration/numerics.hpp>
#include <ergodic_exploration/integrator.hpp>

namespace ergodic_exploration
{
using arma::distr_param;
using arma::ivec;
using arma::randi;
using arma::span;

// TODO: check base 2 or e ???
/**
 * @brief Entropy of a single grid cell
 * @param prob_occu - probability grid cell is occupied represented as a decimal
 * @return entropy
 */
inline double entropy(double p)
{
  // Assign zero information gain
  if (almost_equal(0.0, p) || almost_equal(1.0, p) || p < 0.0)
  {
    return 0.0;
  }

  return -p * std::log(p) - (1.0 - p) * std::log(1.0 - p);
}

/**
 * @brief Time derivatve of the co-state variable
 * @param rho - co-state variable
 * @param kldx - ergodic measure derivative
 * @param fdx - jacobian of the dynamics w.r.t state A = D1(f(x,u))
 * @return  d/dt[rho]
 */
inline vec rhodot(const vec& rho, const vec& kldx, const mat& fdx)
{
  return kldx - fdx.t() * rho;
}

/** @brief Store and smaple past states */
class ReplayBuffer
{
public:
  ReplayBuffer(unsigned int buffer_size, unsigned int batch_size)
    : buffer_size_(buffer_size), batch_size_(batch_size)
  {
    // num_ = 0;
  }

  void append(const vec& x)
  {
    if (memory_.size() < buffer_size_)
    {
      memory_.emplace(memory_.size(), x);
      return;
    }
    std::cout << "WARNING: Memory Buffer is full" << std::endl;

    // if (num_ < batch_size_)
    // {
    //   prev_states_.insert_cols(num_, x);
    //   num_++;
    //   return;
    // }
    //
    // prev_states_.cols(0, num_ - 2) = prev_states_.cols(1, num_ - 1);
    // prev_states_.col(num_ - 1) = x;
  }

  void sampleMemory(mat& xt_total, const mat& xt)
  {
    // if (num_ == 0)
    // {
    //   xt_total = xt;
    //   return;
    // }
    //
    // const auto num_states = xt.n_cols + num_;
    // xt_total.resize(xt.n_rows, num_states);
    //
    // xt_total.cols(0, num_ - 1) = prev_states_;
    // xt_total.cols(num_, num_states - 1) = xt;

    if (memory_.empty())
    {
      xt_total = xt;
    }

    // Concatenate the current store states with predicted trajectory
    else if (memory_.size() <= batch_size_)
    {
      const auto num_stored = memory_.size();
      const auto num_states = xt.n_cols + num_stored;
      xt_total.resize(xt.n_rows, num_states);

      for (unsigned int i = 0; i < num_stored; i++)
      {
        // Index is the key
        xt_total.col(i) = memory_.at(i);
      }

      // Copy predicted trajectory to end
      xt_total.cols(num_stored, num_states - 1) = xt;
    }

    // Randomly sample memory and concatenate with predicted trajectory
    else
    {
      const auto num_states = xt.n_cols + batch_size_;
      xt_total.resize(xt.n_rows, num_states);

      // random ints on interval [a b]
      const ivec rand_ints = randi<ivec>(batch_size_, distr_param(0, memory_.size() - 1));

      for (unsigned int i = 0; i < batch_size_; i++)
      {
        // Index is the key
        xt_total.col(i) = memory_.at(rand_ints(i));
      }

      // Copy predicted trajectory to end
      xt_total.cols(batch_size_, num_states - 1) = xt;
    }
  }

private:
  // mat prev_states_;
  // unsigned int num_;
  unsigned int buffer_size_;                      // total number of past states in memory
  unsigned int batch_size_;                       // number of states sampled from memory
  std::unordered_map<unsigned int, vec> memory_;  // past states
};

/** @brief Receeding horizon ergodic trajectory optimization */
template <class ModelT>
class ErgodicControl
{
public:
  ErgodicControl(const ModelT& model, double dt, double horizon, unsigned int num_samples,
                 unsigned int buffer_size, unsigned int batch_size, mat R, mat Sigma);

  /**
   * @brief Update the control signal
   * @param x - current state [x, y, theta] (column vector)
   * @return first control of the updated control signal
   */
  vec control(const GridMap& grid, const vec& x);

private:
  void sample(const GridMap& grid);

  void trajStat(const mat& xt);

  void ergMeasDeriv(const mat& xt);

  void updateControl(const mat& xt, const mat& rhot);

private:
  ModelT model_;
  double dt_;
  double horizon_;
  unsigned int num_samples_;
  unsigned int steps_;
  mat Rinv_;
  mat Sigmainv_;
  mat ut_;
  mat s_;
  mat edx_;
  vec rhoT_;
  vec p_;
  vec q_;
  RungeKutta rk4_;
  ReplayBuffer buffer_;
  CoStateFunc rhodot_;
};

template <class ModelT>
ErgodicControl<ModelT>::ErgodicControl(const ModelT& model, double dt, double horizon,
                                       unsigned int num_samples, unsigned int buffer_size,
                                       unsigned int batch_size, mat R, mat Sigma)
  : model_(model)
  , dt_(dt)
  , horizon_(horizon)
  , num_samples_(num_samples)
  , steps_(static_cast<unsigned int>(horizon / dt))
  , Rinv_(inv(R))
  , Sigmainv_(inv(Sigma))
  , ut_(model.action_space, steps_, arma::fill::zeros)
  , s_(2, num_samples)
  , edx_(model.state_space, steps_)
  , rhoT_(model.state_space, arma::fill::zeros)
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
                      std::placeholders::_3);
}

template <class ModelT>
vec ErgodicControl<ModelT>::control(const GridMap& grid, const vec& x)
{
  // Shift columns to the left by 1 and set last column to zeros
  ut_.cols(0, ut_.n_cols - 2) = ut_.cols(1, ut_.n_cols - 1);
  ut_.col(ut_.n_cols - 1).fill(0.0);

  // Forward simulation
  const mat xt = rk4_.solve(model_, x, ut_, horizon_);
  // xt.print("xt:");

  // Sample (x,y) space
  sample(grid);
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

  // Backwards pass
  const mat rhot = rk4_.solve(rhodot_, model_, rhoT_, xt, ut_, edx_, horizon_);
  // rhot.print("rho:");

  // Update controls
  updateControl(xt, rhot);
  // ut_.print("Control:");

  // add current state to memory
  buffer_.append(x);

  return ut_.col(0);
}

template <class ModelT>
void ErgodicControl<ModelT>::sample(const GridMap& grid)
{
  // Uniformly generate sample [0 1]
  s_.randu(2, num_samples_);
  // s_.print();

  // Scale based on grid domain
  s_.row(0) = s_.row(0) * (grid.xmax() - grid.xmin()) + grid.xmin();
  s_.row(1) = s_.row(1) * (grid.ymax() - grid.ymin()) + grid.ymin();

  // Sample grid
  for (unsigned int i = 0; i < num_samples_; i++)
  {
    // Probability of occupancy [0 1] and -1 for unknown
    // p_(i) = entropy(grid.getCell(s_(0, i), s_(1, i)));
    p_(i) = grid.getCell(s_(0, i), s_(1, i));
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

    // TODO: add 1e-8 is this is 0 ??
    q_(i) = qval;  // + 1e-8;
  }

  // Normalize the trajectory statistics
  const auto sum = accu(q_);
  if (!almost_equal(0.0, sum))
  {
    q_ /= sum;
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

    // TODO: add 1e-8 is this is 0 ??
    // kldx.rows(0, 1) += 1e-8;
    edx_.col(i) = kldx;
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

    // ut_(0, i) = std::clamp(ut_(0, i), -0.1, 0.1);
    // ut_(1, i) = std::clamp(ut_(1, i), -0.1, 0.1);
    // ut_(2, i) = std::clamp(ut_(2, i), -0.5, 0.5);

    // ut_(0,i) = std::clamp(ut_(0,i), -1.0, 1.0);
    // ut_(1,i) = std::clamp(ut_(1,i), -1.0, 1.0);
    // ut_(2,i) = std::clamp(ut_(2,i), -1.0, 1.0);
    // ut_(3,i) = std::clamp(ut_(3,i), -1.0, 1.0);

    if (any(ut_.col(i) > 1.0) || any(ut_.col(i) < -1.0))
    {
      // euclidean norm
      ut_.col(i) /= norm(ut_.col(i), 2);
    }
  }
}

}  // namespace ergodic_exploration
