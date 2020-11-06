/**
 * @file ergodic_control.hpp
 * @aut_hor Boston Cleek
 * @date 5 Nov 2020
 * @brief Ergodic control strategy for exploration
 */

#pragma once

#include <cmath>
#include <stdexcept>

#include <ergodic_exploration/grid.hpp>
#include <ergodic_exploration/numerics.hpp>
#include <ergodic_exploration/integrator.hpp>

namespace ergodic_exploration
{
using arma::span;

// TODO: check base 2 or e
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

/** @brief Receeding horizon ergodic trajectory optimization */
template <class ModelT>
class ErgodicControl
{
public:
  ErgodicControl(const ModelT& model, double dt, double horizon, unsigned int num_samples,
                 unsigned int buffer_size, mat R, mat Sigma);

  /**
   * @brief Update the control signal
   * @param x - current state [x, y, theta] (column vector)
   * @return first control of the updated control signal
   */
  vec control(const GridMap& grid, const vec& x);

  void sample(const GridMap& grid);

  void trajStat(const mat& xt);

  void ergMeasDeriv(const mat& xt);

private:
  ModelT model_;
  double dt_;
  double horizon_;
  unsigned int num_samples_;
  unsigned int buffer_size_;
  unsigned int steps_;
  mat R_;
  mat Sigmainv_;
  mat ut_;
  mat s_;
  mat edx_;
  vec pT_;
  vec p_;
  vec q_;
  RungeKutta rk4_;
};

template <class ModelT>
ErgodicControl<ModelT>::ErgodicControl(const ModelT& model, double dt, double horizon,
                                       unsigned int num_samples, unsigned int buffer_size,
                                       mat R, mat Sigma)
  : model_(model)
  , dt_(dt)
  , horizon_(horizon)
  , num_samples_(num_samples)
  , buffer_size_(buffer_size)
  , steps_(static_cast<unsigned int>(horizon / dt))
  , R_(R)
  , Sigmainv_(inv(Sigma))
  , ut_(model.action_space_, steps_, arma::fill::zeros)
  , s_(2, num_samples)
  , edx_(model.state_space_, steps_)
  , pT_(model.state_space_, arma::fill::zeros)
  , p_(num_samples)
  , q_(num_samples)
  , rk4_(dt)
{
  if (steps_ == 1)
  {
    throw std::invalid_argument("Need at least two steps in forward simulation. \
                                 Increase the horizon or decrease the time step.");
  }
}

template <class ModelT>
vec ErgodicControl<ModelT>::control(const GridMap& grid, const vec& x)
{
  // Shift columns to the left by 1 and set last column to zeros
  ut_.cols(0, ut_.n_cols - 2) = ut_.cols(1, ut_.n_cols - 1);
  ut_.col(ut_.n_cols - 1).fill(0.0);

  // Forward simulation
  const mat xt = rk4_.solve(model_, x, ut_, horizon_);

  // TODO: add buffer

  // Sample (x,y) space
  sample(grid);

  // TODO: append xt to buffer and pass in here
  // Time averages trajectory statistics
  trajStat(xt);

  // Derivative of the ergodic measure w.r.t the state
  ergMeasDeriv(xt);

  // Backwards pass

  // Update controls

  return ut_.col(0);
}

template <class ModelT>
void ErgodicControl<ModelT>::sample(const GridMap& grid)
{
  // Uniformly generate sample [0 1]
  s_.randu(2, num_samples_);

  // Scale based on grid domain
  s_.row(0) = s_.row(0) * (grid.xmax() - grid.xmin()) + grid.xmin();
  s_.row(1) = s_.row(1) * (grid.ymax() - grid.ymin()) + grid.ymin();

  // Sample grid
  for (unsigned int i = 0; i < num_samples_; i++)
  {
    // Probability of occupancy [0 1] and -1 for unknown
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
    }

    q_(i) = qval;
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

    edx_.col(i) = kldx;
  }
}

}  // namespace ergodic_exploration
