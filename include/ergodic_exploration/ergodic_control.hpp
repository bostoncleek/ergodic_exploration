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
  return kldx - fdx * rho;
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
  mat Rinv_;
  mat Sigmainv_;
  mat ut_;
  mat s_;
  mat edx_;
  vec rhoT_;
  vec p_;
  vec q_;
  RungeKutta rk4_;
  CoStateFunc rhodot_;
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
  , Rinv_(inv(R))
  , Sigmainv_(inv(Sigma))
  , ut_(model.action_space_, steps_, arma::fill::zeros)
  , s_(2, num_samples)
  , edx_(model.state_space_, steps_)
  , rhoT_(model.state_space_, arma::fill::zeros)
  , p_(num_samples)
  , q_(num_samples)
  , rk4_(dt)
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
  // xt.print();

  // TODO: add buffer

  // Sample (x,y) space
  sample(grid);
  // p_.print("p:");

  // TODO: append xt to buffer and pass in here
  // Time averages trajectory statistics
  trajStat(xt);
  // q_.print("q:");

  // Derivative of the ergodic measure w.r.t the state
  ergMeasDeriv(xt);
  // edx_.print("edx:");

  // Backwards pass
  const mat rhot = rk4_.solve(rhodot_, model_, rhoT_, xt, ut_, edx_, horizon_);
  // rhot.print("rho:");

  // Update controls
  for (unsigned int i = 0; i < xt.n_cols; i++)
  {
    // index into rho backwards
    const unsigned int j = (xt.n_cols - 1) - i;
    ut_.col(i) = -Rinv_ * model_.fdu(xt.col(i)).t() * rhot.col(j);

    // TODO: add control limits as a parameter
    if (any(ut_.col(i) > 1.0) || any(ut_.col(i) < -1.0))
    {
      ut_.col(i) /= accu(ut_.col(i));
    }
  }

  // ut_.print("Control:");

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

}  // namespace ergodic_exploration
