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
#include <visualization_msgs/Marker.h>
#include <tf2/LinearMath/Quaternion.h>

#include <ergodic_exploration/buffer.hpp>
#include <ergodic_exploration/integrator.hpp>
#include <ergodic_exploration/target.hpp>

namespace ergodic_exploration
{
using arma::span;

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
  // return kldx - dbar - fdx.t() * rho;
  return kldx - fdx.t() * rho;

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
  vec control(const Collision& collision, const GridMap& grid, const Target& target, const vec& x);

  void path(nav_msgs::Path& path, std::string frame) const;

  void samples(visualization_msgs::Marker& marker, std::string frame) const;

private:
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
  void ergMeasDeriv();

  /**
   * @brief Log loss barrier function to obstacles of the form -log(b-x)
   * @param collision - collision detector
   * @param grid - grid map
   * @param xt - predicted trajectory
   * @details Updates a matrix conatining the log loss barrier function derivatve
   * for each state in xt
   */
  void barrier(const Collision& collision, const GridMap& grid);

  /**
   * @brief Update the control signal
   * @param xt - predicted trajectory
   * @param rhot - co-state variable solution
   * @details rhot is assumed to already be sorted from [t0 tf] with t0 at index 0
   */
  void updateControl(const mat& rhot);

private:
  ModelT model_;              // robot's dynamic model
  double dt_;                 // time step in integration
  double horizon_;            // control horizon
  unsigned int num_samples_;  // number of samples drawn from map
  unsigned int steps_;        // number of steps used in integration
  mat Rinv_;                  // inverse of the control weight matrix
  mat Sigmainv_;              // inverse of the robot's (x,y) position uncertainty
  mat ut_;                    // control signal
  mat s_;                     // pair of (x,y) points sampled from map
  mat edx_;                   // ergodic measure derivatives
  mat bdx_;                   // log loss barrier derivatives
  mat xt_;                    // current trajectory
  vec rhoT_;                  // co-state terminal condition
  vec p_;                     // normalized target distribution samples
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
  , Rinv_(inv(R))
  , Sigmainv_(inv(Sigma))
  , ut_(model.action_space, steps_, arma::fill::zeros)
  , s_(2, num_samples)  // exploration space is set to (x,y)
  , edx_(model.state_space, steps_, arma::fill::zeros)
  , bdx_(model.state_space, steps_, arma::fill::zeros)  // set heading column to zeros
  , xt_(model.state_space, steps_)
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

  // arma::arma_rng::set_seed_random();

  rhodot_ = std::bind(rhodot, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, std::placeholders::_4);
}

template <class ModelT>
vec ErgodicControl<ModelT>::control(const Collision& collision, const GridMap& grid,
                                    const Target& target, const vec& x)
{
  // Shift columns to the left by 1 and set last column to zeros
  ut_.cols(0, ut_.n_cols - 2) = ut_.cols(1, ut_.n_cols - 1);
  ut_.col(ut_.n_cols - 1).fill(0.0);

  // Forward simulation
  rk4_.solve(xt_, model_, x, ut_, horizon_);
  // xt.print("xt:");

  // Sample (x,y) space
  target.sample(s_, p_, grid, num_samples_);
  // p_.print("p:");

  // Sample past states
  mat xt_total;
  buffer_.sampleMemory(xt_total, xt_);

  // std::cout << xt_total.n_cols << std::endl;


  // Time average trajectory statistics using memory and predicted trajectory
  trajStat(xt_total);
  // q_.print("q:");

  // Derivative of the ergodic measure w.r.t the state
  ergMeasDeriv();
  // edx_.print("edx:");

  // Derivative of the barrier function
  // barrier(collision, grid);

  // Backwards pass
  mat rhot;
  rk4_.solve(rhot, rhodot_, model_, rhoT_, xt_, ut_, edx_, bdx_, horizon_);
  // rhot.print("rho:");

  max(rhot, 1).print("max rho");
  min(rhot, 1).print("min rho");

  // std::cout << "max rho x: " << max(rhot.row(0)) << std::endl;
  // std::cout << "max rho y: " << max(rhot.row(1)) << std::endl;
  // std::cout << "max rho th: " << max(rhot.row(2)) << std::endl;
  //
  // std::cout << "min rho x: " << min(rhot.row(0)) << std::endl;
  // std::cout << "min rho y: " << min(rhot.row(1)) << std::endl;
  // std::cout << "min rho th: " << min(rhot.row(2)) << std::endl;



  // Update controls
  updateControl(rhot);
  // ut_.print("Control:");

  // add current state to memory
  buffer_.append(x);

  return ut_.col(0);
}

template <class ModelT>
void ErgodicControl<ModelT>::path(nav_msgs::Path& path, std::string frame) const
{
  path.header.frame_id = frame;
  path.poses.resize(xt_.n_cols);
  for (unsigned int i = 0; i < xt_.n_cols; i++)
  {
    path.poses.at(i).pose.position.x = xt_(0, i);
    path.poses.at(i).pose.position.y = xt_(1, i);

    tf2::Quaternion quat(xt_(2, i), 0.0, 0.0);
    path.poses.at(i).pose.orientation.x = quat.x();
    path.poses.at(i).pose.orientation.y = quat.y();
    path.poses.at(i).pose.orientation.z = quat.z();
    path.poses.at(i).pose.orientation.w = quat.w();
  }
}

template <class ModelT>
void ErgodicControl<ModelT>::samples(visualization_msgs::Marker& marker, std::string frame) const
{
  marker.header.frame_id = frame;
  marker.id = 0;
  marker.type = visualization_msgs::Marker::POINTS;
  marker.action = visualization_msgs::Marker::ADD;
  marker.lifetime = ros::Duration(0.0);
  marker.scale.x = 0.01;
  marker.scale.y = 0.01;

  marker.color.r = 0.39;
  marker.color.g = 1.0;
  marker.color.b = 0.0;
  marker.color.a = 1.0;

  marker.points.resize(s_.n_cols);

  for (unsigned int i = 0; i < s_.n_cols; i++)
  {
    marker.points.at(i).x = s_(0, i);
    marker.points.at(i).y = s_(1, i);
  }
}

template <class ModelT>
void ErgodicControl<ModelT>::trajStat(const mat& xt)
{
  // for (unsigned int i = 0; i < num_samples_; i++)
  // {
  //   vec qvals(xt.n_cols);
  //   for (unsigned int j = 0; j < xt.n_cols; j++)
  //   {
  //     const vec diff = s_.col(i) - xt(span(0, 1), span(j, j));
  //     qvals(j) = -0.5 * dot(diff.t() * Sigmainv_, diff);
  //   }
  //
  //   const auto c = max(qvals);
  //   q_(i) = c + log(sum(exp(qvals - c)));
  //
  //   if(almost_equal(0.0, q_(i)))
  //   {
  //     std::cout << "WARNING: trajectory statistics are zero" << std::endl;
  //     q_(i) = 1e-8;
  //   }
  // }
  //
  // q_ /= sum(q_);


  for (unsigned int i = 0; i < num_samples_; i++)
  {
    auto qval = 0.0;
    for (unsigned int j = 0; j < xt.n_cols; j++)
    {
      const vec diff = s_.col(i) - xt(span(0, 1), span(j, j));
      qval += std::exp(-0.5 * dot(diff.t() * Sigmainv_, diff));
    }

    if (almost_equal(0.0, qval))
    {
      // std::cout << "WARNING: trajectory statistics are zero" << std::endl;
      qval = 1e-8;
    }

    q_(i) = qval;
  }

  // Normalize the trajectory statistics
  const auto total = sum(q_);
  if (!almost_equal(0.0, total))
  {
    q_ /= total;
  }


  // for (unsigned int i = 0; i < num_samples_; i++)
  // {
  //   // (si - xt)
  //   mat m = -xt.rows(0, 1);
  //   m.each_col() += s_.col(i);
  //   m = square(m.t()) * Sigmainv_;
  //
  //   // log pdf of each state in xt
  //   const vec v = -0.5 * sum(m, 1);
  //
  //   const auto c = max(v);
  //   q_(i) = c + log(sum(exp(v - c)));
  //
  //   if(almost_equal(0.0, q_(i)))
  //   {
  //     std::cout << "WARNING: trajectory statistics are zero" << std::endl;
  //     q_(i) = 1e-8;
  //   }
  // }
  //
  // q_ /= sum(q_);
}

template <class ModelT>
void ErgodicControl<ModelT>::ergMeasDeriv()
{
  for (unsigned int i = 0; i < steps_; i++)
  {
    vec kldx(3, arma::fill::zeros);
    for (unsigned int j = 0; j < num_samples_; j++)
    {
      const vec diff = s_.col(j) - xt_(span(0, 1), span(i, i));
      const auto qval = std::exp(-0.5 * dot(diff.t() * Sigmainv_, diff));
      kldx.rows(0, 1) += (p_(j) / q_(j)) * qval * (Sigmainv_ * diff);
    }

    edx_.col(i) = kldx;
  }


  // vec qbar_max(num_samples_);
  // for (unsigned int i = 0; i < num_samples_; i++)
  // {
  //   vec qvals(steps_);
  //   for (unsigned int j = 0; j < steps_; j++)
  //   {
  //     const vec diff = s_.col(i) - xt_(span(0, 1), span(j, j));
  //     qvals(j) = -0.5 * dot(diff.t() * Sigmainv_, diff);
  //   }
  //
  //   const auto c = max(qvals);
  //   // qbar(i) = c + log(sum(exp(qvals - c)));
  //   qbar_max(i) = c;
  // }
  //
  //
  // for (unsigned int i = 0; i < steps_; i++)
  // {
  //   vec kldx(3, arma::fill::zeros);
  //   for (unsigned int j = 0; j < num_samples_; j++)
  //   {
  //     const vec diff = s_.col(j) - xt_(span(0, 1), span(i, i));
  //     const auto qbar = std::exp(-0.5 * dot(diff.t() * Sigmainv_, diff) - qbar_max(j));
  //     // const auto qbar = std::exp(-0.5 * dot(diff.t() * Sigmainv_, diff));
  //
  //     const vec dq = qbar * (Sigmainv_ * diff);
  //
  //     kldx.rows(0, 1) += (p_(j) / q_(j)) * dq;
  //   }
  //
  //   edx_.col(i) = kldx;
  // }
}

template <class ModelT>
void ErgodicControl<ModelT>::barrier(const Collision& collision, const GridMap& grid)
{
  const auto weight = 1e12;
  const auto eps = collision.totalPadding();

  for (unsigned int i = 0; i < steps_; i++)
  {
    vec dbar;
    if (collision.minDirection(dbar, grid, xt_.col(i)) == CollisionMsg::obstacle)
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
void ErgodicControl<ModelT>::updateControl(const mat& rhot)
{
  for (unsigned int i = 0; i < xt_.n_cols; i++)
  {
    // rhot is already sorted from t0 to tf
    ut_.col(i) = -Rinv_ * model_.fdu(xt_.col(i)).t() * rhot.col(i);

    // TODO: add control limits as a parameter
    // TODO: apply filter to control signal (savgol filter)

    // ut_.col(i).print("u:");

    // see armadillo clamp()
    // ut_(0, i) = std::clamp(ut_(0, i), -0.1, 0.1);
    // ut_(1, i) = std::clamp(ut_(1, i), -0.1, 0.1);
    // ut_(2, i) = std::clamp(ut_(2, i), -0.5, 0.5);

    // ut_(0, i) = std::clamp(ut_(0, i), -1.0, 1.0);
    // ut_(1, i) = std::clamp(ut_(1, i), -1.0, 1.0);
    // ut_(2, i) = std::clamp(ut_(2, i), -2.0, 2.0);

    if (any(ut_.col(i) > 1.0) || any(ut_.col(i) < -1.0))
    {
      // euclidean norm
      ut_.col(i) /= norm(ut_.col(i), 2);
    }
  }
}
}  // namespace ergodic_exploration
#endif
