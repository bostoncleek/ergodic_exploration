/**
 * @file target.hpp
 * @author Boston Cleek
 * @date 17 Nov 2020
 * @brief Target distribution
 */
#ifndef TARGET_HPP
#define TARGET_HPP

#include <armadillo>

#include <visualization_msgs/MarkerArray.h>

#include <ergodic_exploration/grid.hpp>
#include <ergodic_exploration/numerics.hpp>

namespace ergodic_exploration
{
using arma::mat;
using arma::vec;
struct Gaussian;
typedef std::vector<Gaussian> GaussianList;

struct Gaussian
{
  Gaussian()
  {
  }

  Gaussian(const vec& mu, const vec& sigmas)
    : mu(mu), cov(arma::diagmat(square(sigmas))), cov_inv(inv(cov))
  {
  }

  double operator()(const vec& pt) const
  {
    // TODO: include 1/(2pi * sqrt(det(cov))) ????
    const vec diff = pt - mu;
    return std::exp(-0.5 * dot(diff.t() * cov_inv, diff));
  }

  vec mu;
  mat cov;
  mat cov_inv;
};

struct Target
{
public:
  Target(const GaussianList& gaussians);

  double evaluate(const vec& pt) const;

  void fill(vec& phi_vals, const mat& phi_grid) const;

  void markers(visualization_msgs::MarkerArray& marker_array, std::string frame) const;

private:
  GaussianList gaussians_;
};

}  // namespace ergodic_exploration
#endif
