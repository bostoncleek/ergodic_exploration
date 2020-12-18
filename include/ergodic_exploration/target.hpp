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

/** @brief 2D gaussian */
struct Gaussian
{
  /** @brief Constructor */
  Gaussian()
  {
  }

  /**
   * @brief Constructor
   * @param mu - mean [mean x, mean y]
   * @param sigmas - standard deviations [sigma x, sigma y]
   */
  Gaussian(const vec& mu, const vec& sigmas)
    : mu(mu), cov(arma::diagmat(square(sigmas))), cov_inv(inv(cov))
  {
  }

  /**
   * @brief Evaluate gaussian
   * @param pt - point [x y]
   * @return evaluated gaussian at pt
   */
  double operator()(const vec& pt) const
  {
    const vec diff = pt - mu;
    return std::exp(-0.5 * dot(diff.t() * cov_inv, diff));
  }

  /**
   * @brief Evaluate gaussian
   * @param pt - point [x y]
   * @param trans - translation from map frame to fourier domain
   * @return evaluated gaussian at pt translated by trans
   * @details the translation is used to translate the mean into the fourier domain
   */
  double operator()(const vec& pt, const vec& trans) const
  {
    // DEBUG
    // if (any(mu - trans) < 0.0)
    // {
    //   std::cout << "WARNING: Targert mean not within fourier domain" << std::endl;
    // }

    // translate mu into frame of fourier domain
    const vec diff = pt - (mu - trans);
    return std::exp(-0.5 * dot(diff.t() * cov_inv, diff));
  }

  vec mu;       // mean
  mat cov;      // covariance
  mat cov_inv;  // inverse of covariance
};

/** @brief Target distribution */
class Target
{
public:
  /** @brief Constructor */
  Target();

  /**
   * @brief Constructor
   * @param gaussians - list of target gaussians
   */
  Target(const GaussianList& gaussians);

  /**
   * @brief Adds gaussian to list
   * @param g - gaussians
   */
  void addGaussian(const Gaussian& g);

  /**
   * @brief Remove gaussian from list
   * @param idx - index of gaussian to remove
   */
  void deleteGaussian(unsigned int idx);

  /**
   * @brief Evaluate the list of gaussians
   * @param pt - point [x y]
   * @param trans - translation from map frame to fourier domain
   * @return value of the list of gaussians evaluated at pt translated by trans
   * @details the translation is used to translate the mean into the fourier domain
   */
  double evaluate(const vec& pt, const vec& trans) const;

  /**
   * @brief Evaluate the target distribution
   * @param trans - translation from map frame to fourier domain
   * @param phi_grid - discretization of fourier domain
   * @return target evaluated at each grid cell in phi_grid
   * @details the translation is used to translate the mean into the fourier domain
   */
  vec fill(const vec& trans, const mat& phi_grid) const;

  /**
   * @brief Visualize target distribution
   * @param frame - target frame
   * @return target is visualized as an ellipse
   */
  visualization_msgs::MarkerArray markers(const std::string& frame) const;

private:
  GaussianList gaussians_;  // list of target gaussians
};
}  // namespace ergodic_exploration
#endif
