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
