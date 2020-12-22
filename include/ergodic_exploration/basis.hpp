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
 * @file basis.hpp
 * @author Boston Cleek
 * @date 20 Nov 2020
 * @brief Fourier cosine basis
 */
#ifndef BASIS_HPP
#define BASIS_HPP

#include <armadillo>

namespace ergodic_exploration
{
using arma::imat;
using arma::mat;
using arma::vec;

/**
 * @brief Fourier cosine basis */
class Basis
{
private:
  template <class ModelT>
  friend class ErgodicControl;

public:
  /**
   * @brief Constructor
   * @param lx - x-axis lenght of domain
   * @param ly - y-axis lenght of domain
   * @param num_basis - number of each basis functions per dimension
   */
  Basis(double lx, double ly, unsigned int num_basis);

  /**
   * @brief Compose cosine basis functions given positon
   * @param x - robot position [x y]
   * @return cosine basis functions (num_basis^2 x 1)
   * @details x must be on domain [0 lx] x [0 ly]
   */
  vec fourierBasis(const vec& x) const;

  /**
   * @brief Compose gradient cosine basis functions
   * @param x - robot position [x y]
   * @return gradient of each basis function (2 x num_basis^2)
   * @details x must be on domain [0 lx] x [0 ly]
   */
  mat gradFourierBasis(const vec& x) const;

  /**
   * @brief Compose trajectory fourier coefficients
   * @param xt - trajectory
   * @return trajectory fourier coefficients (num_basis^2 x 1)
   * @details xt must be on domain [0 lx] x [0 ly] and robot heading is not required
   */
  vec trajCoeff(const mat& xt) const;

  /**
   * @brief Compose spatial fourier coefficients
   * @param phi_vals - target evaluated at each grid cell in phi_grid
   * @param phi_grid - discretization of fourier domain
   * @return spatial fourier coefficients (num_basis^2 x 1)
   * @details phi_grid is not the occupancy grid, the first row contains the
   * x-cordinates and the second row contains the corresponding y-coordinates.
   * The index of the elements in phi_vals must correspond to the columns in phi_grid.
   */
  vec spatialCoeff(const vec& phi_vals, const mat& phi_grid) const;

private:
  double lx_, ly_;            // length of domain
  unsigned int total_basis_;  // number of basis functions
  vec lamdak_;                // frequency coefficients weights
  imat k_;                    // basis number
};
}  // namespace ergodic_exploration
#endif
