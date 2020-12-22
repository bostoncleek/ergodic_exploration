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
 * @file target.cpp
 * @author Boston Cleek
 * @date 17 Nov 2020
 * @brief Target distribution
 */

#include <cmath>

#include <ergodic_exploration/target.hpp>
#include <ergodic_exploration/numerics.hpp>

namespace ergodic_exploration
{
using arma::distr_param;
using arma::ivec;
using arma::randi;

Target::Target()
{
}

Target::Target(const GaussianList& gaussians) : gaussians_(gaussians)
{
}

void Target::addGaussian(const Gaussian& g)
{
  gaussians_.emplace_back(g);
}

void Target::deleteGaussian(unsigned int idx)
{
  gaussians_.erase(gaussians_.begin() + idx);
}

double Target::evaluate(const vec& pt, const vec& trans) const
{
  auto val = 0.0;
  for (const auto& gaussian : gaussians_)
  {
    val += gaussian(pt, trans);
  }
  return val;
}

vec Target::fill(const vec& trans, const mat& phi_grid) const
{
  vec phi_vals(phi_grid.n_cols);
  for (unsigned int i = 0; i < phi_grid.n_cols; i++)
  {
    phi_vals(i) = evaluate(phi_grid.col(i), trans);
  }

  // normalize distribution values to sum to 1
  phi_vals /= sum(phi_vals);
  return phi_vals;
}

visualization_msgs::MarkerArray Target::markers(const std::string& frame) const
{
  visualization_msgs::MarkerArray marker_array;
  marker_array.markers.resize(gaussians_.size());

  for (unsigned int i = 0; i < gaussians_.size(); i++)
  {
    const vec eigval = eig_sym(gaussians_.at(i).cov);

    marker_array.markers.at(i).header.frame_id = frame;
    marker_array.markers.at(i).id = i;
    marker_array.markers.at(i).type = visualization_msgs::Marker::SPHERE;
    marker_array.markers.at(i).action = visualization_msgs::Marker::ADD;
    marker_array.markers.at(i).pose.position.x = gaussians_.at(i).mu(0);
    marker_array.markers.at(i).pose.position.y = gaussians_.at(i).mu(1);
    marker_array.markers.at(i).pose.orientation.w = 1.0;
    marker_array.markers.at(i).scale.x = 2.0 * eigval(0);
    marker_array.markers.at(i).scale.y = 2.0 * eigval(1);
    marker_array.markers.at(i).scale.z = 0.01;
    marker_array.markers.at(i).lifetime = ros::Duration(0.0);
    ;
    marker_array.markers.at(i).color.r = 1.0;
    marker_array.markers.at(i).color.g = 1.0;
    marker_array.markers.at(i).color.b = 0.6;
    marker_array.markers.at(i).color.a = 0.5;
  }

  return marker_array;
}
}  // namespace ergodic_exploration
