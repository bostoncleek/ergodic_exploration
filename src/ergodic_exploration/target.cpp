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

Target::Target(const GaussianList& gaussians)
  : gaussians_(gaussians)
{
}

double Target::evaluate(const vec& pt) const
{
  auto val = 0.0;
  for (const auto& gaussian : gaussians_)
  {
    val += gaussian(pt);
  }
  return val;
}

void Target::sample(mat& points, vec& phi, const GridMap& grid, unsigned int num_samples) const
{
  // Uniformly generate sample [0 1]
  points.randu(2, num_samples);

  // Scale based on grid domain
  points.row(0) = points.row(0) * (grid.xmax() - grid.xmin()) + grid.xmin();
  points.row(1) = points.row(1) * (grid.ymax() - grid.ymin()) + grid.ymin();


  for (unsigned int i = 0; i < num_samples; i++)
  {
    // phi(i) = entropy(grid.getCell(points(0, i), points(1, i)));
    phi(i) = evaluate(points.col(i));
    // phi(i) = evaluate(points.col(i)) + entropy(grid.getCell(points(0, i), points(1, i)));
  }

  // const ivec rand_ints = randi<ivec>(num_samples, distr_param(0, gaussians_.size() - 1));
  //
  // vec pt(2);
  // unsigned int i = 0;
  // while (i < num_samples)
  // {
  //   if (!mvnrnd(pt, gaussians_.at(rand_ints(i)).mu, gaussians_.at(rand_ints(i)).cov))
  //   {
  //     std::cout << "FAILURE! Unable to sample target distribution. " << std::endl;
  //   }
  //
  //   if ((pt(0) > grid.xmax()) || (pt(0) < grid.xmin()) || (pt(1) > grid.ymax()) ||
  //       (pt(1) < grid.ymin()))
  //   {
  //     continue;
  //   }
  //
  //   // pt.print("sample:");
  //
  //   points.col(i) = pt;
  //
  //   // evaluate all gaussians
  //   phi(i) = evaluate(pt);
  //
  //   // evaluate the gaussian it was drawn from
  //   // phi(i) = gaussians_.at(rand_ints(i))(pt);
  //   // phi(i) = gaussians_.at(0)(pt) + entropy(grid.getCell(points(0, i), points(1, i)));
  //   i++;
  // }

  // Normalize the sampled values
  const auto total = sum(phi);
  if (!almost_equal(0.0, total))
  {
    phi /= total;
  }
}

void Target::markers(visualization_msgs::MarkerArray& marker_array, std::string frame)
{
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
    marker_array.markers.at(i).lifetime = ros::Duration(0.0);;
    marker_array.markers.at(i).color.r = 1.0;
    marker_array.markers.at(i).color.g = 1.0;
    marker_array.markers.at(i).color.b = 0.6;
    marker_array.markers.at(i).color.a = 0.5;
  }
}
}  // namespace ergodic_exploration
