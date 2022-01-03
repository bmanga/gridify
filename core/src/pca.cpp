#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include "core/common.h"

using Mat = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

static const auto view_as_eigen_mat(const std::vector<Point_3> &points)
{
  static_assert(sizeof(Point_3) == sizeof(double) * 3);
  const auto m =
      Eigen::Map<const Mat>((const double *)points.data(), points.size(), 3);
  return m;
}

static auto view_as_eigen_mat(std::vector<Point_3> &points)
{
  static_assert(sizeof(Point_3) == sizeof(double) * 3);
  auto m = Eigen::Map<Mat>((double *)points.data(), points.size(), 3);
  return m;
}

std::vector<Point_3> pca_aligned_points(const std::vector<Point_3> &points)
{
  if (points.empty()) {
    return {};
  }
  const auto m = view_as_eigen_mat(points);
  auto result = std::vector<Point_3>(points.size());
  auto res_eigen = view_as_eigen_mat(result);
  Eigen::VectorXd means = m.colwise().mean();
  auto centered = m.rowwise() - means.transpose();

  // Covariance matrix
  auto cov = centered.adjoint() * centered / double(centered.rows() - 1);

  auto eig = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(cov);

  res_eigen = centered * Eigen::Reverse<Eigen::MatrixXd, Eigen::Horizontal>(
                      eig.eigenvectors());
  return result;
}

bounds get_bounds(const std::vector<Point_3>& points) {
  auto mat = view_as_eigen_mat(points);
  auto cols = mat.colwise();
  auto mins = cols.minCoeff();
  auto maxs = cols.maxCoeff();

  return
  {
    .x = {mins[0], maxs[0]}, .y = {mins[1], maxs[1]}, .z = {mins[2], maxs[2]}
  };
}
