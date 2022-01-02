#ifndef GRIDIFY_PCA_H
#define GRIDIFY_PCA_H

#include <vector>
#include "common.h"

bounds get_bounds(const std::vector<Point_3> &points);
std::vector<Point_3> pca_aligned_points(const std::vector<Point_3> &points);

#endif  // GRIDIFY_PCA_H
