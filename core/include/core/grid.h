#include <vector>
#include "core/common.h"
#include "core/aabb_tree.hpp"
#include "core/points_checker.h"

std::vector<Point_3> gen_grid(const points_checker &checker,
                              const Surface_Mesh &poly,
                              const bounds &bounds,
                              float spacing,
                              float pt_radius,
                              bool dense_packing);

void rm_low_connectivity_points(abt::tree3d &points,
                                double dist_cutoff,
                                double min_conn_score,
																double tangent_weight,
																double proximity_weight);

void keep_largest_cluster_only(abt::tree3d &grid);