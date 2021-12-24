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

std::vector<Point_3> generate_grid_points(const pdb_frame &frame,
                                          const std::vector<int> &site_residues,
                                          bool rm_atom_overlaps,
                                          bool largest_cluster_only,
                                          bool dense_packing,
                                          double rm_lc_cutoff,
                                          double scale_radius,
                                          double spacing,
                                          double point_radius,
                                          double rm_lc_tangent_weight,
                                          double rm_lc_proximity_weight);