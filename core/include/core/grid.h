#include <vector>
#include "core/common.h"
#include "core/aabb_tree.hpp"
#include "core/points_checker.h"

struct binding_site {
  std::vector<pdb_atom_entry> atoms;
  bounds bounds;
};

binding_site get_binding_site(const std::vector<pdb_atom_entry> &protein,
                              const std::vector<int> &residues);

std::vector<Point_3> surface_only(const std::vector<Point_3> &points,
                                  double radius);

std::vector<Point_3> gen_site_grid(
    const std::vector<pdb_atom_entry> &protein,
    const binding_site &site,
    bool rm_atom_overlaps,
    bool largest_cluster_only,
    bool dense_packing,
    double rm_lc_cutoff,
    double scale_radius,
    double spacing,
    double point_radius,
    double rm_lc_tangent_weight,
    double rm_lc_proximity_weight);