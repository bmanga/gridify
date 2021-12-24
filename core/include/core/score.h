#include "core/pdb.h"
#include "core/common.h"

Surface_Mesh gen_ligand_geometry(const pdb_frame &f,
                                 double scale_radius,
                                 bool pca_align);

struct site_ligand_stats {
  double site_volume;
  double ligand_volume;
  double intersection_volume;
  double union_volume;
};

site_ligand_stats calc_grid_ligand_stats(
    const std::vector<Point_3> &grid_points,
    const Surface_Mesh &ligand,
    double site_r);