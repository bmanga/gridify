#include "core/pdb.h"
#include "core/common.h"

#include <vector>

struct site_ligand_stats {
  double site_volume;
  double ligand_volume;
  double intersection_volume;
  double union_volume;
};

Surface_Mesh gen_ligand_geometry(const std::vector<pdb_atom_entry> &ligand,
                                 double scale_radius,
                                 bool pca_align);

site_ligand_stats calc_grid_ligand_stats(const std::vector<Point_3> &grid,
                                         const Surface_Mesh &ligand_geom,
                                         double grid_radius,
                                         bool is_grid_packed);