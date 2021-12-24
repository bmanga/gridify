#include <unordered_map>
#include "core/aabb_tree.hpp"
#include "core/sitefind.h"

std::unordered_set<ligand> discover_ligands(const pdb_frame &frame)
{
  std::unordered_set<ligand> ligands;

  for (const auto &atom : frame.atoms) {
    if (atom.kind == "HETATM" && atom.residue != "HOH") {
      ligands.insert({atom.chain, atom.residue});
    }
  }
  return ligands;
}

static bool has_intersections(const abt::tree3d &atoms,
                       const std::vector<abt::aabb3d> &ligand_points)
{
  for (const auto &lp : ligand_points) {
    if (atoms.any_overlap(lp, [&](unsigned int id) {
          auto bb = atoms.get_aabb(id);
          auto r1 = (bb.upperBound[0] - bb.lowerBound[0]) / 2;
          auto r2 = (lp.upperBound[0] - lp.lowerBound[0]) / 2;
          auto dist = r1 + r2;
          auto dx = bb.centre[0] - lp.centre[0];
          auto dy = bb.centre[1] - lp.centre[1];
          auto dz = bb.centre[2] - lp.centre[2];
          return dx * dx + dy * dy + dz * dz < dist * dist;
        })) {
      return true;
    }
  }
  return false;
}

std::vector<int> get_protein_residues_near_ligand(
    const pdb_frame &frame,
    const ligand &ligand,
    double distance,
    bool ignore_radii
)
{
  std::unordered_map<int, abt::tree3d> residues;
  std::vector<abt::aabb3d> ligand_points;
  const auto &radmatch = radius_matcher::get();
  for (const auto &atom : frame.atoms) {
    if (atom.chain == ligand.chain) {
      auto pos = atom.pos;
      double radius = distance / 2;
      if (!ignore_radii) {
        radius += radmatch.radius(atom);
      }
      auto bb = abt::aabb3d::of_sphere({pos.x(), pos.y(), pos.z()}, radius);
      if (atom.kind == "ATOM") {
        residues[atom.residue_id].insert(bb);
      }
      else if (atom.kind == "HETATM" && atom.residue == ligand.resid) {
        ligand_points.push_back(bb);
      }
    }
  }

  std::vector<int> intersecting_residues;
  for (const auto &[id, atoms] : residues) {
    if (has_intersections(atoms, ligand_points)) {
      intersecting_residues.push_back(id);
    }
  }
  return intersecting_residues;
}
