#include <string>
#include <unordered_map>
#include "core/pdb.h"
#include "core/radius.h"

struct ligand {
  std::string chain;
  std::string resid;

  bool operator==(const ligand &) const = default;
};

namespace std {
template <>
struct hash<ligand> {
  size_t operator()(const ligand &ligand) const
  {
    return std::hash<std::string>{}(ligand.chain) +
           std::hash<std::string>{}(ligand.resid);
  }
};
}  // namespace std

std::unordered_map<ligand, std::vector<pdb_atom_entry>> discover_ligands(
    const pdb_frame &frame);

std::vector<int> get_protein_residues_near_ligand(
    const pdb_frame &frame,
    const ligand &ligand,
    double distance,
    bool ignore_radii);

std::vector<int> get_protein_residues_near_ligand(
    const std::vector<pdb_atom_entry> &protein,
    const std::vector<pdb_atom_entry> &ligand,
    double distance,
    bool ignore_radii);