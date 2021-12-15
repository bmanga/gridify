#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include "common/pdb.h"
#include "common/processing.h"
#include "common/radius.h"
#include "common/aabb_tree.hpp"

using frames = std::vector<pdb_frame>;

#define ALL_SETTINGS            \
  X(std::string, pdb_file, "", "The file to analyze") \
  X(float, distance, 4, "The minimum distance between the ligand and any residue atom to be considered part of the binding site")\
  X(std::string, chain, "", "The chain to analyze")\
  X(std::string, ligand, "", "The name of the ligand") \
  X(bool, ignore_radii, false, "Calculate the distances from the centers only, ignoring the atomic radii") \
  X(bool, all_ligands, false, "Find all ligands and calculate nearby residues to each")
#include "common/cmdline.inc"

bool g_verbose;

bool has_intersections(const abt::tree3d &atoms, const std::vector<abt::aabb3d> &ligand_points) {
  for (const auto & lp : ligand_points) {
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

struct ligand {
  std::string chain;
  std::string resid;

  bool operator==(const ligand &) const = default;
};

namespace std {
  template <>
  struct hash<ligand> {
    size_t operator()(const ligand &ligand) const {
      return std::hash<std::string>{}(ligand.chain) + std::hash<std::string>{}(ligand.resid);
    }
  };
}

std::unordered_set<ligand> discover_ligands(const pdb_frame &frame) {
  std::unordered_set<ligand> ligands;

  for (const auto &atom : frame.atoms) {
    if (atom.kind == "HETATM" && atom.residue != "HOH") {
      ligands.insert({atom.chain, atom.residue});
    }
  }
  return ligands;
}

std::vector<int> get_protein_residues_near_ligand(const config &config,
                                                  const radius_matcher &radmatch,
                                                  const pdb_frame &frame,
                                                  const ligand &ligand)
{
  std::unordered_map<int, abt::tree3d> residues;
  std::vector<abt::aabb3d> ligand_points;
  for (const auto &atom : frame.atoms) {
    if (atom.chain == ligand.chain) {
      auto pos = atom.pos;
      double radius = config.ignore_radii ? config.distance / 2 : radmatch.radius(atom);
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

int main(int argc, char *argv[])
{
  auto config = parse_cmd_line_settings("sitefind",
                                        "Find the residues that are close to the ligand",
                                        argc, argv, [](cxxopts::Options &){});

  auto pdb_file = std::ifstream(config.pdb_file);
  auto radii_file = std::ifstream("radii.json");

  auto radmatch = radius_matcher(radii_file, 1, config.distance / 2);

  producer_consumer_queue queue;
  parse_pdb(pdb_file, queue);


  pdb_frame frame;
  while (queue.frames.try_dequeue(frame)) {
    std::cout << "New frame" << std::endl;
    if (config.all_ligands) {
      auto all_ligands = discover_ligands(frame);
      for (const auto &ligand : all_ligands) {
        auto resids = get_protein_residues_near_ligand(config, radmatch, frame, ligand);
        std::cout << "[" << ligand.resid << " chain \"" << ligand.chain << "\"]: ";
        for (const auto &resid : resids) {
          std::cout << resid << " ";
        }
        std::cout << std::endl;
      }
    }
    else {
      auto resids = get_protein_residues_near_ligand(config, radmatch, frame, {config.chain, config.ligand});
      for (const auto &resid : resids) {
          std::cout << resid << " ";
        }
        std::cout << std::endl;
    }
  }
  return 0;
}