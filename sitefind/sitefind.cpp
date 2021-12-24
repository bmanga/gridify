#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include "core/pdb.h"
#include "core/processing.h"
#include "core/radius.h"
#include "core/aabb_tree.hpp"
#include "core/sitefind.h"

using frames = std::vector<pdb_frame>;

#define ALL_SETTINGS            \
  X(std::string, pdb_file, "", "The file to analyze") \
  X(float, distance, 4, "The minimum distance between the ligand and any residue atom to be considered part of the binding site")\
  X(std::string, chain, "", "The chain to analyze")\
  X(std::string, ligand, "", "The name of the ligand") \
  X(bool, ignore_radii, false, "Calculate the distances from the centers only, ignoring the atomic radii") \
  X(bool, all_ligands, false, "Find all ligands and calculate nearby residues to each")
#include "core/cmdline.inc"

bool g_verbose;

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
        auto resids = get_protein_residues_near_ligand(radmatch, frame, ligand, config.distance, config.ignore_radii);
        std::cout << "[" << ligand.resid << " chain \"" << ligand.chain << "\"]: ";
        for (const auto &resid : resids) {
          std::cout << resid << " ";
        }
        std::cout << std::endl;
      }
    }
    else {
      auto resids = get_protein_residues_near_ligand(radmatch, frame, {config.chain, config.ligand}, config.distance, config.ignore_radii);
      for (const auto &resid : resids) {
          std::cout << resid << " ";
        }
        std::cout << std::endl;
    }
  }
  return 0;
}