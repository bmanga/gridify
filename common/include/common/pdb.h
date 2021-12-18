#ifndef GRIDIFY_PDB_H
#define GRIDIFY_PDB_H

#include "common/common.h"
#include <vector>
#include <string>


struct pdb_atom_entry {
  Point_3 pos;
  std::string residue;
  std::string atom;
  int atom_id;
  int residue_id;
  std::string chain;
  std::string kind;
};

struct pdb_frame {
  int frame_idx;
  std::vector<pdb_atom_entry> atoms;
};

struct pdb {
  std::vector<pdb_frame> frames;
};



void parse_pdb(std::ifstream &ifs, struct producer_consumer_queue &queue);

std::vector<pdb_frame> parse_pdb(const std::string &file);

double parse_pdb_gridify_spacing(std::ifstream &ifs);

#endif  // GRIDIFY_PDB_H
