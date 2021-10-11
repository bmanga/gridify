#ifndef GRIDIFY_PDB_H
#define GRIDIFY_PDB_H

#include "common.h"
#include <moodycamel/concurrentqueue.h>



struct pdb_atom_entry {
  Point_3 pos;
  std::string residue;
  std::string atom;
  int atom_id;
  int residue_id;
};

struct pdb_frame {
  int frame_idx;
  std::vector<pdb_atom_entry> atoms;
};

struct pdb {
  std::vector<pdb_frame> frames;
};

using frame_queue = moodycamel::ConcurrentQueue<pdb_frame>;

struct producer_consumer_queue {
  std::atomic_bool producer_done = false;
  std::atomic<int> consumers_done = 0;

  frame_queue frames;
};

void parse_pdb(std::ifstream &ifs, producer_consumer_queue &queue);

#endif  // GRIDIFY_PDB_H
