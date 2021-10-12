#include <fstream>
#include <thread>
#include <sstream>
#include <fmt/printf.h>
#include "common/pdb.h"
#include "common/processing.h"

std::string trim(const std::string &s)
{
  auto start = s.begin();
  while (start != s.end() && std::isspace(*start)) {
    start++;
  }

  auto end = s.end();
  do {
    end--;
  } while (std::distance(start, end) > 0 && std::isspace(*end));

  return std::string(start, end + 1);
}


void parse_pdb(std::ifstream &ifs, producer_consumer_queue &queue)
{
  std::string str;
  int frame_idx = 0;
  int num_atoms = 0;
  auto frame = pdb_frame{frame_idx++};
  while (std::getline(ifs, str)) {
    str = trim(str);
    if (str.starts_with("CRYST")) {
      continue;
    }
    if (str.starts_with("REMARK")) {
      continue;
    }
    if (str == "END") {
      // Slow down if the consumers can't keep up to avoid memory bloat.
      while (queue.frames.size_approx() > 200) {
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
      }
      queue.frames.enqueue(std::move(frame));

      frame = pdb_frame{frame_idx++};
      continue;
    }
    std::istringstream iss(str);
    std::string unused, atom, residue;
    int residue_id;
    int atom_id;
    // NOTE: We could read lines at specific offsets, but I think this is more
    // robust.
    iss >> unused >> atom_id >> atom >> residue;
    iss >> residue_id;
    if (iss.fail()) {
      iss.clear();
      iss >> unused;
      iss >> residue_id;
    }

    float x, y, z;
    iss >> x >> y >> z;
    auto entry = pdb_atom_entry{{x, y, z}, residue, atom, atom_id, residue_id};
    ++num_atoms;
    frame.atoms.push_back(std::move(entry));
  }

  queue.producer_done = true;

  if (g_verbose) {
    fmt::print("Parsed pdb file: {} frames with {} atoms each\n",
                             frame_idx, num_atoms / frame_idx);
  }
}

double parse_pdb_gridify_spacing(std::ifstream &ifs) {
  std::string str;
  while (std::getline(ifs, str)) {
    if (!str.starts_with("REMARK 100")) {
      break;
    }
    if (str.starts_with("REMARK 100 spacing")) {
      return std::stod(str.substr(21));
    }
  }
  return -1;
}