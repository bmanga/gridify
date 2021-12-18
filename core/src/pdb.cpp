#include <fstream>
#include <thread>
#include <sstream>
#include <fmt/printf.h>
#include "core/pdb.h"
#include "core/processing.h"

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

template <class FnOnFrame, class FnOnDone>
void parse_pdb(std::ifstream &ifs, FnOnFrame &&on_frame, FnOnDone &&on_done)
{
  std::string str;
  int frame_idx = 0;
  int num_atoms = 0;
  auto frame = pdb_frame{frame_idx++};
  while (std::getline(ifs, str)) {
    str = trim(str);
    if (str == "END") {
      std::forward<FnOnFrame>(on_frame)(std::move(frame));
      frame = pdb_frame{frame_idx++};
      continue;
    }
    std::istringstream iss(str);
    std::string start;

    iss >> start;

    if (start != "ATOM" && start != "HETATM") {
      continue;
    }

    std::string kind, atom, residue, chain;
    kind = start;
    int residue_id;
    int atom_id;
    // NOTE: We could read lines at specific offsets, but I think this is more
    // robust.
    iss >> atom_id >> atom >> residue;
    iss >> residue_id;
    if (iss.fail()) {
      iss.clear();
      iss >> chain;
      iss >> residue_id;
    }

    float x, y, z;
    iss >> x >> y >> z;
    auto entry = pdb_atom_entry{{x, y, z}, residue, atom, atom_id, residue_id, chain, kind};
    ++num_atoms;
    frame.atoms.push_back(std::move(entry));
  }

  std::forward<FnOnDone>(on_done)();

  int num_frames = frame_idx - 1;
  if (g_verbose) {
    fmt::print("Parsed pdb file: {} frames with {} atoms each\n",
                             num_frames, num_atoms / num_frames);
  }
}

void parse_pdb(std::ifstream &ifs, producer_consumer_queue &queue)
{
  auto on_frame = [&](pdb_frame &&frame) {
    // Slow down if the consumers can't keep up to avoid memory bloat.
      while (queue.frames.size_approx() > 200) {
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
      }
      queue.frames.enqueue(std::move(frame));
  };
  auto on_done = [&]{
    queue.producer_done = true;
  };
  parse_pdb(ifs, on_frame, on_done);
}

std::vector<pdb_frame> parse_pdb(const std::string &file) {
  std::vector<pdb_frame> frames;
  auto on_frame = [&](pdb_frame &&frame) {
    frames.push_back(std::move(frame));
  };
  auto on_done = []{
  };
  auto ifs = std::ifstream(file);
  parse_pdb(ifs, on_frame, on_done);

  return frames;
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