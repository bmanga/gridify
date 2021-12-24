#ifndef GRIDIFY_RADIUS_H
#define GRIDIFY_RADIUS_H

#include <fstream>
#include <nlohmann/json.hpp>
#include <fmt/printf.h>

using json = nlohmann::json;

extern bool g_verbose;

template <class... Ts>
void vlog(Ts &&...ts)
{
  if (g_verbose) {
    fmt::print(std::forward<Ts>(ts)...);
  }
}

struct radius_matcher {
 private:
  struct private_token {
  };

 public:
  static const radius_matcher& get() {
    static auto rm = radius_matcher(private_token{});
    return rm;
  }
  radius_matcher(private_token)
  {
    auto descriptor = std::ifstream("radii.json");
    if (!descriptor.good()) {
      fmt::print("WARN: no radius file found (radii.json).\n");
      radii["default"] = 2;
    }
    else {
      descriptor >> radii;
    }
  }

  double radius(pdb_atom_entry entry) const
  {
    const auto &resid = entry.residue;
    auto atom = entry.atom;
    while (std::isdigit(atom.back())) {
      atom.pop_back();
    }

    auto vlog_radius = [&entry, &atom](const char *where, double radius) {
      vlog(
          "RADIUS - value for atom {} ({}.{}, using {}.{}) "
          "= {}\n",
          entry.atom_id, entry.residue, entry.atom, where, atom, radius);
    };

    if (radii.contains(resid)) {
      if (radii[resid].contains(atom)) {
        auto radius = radii[resid][atom].get<double>();
        vlog_radius("residue", radius);
        return radius;
      }
    }
    while (atom.size() > 0) {
      if (radii["defaults"].contains(atom)) {
        auto radius = radii["defaults"][atom].get<double>();
        vlog_radius("defaults", radius);
        return radius;
      }
      atom.pop_back();
    }

    auto def = radii["default"].get<double>();
    vlog(
        "WARN: RADIUS - no suitable radius found for atom {} ({}.{}): Using "
        "specified "
        "default of {}.\n",
        entry.atom_id, entry.atom, entry.residue, def);
    return def;
  }

  json radii;
};

#endif  // GRIDIFY_RADIUS_H
