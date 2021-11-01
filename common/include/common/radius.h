#ifndef GRIDIFY_RADIUS_H
#define GRIDIFY_RADIUS_H

#include <istream>
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
  radius_matcher(std::istream &descriptor, double scale_factor):
        scale_factor(scale_factor)
  {
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

    auto vlog_radius = [&entry, &atom, scale_factor=this->scale_factor](const char *where, double radius) {
      vlog(
          "RADIUS - value for atom {} ({}.{}, using {}.{}, sf {}) "
          "= {}\n",
          entry.atom_id, entry.residue, entry.atom, where, atom, scale_factor, radius);
    };

    if (radii.contains(resid)) {
      if (radii[resid].contains(atom)) {
        auto radius = radii[resid][atom].get<double>() * scale_factor;
        vlog_radius("residue", radius);
        return radius;
      }
    }
    while (atom.size() > 0) {
      if (radii["defaults"].contains(atom)) {
        auto radius = radii["defaults"][atom].get<double>() * scale_factor;
        vlog_radius("defaults", radius);
        return radius;
      }
      atom.pop_back();
    }

    auto def = radii["default"].get<double>() * scale_factor;
    vlog(
        "WARN: RADIUS - no suitable radius found for atom {} ({}.{}): Using "
        "specified "
        "default of {} (sf {}).\n",
        entry.atom_id, entry.atom, entry.residue, def, scale_factor);
    return def;
  }

  json radii;
  double scale_factor;
};

#endif  // GRIDIFY_RADIUS_H
