#pragma once

#include <cassert>
#include <optional>
#include <string>
#include <vector>


#include <nlohmann/json.hpp>
#include <fmt/format.h>
#include <yaml-cpp/yaml.h>

#include "aabb_tree.hpp"

#include <common/common.h>
#include <common/pdb.h>

using json = nlohmann::json;



extern bool g_verbose;

template <class... Ts>
void vlog(Ts &&...ts)
{
  if (g_verbose) {
    fmt::print(std::forward<Ts>(ts)...);
  }
}

struct bounds {
  static constexpr std::pair<float, float> INIT = {
      std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest()};

  std::pair<float, float> x = INIT, y = INIT, z = INIT;

  bool is_inside(Point_3 point) const
  {
    return point.x() <= x.second && point.x() >= x.first &&
           point.y() <= y.second && point.y() >= y.first &&
           point.z() <= z.second && point.z() >= z.first;
  }
};

struct radius_matcher {
  radius_matcher(std::istream &descriptor)
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

struct points_checker {
  points_checker(const Surface_Mesh &poly)
      : tree(faces(poly).first, faces(poly).second, poly), inside_tester(tree)
  {
    tree.accelerate_distance_queries();
  }

  bool is_inside_ch(Point_3 point) const
  {
    return inside_tester(point) != CGAL::ON_UNBOUNDED_SIDE;
  }

  bool has_no_atomic_clashes(Point_3 point, float radius) const
  {
    if (!opt_atoms_tree.has_value()) {
      return true;
    }
    auto &atoms_tree = opt_atoms_tree.value();

    std::vector<unsigned> intersections;

    auto center = abt::point3d(point.x(), point.y(), point.z());
    if (radius == 0) {
      intersections = atoms_tree.get_overlaps(center);
    }
    else {
      intersections =
          atoms_tree.get_overlaps(abt::aabb3d::of_sphere(center, radius));
    }
    if (intersections.empty()) {
      return true;
    }
    for (const auto id : intersections) {
      auto bb = opt_atoms_tree->get_aabb(id);
      auto radius = bb.upperBound.x() - bb.centre.x();
      auto dist_to_center = CGAL::squared_distance(
          point, Point_3(bb.centre.x(), bb.centre.y(), bb.centre.z()));

      if (dist_to_center <= radius * radius) {
        return false;
      }
    }
    return true;
  }

  void enable_check_atoms(const pdb_frame &frame,
                          const radius_matcher &radius_matcher,
                          const bounds &bs_bounds)
  {
    abt::tree3d atoms_tree;
    int cnt = 0;

    for (const auto &atom : frame.atoms) {
      // TODO: use csv.
      float radius = radius_matcher.radius(atom);
      if (radius == 0) {
        // Radius of 0 means the atom can be ignored.
        continue;
      }
      auto pos = atom.pos;
      atoms_tree.insert(abt::aabb3d::of_sphere({pos.x(), pos.y(), pos.z()}, radius));
    }
    opt_atoms_tree = std::move(atoms_tree);
  }

  Tree tree;
  std::optional<abt::tree3d> opt_atoms_tree;
  Point_inside inside_tester;
};
