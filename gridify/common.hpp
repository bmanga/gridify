#pragma once

#include <cassert>
#include <optional>
#include <string>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <nlohmann/json.hpp>

#include "aabb_tree.hpp"

using json = nlohmann::json;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;
using Surface_Mesh = CGAL::Surface_mesh<Point_3>;

using Primitive = CGAL::AABB_face_graph_triangle_primitive<Surface_Mesh>;
using Traits = CGAL::AABB_traits<K, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;
using Point_inside = CGAL::Side_of_triangle_mesh<Surface_Mesh, K>;

extern bool g_verbose;

template <class... Ts>
void vlog(Ts &&...ts)
{
  if (g_verbose) {
    fmt::print(std::forward<Ts>(ts)...);
  }
}

struct pdb_entry {
  Point_3 pos;
  std::string residue;
  std::string atom;
  int atom_id;
  int residue_id;
};

struct pdb {
  std::vector<pdb_entry> atoms;
};

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

  double radius(const pdb_entry &entry) const
  {
    const auto &resid = entry.residue;
    auto atom = entry.atom;
    while (std::isdigit(atom.back())) {
      atom.pop_back();
    }

    auto vlog_radius = [&entry, &atom](const char *where, double radius) {
      vlog(
          "RADIUS - using {} value for atom {} ({}.{}, using atomname {}) "
          "= {}\n",
          where, entry.atom_id, entry.residue, entry.atom, atom, radius);
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

  bool is_inside(Point_3 point) const
  {
    bool inside_poly = inside_tester(point) != CGAL::ON_UNBOUNDED_SIDE;
    if (inside_poly && opt_atoms_tree.has_value()) {
      return has_no_atomic_clashes(point);
    }
    return inside_poly;
  }

  bool has_no_atomic_clashes(Point_3 point) const
  {
    auto intersections =
        opt_atoms_tree->query(abt::point3d(point.x(), point.y(), point.z()));
    if (intersections.empty()) {
      return true;
    }
    for (const auto id : intersections) {
      auto bb = opt_atoms_tree->getAABB(id);
      auto radius = bb.upperBound.x() - bb.centre.x();
      auto dist_to_center = CGAL::squared_distance(
          point, Point_3(bb.centre.x(), bb.centre.y(), bb.centre.z()));

      if (dist_to_center <= radius * radius) {
        return false;
      }
    }
    return true;
  }

  void enable_check_atoms(const pdb &pdb,
                          const radius_matcher &radius_matcher,
                          const bounds &bs_bounds)
  {
    abt::tree3d atoms_tree;
    int cnt = 0;

    for (const auto &atom : pdb.atoms) {
      // TODO: use csv.
      float radius = radius_matcher.radius(atom);
      if (radius == 0) {
        // Radius of 0 means the atom can be ignored.
        continue;
      }
      auto pos = atom.pos;
      if (bs_bounds.is_inside(pos) && is_inside(pos)) {
        auto pos = atom.pos;
        atoms_tree.insertParticle(++cnt, {pos.x(), pos.y(), pos.z()}, radius);
      }
    }
    opt_atoms_tree = std::move(atoms_tree);
  }

  Tree tree;
  std::optional<abt::tree3d> opt_atoms_tree;
  Point_inside inside_tester;
};
