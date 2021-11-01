#pragma once

#include <cassert>
#include <optional>
#include <string>
#include <vector>

#include <fmt/format.h>

#include "aabb_tree.hpp"

#include <common/common.h>
#include <common/pdb.h>
#include "common/radius.h"

extern bool g_verbose;

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

    auto center = abt::point3d(point.x(), point.y(), point.z());
    bool no_intersections = true;

    auto visitor = [&](unsigned id, const abt::aabb3d &bb) {
      auto radius = bb.upperBound.x() - bb.centre.x();
      auto dist_to_center = CGAL::squared_distance(
          point, Point_3(bb.centre.x(), bb.centre.y(), bb.centre.z()));
      if (dist_to_center <= radius * radius) {
        no_intersections = false;
        return abt::visit_stop;
      }
      else {
        return abt::visit_continue;
      }
    };
    if (radius == 0) {
      atoms_tree.visit_overlaps(center, visitor);
    }
    else {
      atoms_tree.visit_overlaps(abt::aabb3d::of_sphere(center, radius), visitor);
    }

    return no_intersections;
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
