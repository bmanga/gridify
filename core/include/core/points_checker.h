#include "core/common.h"
#include "core/aabb_tree.hpp"
#include "core/pdb.h"
#include "core/radius.h"

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

  void enable_check_atoms(const std::vector<pdb_atom_entry> &atoms,
                          double radius_scale_factor = 1)
  {
    abt::tree3d atoms_tree;
    int cnt = 0;

    const auto &radius_matcher = radius_matcher::get();
    for (const auto &atom : atoms) {
      // TODO: use csv.
      float radius = radius_matcher.radius(atom) * radius_scale_factor;
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
