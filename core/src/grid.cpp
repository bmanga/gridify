#include "core/grid.h"

std::vector<Point_3> get_grid_points(const Surface_Mesh &poly,
                                     const points_checker &checker,
                                     const bounds &bounds,
                                     float spacing_x,
                                     float spacing_y,
                                     float spacing_z,
                                     float off_xy,
                                     float point_radius)
{
  std::vector<Point_3> result;
  auto [x, y, z] = bounds;
  double cum_off_xy = 0;
  for (auto z1 = z.first; z1 <= z.second; z1 += spacing_z) {
    for (auto y1 = y.first - cum_off_xy; y1 <= y.second; y1 += spacing_y) {
      for (auto x1 = x.first - cum_off_xy; x1 <= x.second; x1 += spacing_x) {
        auto pt = Point_3(x1, y1, z1);
        if (checker.is_inside_ch(pt) &&
            checker.has_no_atomic_clashes(pt, point_radius)) {
          result.push_back(pt);
        }
      }
    }
    cum_off_xy += off_xy;
  }
  return result;
}

std::vector<Point_3> gen_grid(const points_checker &checker,
                              const Surface_Mesh &poly,
                              const bounds &bounds,
                              float spacing,
                              float pt_radius,
                              bool dense_packing)
{
  float spacing_x = spacing;
  float spacing_y = spacing;
  float spacing_z = spacing;
  float off_xy = 0;

  if (dense_packing) {
    spacing_z = std::sqrt(2) * pt_radius;
    off_xy = pt_radius;
  }
  return get_grid_points(poly, checker, bounds, spacing_x, spacing_y, spacing_z,
                         off_xy, pt_radius);
}

void rm_low_connectivity_points(abt::tree3d &points,
                                double dist_cutoff,
                                double min_conn_score,
																double tangent_weight,
																double proximity_weight)
{
  auto tangent_addition = tangent_weight - proximity_weight;

  auto overlap_score = [&](abt::aabb3d bb) {
    double cnt = 0;
    points.visit_overlaps(bb, [&](unsigned, const abt::aabb3d &bb2) {
      cnt += proximity_weight;
      if (std::hypot(bb2.centre.x() - bb.centre.x(),
                     bb2.centre.y() - bb.centre.y(),
                     bb2.centre.z() - bb.centre.z()) <= dist_cutoff) {
        cnt += tangent_addition;
      }
    });
    // Ignore the self overlap.
    return cnt - tangent_weight;
  };

  bool stable = false;
  while (!stable) {
    stable = true;
    points.for_each([&](unsigned id, const auto &bb) {
      if (overlap_score(bb) < min_conn_score) {
        stable = false;
        points.remove(id);
      }
    });
  }
}

template <class Fn>
void dfs_search(unsigned node_id,
                const abt::tree3d &tree,
                std::unordered_set<unsigned> &visited,
                unsigned group,
                Fn &&fn)
{
  std::stack<unsigned> unvisited_ids;
  unvisited_ids.push(node_id);
  while (!unvisited_ids.empty()) {
    auto unvisited_id = unvisited_ids.top();
    unvisited_ids.pop();
    visited.insert(unvisited_id);
    tree.visit_overlaps(tree.get_aabb(unvisited_id), [&](unsigned id) {
      if (visited.count(id) == 0) {
        unvisited_ids.push(id);
      }
    });
    std::forward<Fn>(fn)(group, unvisited_id);
  }
}

template <class Fn>
void visit_connected_components(const abt::tree3d &tree, Fn &&fn)
{
  std::unordered_set<unsigned> visited;
  visited.reserve(tree.size());
  unsigned group = 0;
  tree.for_each([&](unsigned id, const auto &) {
    if (visited.count(id) == 0) {
      dfs_search(id, tree, visited, group++, std::forward<Fn>(fn));
    }
  });
}

void keep_largest_cluster_only(abt::tree3d &grid)
{
  std::vector<std::vector<unsigned>> groups(grid.size());
  visit_connected_components(grid, [&](unsigned group_id, unsigned node_id) {
    groups[group_id].push_back(node_id);
  });
  auto biggest_group_it = std::max_element(
      groups.begin(), groups.end(),
      [](const auto &a, const auto &b) { return a.size() < b.size(); });
  for (auto it = groups.begin(); it != groups.end(); ++it) {
    if (it == biggest_group_it) {
      continue;
    }
    for (unsigned id : *it) {
      grid.remove(id);
    }
  }
}