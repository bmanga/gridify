#include "core/grid.h"
#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_3.h>

static std::vector<Point_3> get_grid_points(const Surface_Mesh &poly,
                                            const points_checker &checker,
                                            const bbox_bounds &bounds,
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

static std::vector<Point_3> gen_grid(const points_checker &checker,
                                     const Surface_Mesh &poly,
                                     const bbox_bounds &bounds,
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

template <class Fn>
static void visit_points_with_overlap_score(abt::tree3d &points,
                                            Fn &&fn,
                                            double dist_cutoff,
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

  points.for_each([&](unsigned id, const auto &bb) {
    std::forward<Fn>(fn)(id, overlap_score(bb));
  });
}


static void rm_low_connectivity_points(abt::tree3d &points,
                                       double dist_cutoff,
                                       double min_conn_score,
                                       double tangent_weight,
                                       double proximity_weight)
{
  bool stable = false;
  while (!stable) {
    stable = true;
    visit_points_with_overlap_score(
        points,
        [&](unsigned id, double score) {
          if (score < min_conn_score) {
            stable = false;
            points.remove(id);
          }
        },
        dist_cutoff, tangent_weight, proximity_weight);
  }
}

static void keep_surface_only(abt::tree3d& points, double dist_cutoff) {
  double tangent_weight = 1;
  double proximity_weight = 0;

  int cnt_removed = 0;
  std::vector<unsigned> to_remove;
  visit_points_with_overlap_score(
      points,
      [&](unsigned id, double score) {
        if (score == 6) {
          to_remove.push_back(id);
          ++cnt_removed;
        }
      },
      dist_cutoff, tangent_weight, proximity_weight);

  for (auto id : to_remove) {
    points.remove(id);
  }
  std::cout << "removed: " << cnt_removed << std::endl;
}

std::vector<Point_3> surface_only(const std::vector<Point_3>& points, double radius) {
  abt::tree3d tree;
  for (const auto &p : points) {
    tree.insert(abt::aabb3d::of_sphere({p.x(), p.y(), p.z()}, radius));
  }
  keep_surface_only(tree, radius * 2);
  std::vector<Point_3> res;
  tree.for_each([&](unsigned, abt::aabb3d bb) {
    auto c = bb.centre;
    res.emplace_back(c.x(), c.y(), c.z());
  });
  return res;
}



template <class Fn>
static void dfs_search(unsigned node_id,
                       const abt::tree3d &tree,
                       std::unordered_set<unsigned> &visited,
                       unsigned group,
                       Fn &&fn)
{
  std::unordered_set<unsigned> unvisited_ids;
  unvisited_ids.insert(node_id);
  while (!unvisited_ids.empty()) {
    auto unvisited_id = *unvisited_ids.begin();
    unvisited_ids.erase(unvisited_ids.begin());
    visited.insert(unvisited_id);
    tree.visit_overlaps(tree.get_aabb(unvisited_id), [&](unsigned id) {
      if (visited.count(id) == 0) {
        unvisited_ids.insert(id);
      }
    });
    std::forward<Fn>(fn)(group, unvisited_id);
  }
}

template <class Fn>
static void visit_connected_components(const abt::tree3d &tree, Fn &&fn)
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

static void keep_largest_cluster_only(abt::tree3d &grid)
{
  std::vector<std::vector<unsigned>> groups;
  groups.reserve(5);
  visit_connected_components(grid, [&](unsigned group_id, unsigned node_id) {
    if (group_id <= groups.size()) {
      groups.resize(group_id + 1);
    }
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

binding_site get_binding_site(const std::vector<pdb_atom_entry> &protein,
                              const std::vector<int> &residues)
{
  std::vector<pdb_atom_entry> site_atoms;
  bbox_bounds bounds;

  struct Comp {
    bool operator()(const pdb_atom_entry &s, int i) const
    {
      return s.residue_id < i;
    }
    bool operator()(int i, const pdb_atom_entry &s) const
    {
      return i < s.residue_id;
    }
  };

  for (int id : residues) {
    auto [begin, end] =
        std::equal_range(protein.begin(), protein.end(), id, Comp{});

    for (auto it = begin; it != end; ++it) {
      site_atoms.push_back(*it);

      auto pos = it->pos;
      bounds.x.first = std::min<float>(bounds.x.first, pos.x());
      bounds.x.second = std::max<float>(bounds.x.second, pos.x());

      bounds.y.first = std::min<float>(bounds.y.first, pos.y());
      bounds.y.second = std::max<float>(bounds.y.second, pos.y());

      bounds.z.first = std::min<float>(bounds.z.first, pos.z());
      bounds.z.second = std::max<float>(bounds.z.second, pos.z());
    }
  }

  if (g_verbose) {
    std::cout << fmt::format(R"(
--- Binding site info ---
Number of atoms specified: {}
Bounds are: 
     MIN     MAX
x    {:.3f}  {:.3f}
y    {:.3f}  {:.3f}
z    {:.3f}  {:.3f}
)",
                             site_atoms.size(), bounds.x.first, bounds.x.second,
                             bounds.y.first, bounds.y.second, bounds.z.first,
                             bounds.z.second);
  }

  return {site_atoms, bounds};
}

template <class It>
struct point_iterator {
  using value_type = Point_3;
  using pointer = Point_3 *;
  using reference = Point_3 &;
  using iterator_category = std::random_access_iterator_tag;
  using difference_type = std::ptrdiff_t;

  point_iterator(It it) : it(it) {}
  value_type operator*() const { return it->pos; }
  bool operator==(point_iterator rhs) { return it == rhs.it; }
  bool operator!=(point_iterator rhs) { return it != rhs.it; }
  point_iterator& operator=(const point_iterator& rhs) { it = rhs.it;
    return *this;
  }
  point_iterator &operator++()
  {
    ++it;
    return *this;
  }
  point_iterator operator++(int) {
    auto it2 = it + 1;
    return point_iterator(it2);
  }
  It it;
};
static_assert(CGAL::is_iterator<point_iterator<std::vector<pdb_atom_entry>::iterator>>::value);

std::vector<Point_3> gen_site_grid(
    const std::vector<pdb_atom_entry> &protein,
    const binding_site &site,
    bool rm_atom_overlaps,
    bool largest_cluster_only,
    bool dense_packing,
    double rm_lc_cutoff,
    double scale_radius,
    double spacing,
    double point_radius,
    double rm_lc_tangent_weight,
    double rm_lc_proximity_weight)
{
  if (site.atoms.empty()) {
    return {};
  }
  Surface_Mesh poly;
  CGAL::convex_hull_3(point_iterator(site.atoms.begin()),
                      point_iterator(site.atoms.end()), poly);

  auto checker = points_checker(poly);
  if (rm_atom_overlaps) {
    checker.enable_check_atoms(protein, scale_radius);
  }

  auto grid_points = gen_grid(checker, poly, site.bounds, spacing, point_radius,
                              dense_packing);

  if (rm_lc_cutoff > 0 || largest_cluster_only) {
    abt::tree3d points;
    double conn_radius = spacing / 2.0 + std::numeric_limits<double>::epsilon();
    for (const auto &pt : grid_points) {
      points.insert(
          abt::aabb3d::of_sphere({pt.x(), pt.y(), pt.z()}, conn_radius));
    }

    if (rm_lc_cutoff > 0) {
      rm_low_connectivity_points(points, spacing, rm_lc_cutoff,
                                 rm_lc_tangent_weight, rm_lc_proximity_weight);
    }
    if (largest_cluster_only) {
      keep_largest_cluster_only(points);
    }

    grid_points.clear();
    points.for_each([&](unsigned id, const auto &bb) {
      auto pt = bb.centre;
      grid_points.emplace_back(pt.x(), pt.y(), pt.z());
    });
  }

  return grid_points;
}
