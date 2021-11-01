#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <cxxopts.hpp>

#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_3.h>

#include "aabb_tree.hpp"
#include "common.hpp"
#include "common/processing.h"
#include "common/radius.h"
#include "common/pca.h"

#define ALL_SETTINGS            \
  X(std::string, in_file, "", "input file")       \
  X(std::string, out_file, "", "output file (stdout if omitted)")     \
  X(unsigned, jobs, 0, "number of worker threads (0 for automatic)")             \
  X(bool, verbose, false, "display extra information")              \
  X(bool, dense_packing, false, "enable dense packing for the grid point (point_radius must be > 0, spacing is set to 2 times point_radius)")        \
  X(bool, rm_atom_overlaps, false, "remove grid points that overlap with atoms")     \
  X(bool, largest_cluster_only, false, "only keep the largest cluster of close points (uses spacing to determine connectivity)") \
  X(bool, pca_align, false, "align the grid using PCA, with maximum variance on the X axis")            \
  X(bool, calculate_properties, false, "calculate additional properties and output them to stdout or <out_file>.props") \
  X(double, spacing, 1.0, "grid spacing")            \
  X(double, point_radius, 0, "if not 0, the grid points are considered spheres with the given radius")       \
  X(double, rm_lc_cutoff, 0, "if > 0, enables low connectivity grid points cutoff (uses spacing to determine connectivity)")       \
  X(std::vector<int>, site_residues, {}, "list of residues ids that make up the binding site")                                                       \
  X(double, scale_radius, 1.0, "The scale factor for atomic radii (default 1)")
#include "common/cmdline.inc"


struct site_properties {
  double volume;
  bounds site_bounds;
  bounds pca_bounds;
};

struct processed_data {
  std::vector<Point_3> out_grid_points;
  std::optional<site_properties> properties;
};

bool g_verbose = false;




std::pair<std::vector<Point_3>, bounds> get_binding_site(
    const pdb_frame &frame,
    const std::vector<int> &residues)
{
  std::vector<Point_3> points;
  bounds bounds;

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
        std::equal_range(frame.atoms.begin(), frame.atoms.end(), id, Comp{});

    for (auto it = begin; it != end; ++it) {
      const auto &pos = it->pos;
      points.push_back(pos);

      bounds.x.first = std::min<float>(bounds.x.first, pos.x());
      bounds.x.second = std::max<float>(bounds.x.second, pos.x());

      bounds.y.first = std::min<float>(bounds.y.first, pos.y());
      bounds.y.second = std::max<float>(bounds.y.second, pos.y());

      bounds.z.first = std::min<float>(bounds.z.first, pos.z());
      bounds.z.second = std::max<float>(bounds.z.second, pos.z());
    }
  }

  return {points, bounds};
}

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
                                double min_conn_score)
{
  auto overlap_score = [&](abt::aabb3d bb) {
    double cnt = 0;
    points.visit_overlaps(bb, [&](unsigned, const abt::aabb3d &bb2) {
      cnt += 0.5;
      if (std::hypot(bb2.centre.x() - bb.centre.x(),
                     bb2.centre.y() - bb.centre.y(),
                     bb2.centre.z() - bb.centre.z()) <= dist_cutoff) {
        cnt += 0.5;
      }
    });
    // Ignore the self overlap.
    return cnt - 1;
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
  visited.insert(node_id);
  std::forward<Fn>(fn)(group, node_id);

  tree.visit_overlaps(tree.get_aabb(node_id), [&](unsigned id) {
    if (visited.count(id) == 0) {
      dfs_search(id, tree, visited, group, std::forward<Fn>(fn));
    }
  });
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

struct remark_group {
  remark_group(std::FILE *out) : out(out) { fmt::print(out, "REMARK 100\n"); }

  template <class T>
  void operator()(std::string_view name, const T &value)
  {
    fmt::print(out, "REMARK 100 {} : {}\n", name, value);
  }

  void operator()(std::string_view name)
  {
    fmt::print(out, "REMARK 100 {}\n", name);
  }
  std::FILE *out;
};

void write_out_config_remarks(std::FILE *out, const config &c) {
  remark_group props(out);
  props("GENERATED FROM THE FOLLOWING CONFIGURATION");

#define X(type, name, ...)             \
  props(#name, c.name);

  ALL_SETTINGS

#undef X
}

void write_out_pdb_header(std::FILE *out, const config &c)
{
  write_out_config_remarks(out, c);
  fmt::print(
      out,
      "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           "
      "1\n");
}

void write_out_pdb_frame_remarks(std::FILE *out,
                                 const site_properties &sp,
                                 int frame_idx)
{
  const auto &[volume, site_bounds, pca_bounds] = sp;
  remark_group props(out);
  props("FRAME", std::to_string(frame_idx));
  props("VOLUME", fmt::format("{:3f}", volume));

  // SITE BOUNDS
  props("SITE_BOUNDS_X", fmt::format("{:.3f} {:.3f}", site_bounds.x.first,
                                     site_bounds.x.second));
  props("SITE_BOUNDS_Y", fmt::format("{:.3f} {:.3f}", site_bounds.y.first,
                                     site_bounds.y.second));
  props("SITE_BOUNDS_Z", fmt::format("{:.3f} {:.3f}", site_bounds.z.first,
                                     site_bounds.z.second));
  props("SITE_BOUNDS_DIMS",
        fmt::format("{:.3f} {:.3f} {:.3f}",
                    site_bounds.x.second - site_bounds.x.first,
                    site_bounds.y.second - site_bounds.y.first,
                    site_bounds.z.second - site_bounds.z.first));

  // PCA_BOUNDS
  props("PCA_BOUNDS_X",
        fmt::format("{:.3f} {:.3f}", pca_bounds.x.first, pca_bounds.x.second));
  props("PCA_BOUNDS_Y",
        fmt::format("{:.3f} {:.3f}", pca_bounds.y.first, pca_bounds.y.second));
  props("PCA_BOUNDS_Z",
        fmt::format("{:.3f} {:.3f}", pca_bounds.z.first, pca_bounds.z.second));
  props("PCA_BOUNDS_DIMS",
        fmt::format("{:.3f} {:.3f} {:.3f}",
                    pca_bounds.x.second - pca_bounds.x.first,
                    pca_bounds.y.second - pca_bounds.y.first,
                    pca_bounds.z.second - pca_bounds.z.first));
}

void write_out_pdb_frame(std::FILE *out, int frame_idx, const processed_data &data)
{
  int cnt = 0;
  int resID = 0;
  for (const auto &pt : data.out_grid_points) {
    fmt::print(
        out,
        "HETATM{:>5}  C   PTH {:>5}{:>12.3f}{:>8.3f}{:>8.3f}  0.00  0.00\n",
        ++cnt, ++resID, pt.x(), pt.y(), pt.z());
  }

  if (data.properties.has_value()) {
    write_out_pdb_frame_remarks(out, *data.properties, frame_idx);
  }
  fmt::print(out, "END\n");
}

void validate_config(config &c)
{
  if (c.in_file.empty()) {
    std::cerr << "ERROR: You must specify input file\n";
    std::exit(-1);
  }

  if (c.site_residues.empty()) {
    std::cerr << "ERROR: You must specify the atoms of the binding site\n";
    std::exit(-1);
  }

  if (c.dense_packing && c.point_radius == 0) {
    std::cerr << "With dense packing, you must specify a point_radius > 0\n";
    std::exit(-1);
  }

  if (c.dense_packing) {
    c.spacing = 2 * c.point_radius;
  }

  if (c.point_radius > 0 && c.spacing < 2 * c.point_radius) {
    std::cerr << "If point_radius is set, the spacing must be at least 2 * "
                 "point_radius\n";
    std::exit(-1);
  }

  if (c.jobs == 0) {
    c.jobs = std::thread::hardware_concurrency();
  }
}

template <>
struct fmt::formatter<std::vector<int>> {
  template <typename ParseContextT>
  typename ParseContextT::iterator parse(ParseContextT &ctx)
  {
    return ctx.begin();
  }

  template <typename FormatContextT>
  auto format(const std::vector<int> &a, FormatContextT &ctx)
  {
    return fmt::format_to(ctx.out(), "{}", fmt::join(a, " "));
  }
};

std::vector<Point_3> generate_grid_points(const config &config,
                                          const pdb_frame &frame)
{
#define X(type, name, ...) auto &name = config.name;
  ALL_SETTINGS
#undef X

  std::vector<std::vector<Point_3>> result;
  auto [points, bounds] = get_binding_site(frame, site_residues);

  if (g_verbose) {
    std::cout << fmt::format(R"(
--- Binding site info [FRAME {}] ---
Number of atoms specified: {}
Bounds are: 
     MIN     MAX
x    {:.3f}  {:.3f}
y    {:.3f}  {:.3f}
z    {:.3f}  {:.3f}
)",
                             frame.frame_idx, points.size(), bounds.x.first,
                             bounds.x.second, bounds.y.first, bounds.y.second,
                             bounds.z.first, bounds.z.second);
  }

  Surface_Mesh poly;
  CGAL::convex_hull_3(points.begin(), points.end(), poly);

  auto checker = points_checker(poly);
  if (rm_atom_overlaps) {
    auto radii_file = std::ifstream("radii.json");
    auto radmatch = radius_matcher(radii_file, config.scale_radius);
    checker.enable_check_atoms(frame, radmatch, bounds);
  }

  auto grid_points =
      gen_grid(checker, poly, bounds, spacing, point_radius, dense_packing);

  if (rm_lc_cutoff > 0 || largest_cluster_only) {
    abt::tree3d points;
    double conn_radius =
        config.spacing / 2.0 + std::numeric_limits<double>::epsilon();
    for (const auto &pt : grid_points) {
      points.insert(
          abt::aabb3d::of_sphere({pt.x(), pt.y(), pt.z()}, conn_radius));
    }

    if (rm_lc_cutoff > 0) {
      rm_low_connectivity_points(points, spacing, rm_lc_cutoff);
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

double calc_site_volume(const config &c, const std::vector<Point_3> &points)
{
  constexpr double PI = 3.14159265359;
  double sphere_vol = 4 / 3.0 * PI * std::pow((c.spacing / 2.0), 3);
  double density = PI / 6;
  if (c.dense_packing) {
    // See https://en.wikipedia.org/wiki/Close-packing_of_equal_spheres
    density = 0.74048;
  }
  double volume = sphere_vol / density * points.size();
  return volume;
}


site_properties calc_site_properties(const config &c,
                                     const std::vector<Point_3> &site_points,
                                     const std::vector<Point_3> &pca_points)
{
  return {calc_site_volume(c, site_points), get_bounds(site_points),
          get_bounds(pca_points)};
}

processed_data process_frame(const config &config, const pdb_frame &frame)
{
  fmt::print("-- Processing frame {}\n", frame.frame_idx);
  auto result = processed_data{};

  auto grid_points = generate_grid_points(config, frame);

  std::vector<Point_3> pca_points;
  std::vector<Point_3> *out_points = &grid_points;

  if (config.calculate_properties || config.pca_align) {
    pca_points = pca_aligned_points(grid_points);
  }

  if (config.pca_align) {
    out_points = &pca_points;
  }

  if (config.calculate_properties) {
    result.properties = calc_site_properties(config, grid_points, pca_points);
  }

  result.out_grid_points = std::move(*out_points);
  return result;
}

int main(int argc, char *argv[])
{
  auto config = parse_cmd_line_settings("gridify",
                                        "Generates a grid of points within the convex hull of the input points",
                                        argc, argv,
                                        [](auto &opts) {opts.parse_positional("in_file");});
  validate_config(config);

  g_verbose = config.verbose;
  if (g_verbose) {
    print_configuration(config);
  }

  auto pdb_file = std::ifstream(config.in_file);
  if (!pdb_file.is_open()) {
    std::cerr << fmt::format("File '{}' does not exist", config.in_file);
    std::exit(-1);
  }

  producer_consumer_queue queue;

  std::thread producer([&] { parse_pdb(pdb_file, queue); });

  processed_queue<processed_data> processed_frames;

  std::vector<std::thread> consumers;
  int num_consumers = config.jobs;
  for (int j = 0; j < num_consumers; ++j) {
    consumers.emplace_back(process_frame_loop(num_consumers, queue, processed_frames,
                           [&](const pdb_frame &frame) {
      return process_frame(config, frame);
    }));
  }

  std::thread writer([&] {
    auto out_file = config.out_file;
    std::FILE *out = stdout;

    if (!out_file.empty()) {
      out = std::fopen(out_file.c_str(), "w");
    }

    write_out_pdb_header(out, config);

    process_serialized_results(num_consumers, queue, processed_frames,
                               [&](int frame_idx, processed_data &data) {
        fmt::print("-- writing frame {}\n", frame_idx);
        write_out_pdb_frame(out, frame_idx, data);
      });

    if (!out_file.empty()) {
      std::fclose(out);
    }
  });

  producer.join();
  for (auto &c : consumers) {
    c.join();
  }
  writer.join();
}
