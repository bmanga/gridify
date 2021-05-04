#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <fmt/format.h>

#include <cxxopts.hpp>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_3.h>

#include "aabb_tree.hpp"
#include "common.hpp"

bool g_verbose = false;

pdb parse_pdb(std::ifstream &ifs)
{
  pdb result;
  bounds bounds;
  std::string str;

  while (std::getline(ifs, str)) {
    if (str == "END") {
      break;
    }
    std::istringstream iss(str);
    std::string unused, atom, residue;
    int residue_id;
    int atom_id;
    // NOTE: We could read lines at specific offsets, but I think this is more
    // robust.
    iss >> unused >> atom_id >> atom >> residue;
    iss >> residue_id;
    if (iss.fail()) {
      iss.clear();
      iss >> unused;
      iss >> residue_id;
    }

    float x, y, z;
    iss >> x >> y >> z;
    auto entry = pdb_entry{{x, y, z}, residue, atom, atom_id, residue_id};
    result.atoms.push_back(std::move(entry));
  }
  return result;
}

std::pair<std::vector<Point_3>, bounds> get_binding_site(
    const pdb &pdb,
    const std::vector<int> &residues)
{
  std::vector<Point_3> points;
  bounds bounds;

  struct Comp {
    bool operator()(const pdb_entry &s, int i) const
    {
      return s.residue_id < i;
    }
    bool operator()(int i, const pdb_entry &s) const
    {
      return i < s.residue_id;
    }
  };

  for (int id : residues) {
    auto [begin, end] =
        std::equal_range(pdb.atoms.begin(), pdb.atoms.end(), id, Comp{});

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
    spacing_x = 2 * pt_radius;
    spacing_y = 2 * pt_radius;
    spacing_z = std::sqrt(2) * pt_radius;
    off_xy = pt_radius;
  }
  return get_grid_points(poly, checker, bounds, spacing_x, spacing_y, spacing_z,
                         off_xy, pt_radius);
}

void rm_low_connectivity_points(std::vector<Point_3> &grid,
                                double dist_cutoff,
                                double min_conn_score)
{
  abt::tree3d points;
  unsigned idx = 0;
  for (const auto &pt : grid) {
    points.insert(idx++,
                  abt::aabb3d::of_sphere({pt.x(), pt.y(), pt.z()}, dist_cutoff));
  }
  auto count_overlaps = [&](abt::point3d pt) {
    float cnt = 0;
    points.visit_overlaps(
        pt,
        [&](unsigned, const abt::aabb3d &bb) {
          cnt += 0.5;
          if (std::hypot(pt.x() - bb.centre.x(), pt.y() - bb.centre.y(),
                         pt.z() - bb.centre.z()) <= dist_cutoff) {
            cnt += 0.5;
          }
        });
    return cnt < min_conn_score;
  };
  bool stable = false;
  while (!stable) {
    stable = true;

    points.for_each([&](unsigned id, const auto &bb) {
      if (count_overlaps(bb.centre)) {
        stable = false;
        points.remove(id);
      }
    });
  }

  grid.clear();
  points.for_each([&](unsigned id, const auto &bb) {
    auto pt = bb.centre;
    grid.emplace_back(pt.x(), pt.y(), pt.z());
  });
}

void write_out_pdb(std::ostream &out, const std::vector<Point_3> &points)
{
  out << "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           "
         "1\n";
  int cnt = 0;
  int resID = 0;
  for (const auto &pt : points) {
    out << fmt::format(
        "HETATM{:>5}  C   PTH {:>5}{:>12.3f}{:>8.3f}{:>8.3f}  0.00  0.00\n",
        ++cnt, ++resID, pt.x(), pt.y(), pt.z());
  }
  out << "END\n";
}

int main(int argc, char *argv[])
{
  auto opts = cxxopts::Options(
      "gridify",
      "Generates a grid of points within the convex hull of the input points");

  // clang-format off
  opts.add_options()("i,in_file", "input file", cxxopts::value<std::string>())
    ("o,out_file", "output file (stdout if omitted)",cxxopts::value<std::string>()->default_value(""))
    ("s,spacing", "grid spacing", cxxopts::value<float>()->default_value("1.0f"))
    ("v,verbose", "display extra information", cxxopts::value<bool>()->default_value("false"))
    ("r,site_residues", "list of residues ids that make up the binding site", cxxopts::value<std::vector<int>>())
    ("rm_atom_overlaps", "remove grid points that overlap with atoms", cxxopts::value<bool>()->default_value("false"))
    ("point_radius", "if not 0, the grid points are considered spheres with the given radius", cxxopts::value<double>()->default_value("0"))
    ("dense_packing", "enable dense packing for the grid point (point_radius must be > 0, spacing is ignored)", cxxopts::value<bool>()->default_value("false"))
    ("rm_lc_cutoff", "if > 0, enables low connectivity grid points cutoff (uses spacing to determine connectivity)", cxxopts::value<double>()->default_value("0"));
  // clang-format on

  opts.parse_positional("in_file");

  auto parsed_opts = opts.parse(argc, argv);

  if (parsed_opts.count("in_file") == 0) {
    std::cerr << "ERROR: You must specify input file\n";
    std::cout << opts.help();
    return -1;
  }

  if (parsed_opts.count("site_residues") == 0) {
    std::cerr << "ERROR: You must specify the atoms of the binding site\n";
    std::cout << opts.help();
    return -1;
  }

  auto in_file = parsed_opts["in_file"].as<std::string>();
  auto out_file = parsed_opts["out_file"].as<std::string>();
  auto spacing = parsed_opts["spacing"].as<float>();
  g_verbose = parsed_opts["verbose"].as<bool>();
  auto residues = parsed_opts["site_residues"].as<std::vector<int>>();
  bool check_atoms = parsed_opts["rm_atom_overlaps"].as<bool>();
  auto point_radius = parsed_opts["point_radius"].as<double>();
  bool dense_packing = parsed_opts["dense_packing"].as<bool>();
  auto rm_lc_cutoff = parsed_opts["rm_lc_cutoff"].as<double>();

  if (dense_packing && point_radius == 0) {
    std::cerr << "With dense packing, you must specify a point_radius > 0\n";
    return -1;
  }

  auto pdb_file = std::ifstream(in_file);
  if (!pdb_file.is_open()) {
    std::cerr << fmt::format("File '{}' does not exist", in_file);
    return -1;
  }

  auto pdb = parse_pdb(pdb_file);
  auto [points, bounds] = get_binding_site(pdb, residues);

  if (g_verbose) {
    std::cout << fmt::format("Parsed pdb file '{}': {} atoms\n", in_file,
                             pdb.atoms.size());
    std::cout << fmt::format("Binding site info: Num residues {}, : {} atoms\n",
                             in_file, pdb.atoms.size());
    std::cout << fmt::format(R"(
--- Binding site info ---
Number of atoms specified: {}
Bounds are: 
     MIN     MAX
x    {:.3f}  {:.3f}
y    {:.3f}  {:.3f}
z    {:.3f}  {:.3f}
)",
                             points.size(), bounds.x.first, bounds.x.second,
                             bounds.y.first, bounds.y.second, bounds.z.first,
                             bounds.z.second);
  }

  Surface_Mesh poly;
  CGAL::convex_hull_3(points.begin(), points.end(), poly);

  auto checker = points_checker(poly);
  if (check_atoms) {
    auto radii_file = std::ifstream("radii.json");
    auto radmatch = radius_matcher(radii_file);
    checker.enable_check_atoms(pdb, radmatch, bounds);
  }

  auto grid_points =
      gen_grid(checker, poly, bounds, spacing, point_radius, dense_packing);

  if (rm_lc_cutoff > 0) {
    rm_low_connectivity_points(grid_points, spacing, rm_lc_cutoff);
  }

  if (out_file.empty()) {
    write_out_pdb(std::cout, grid_points);
  }
  else {
    auto ofs = std::ofstream(out_file);
    write_out_pdb(ofs, grid_points);
  }

  if (g_verbose) {
    std::cout << fmt::format(
        "Generated grid to output '{}': {} grid points with spacing of "
        "{:.4f}\n",
        out_file, grid_points.size(), spacing);
  }
}
