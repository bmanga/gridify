#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <fmt/format.h>

#include <cxxopts.hpp>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron_3 = CGAL::Polyhedron_3<K>;
using Point_3 = K::Point_3;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Polyhedron_3>;
using Traits = CGAL::AABB_traits<K, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;
using Point_inside = CGAL::Side_of_triangle_mesh<Polyhedron_3, K>;

struct bounds {
  static constexpr std::pair<float, float> INIT = {
      std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest()};

  std::pair<float, float> x = INIT, y = INIT, z = INIT;
};

struct points_checker {
  points_checker(const Polyhedron_3 &poly)
      : tree(faces(poly).first, faces(poly).second, poly), inside_tester(tree)
  {
    tree.accelerate_distance_queries();
  }

  bool is_inside(Point_3 point) const
  {
    return inside_tester(point) != CGAL::ON_UNBOUNDED_SIDE;
  }

  Tree tree;
  Point_inside inside_tester;
};

struct pdb_entry {
  Point_3 pos;
  std::string residue;
  std::string atom;
  int residue_id;
};

struct pdb {
  std::vector<pdb_entry> atoms;
};

pdb parse_pdb(std::ifstream &ifs)
{
  pdb result;
  bounds bounds;
  std::string str;

  while (std::getline(ifs, str)) {
    std::istringstream iss(str);
    std::string unused, atom, residue;
    int residue_id;
    // NOTE: We could read lines at specific offsets, but I think this is more
    // robust.
    iss >> unused >> unused >> atom >> residue;
    iss >> residue_id;
    if (iss.fail()) {
      iss.clear();
      iss >> unused;
      iss >> residue_id;
    }

    float x, y, z;
    iss >> x >> y >> z;
    auto entry = pdb_entry{{x, y, z}, residue, atom, residue_id};
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
    bool operator()(const pdb_entry &s, int i) const { return s.residue_id < i; }
    bool operator()(int i, const pdb_entry &s) const { return i < s.residue_id; }
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

template <class CBFun>
void for_point_in_poly(const Polyhedron_3 &poly,
                       const bounds &bounds,
                       float spacing,
                       CBFun &&callback)
{
  auto checker = points_checker(poly);
  for (auto [z1, z2] = bounds.z; z1 <= z2; z1 += spacing) {
    for (auto [y1, y2] = bounds.y; y1 <= y2; y1 += spacing) {
      for (auto [x1, x2] = bounds.x; x1 <= x2; x1 += spacing) {
        auto pt = Point_3(x1, y1, z1);
        if (checker.is_inside(pt)) {
          callback(pt);
        }
      }
    }
  }
}

int gen_grid_pdb(std::ostream &out,
                 const Polyhedron_3 &poly,
                 const bounds &bounds,
                 float spacing)
{
  out << "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           "
         "1\n";

  int cnt = 0;
  int resID = 0;

  for_point_in_poly(poly, bounds, spacing, [&](const Point_3 &pt) {
    out << fmt::format(
        "HETATM{:>5}  C   PTH {:>5}{:>12.3f}{:>8.3f}{:>8.3f}  0.00  0.00\n",
        ++cnt, ++resID, pt.x(), pt.y(), pt.z());
  });

  out << "END\n";
  return cnt;
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
    ("r,site_residues", "list of atom ids that make up the binding site", cxxopts::value<std::vector<int>>());
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
  bool verbose = parsed_opts["verbose"].as<bool>();
  auto residues = parsed_opts["site_residues"].as<std::vector<int>>();

  auto pdb_file = std::ifstream(in_file);
  if (!pdb_file.is_open()) {
    std::cerr << fmt::format("File '{}' does not exist", in_file);
    return -1;
  }

  auto pdb = parse_pdb(pdb_file);
  auto [points, bounds] = get_binding_site(pdb, residues);

  if (verbose) {
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

  Polyhedron_3 poly;
  CGAL::convex_hull_3(points.begin(), points.end(), poly);

  int num_grid_points = 0;

  if (out_file.empty()) {
    num_grid_points = gen_grid_pdb(std::cout, poly, bounds, spacing);
  }
  else {
    auto ofs = std::ofstream(out_file);
    num_grid_points = gen_grid_pdb(ofs, poly, bounds, spacing);
  }

  if (verbose) {
    std::cout << fmt::format(
        "Generated grid to output '{}': {} grid points with spacing of "
        "{:.4f}\n",
        out_file, num_grid_points, spacing);
  }
}
