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
#include <CGAL/algorithm.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Triangulation_3.h>


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

std::pair<std::vector<Point_3>, bounds> parse_pdb(std::ifstream &ifs)
{
  std::vector<Point_3> points;
  bounds bounds;
  std::string str;

  while (std::getline(ifs, str)) {
    std::istringstream iss(str);
    std::string unused;
    // NOTE: We could read lines at specific offsets, but I think this is more robust.
    iss >> unused >> unused >> unused >> unused >> unused;
    float x, y, z;
    iss >> x >> y >> z;
    points.emplace_back(x, y, z);
    bounds.x.first = std::min(bounds.x.first, x);
    bounds.x.second = std::max(bounds.x.second, x);

    bounds.y.first = std::min(bounds.y.first, y);
    bounds.y.second = std::max(bounds.y.second, y);

    bounds.z.first = std::min(bounds.z.first, z);
    bounds.z.second = std::max(bounds.z.second, z);
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

int gen_grid_pdb(std::ostream &out, const Polyhedron_3 &poly, const bounds &bounds, float spacing)
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
  auto opts = cxxopts::Options("gridify", "Generates a grid of points within the convex hull of the input points");

  // clang-format off
  opts.add_options()("i,in_file", "input file", cxxopts::value<std::string>())
    ("o,out_file", "output file (stdout if omitted)",cxxopts::value<std::string>()->default_value(""))
    ("s,spacing", "grid spacing", cxxopts::value<float>()->default_value("1.0f"))
    ("v,verbose", "display extra information", cxxopts::value<bool>()->default_value("false"));
  // clang-format on
 
  opts.parse_positional("in_file");

  auto parsed_opts = opts.parse(argc, argv);

  if (parsed_opts.count("in_file") == 0) {
    std::cerr << "ERROR: You must specify input file\n";
    std::cout << opts.help();
    return -1;
  }

  auto in_file = parsed_opts["in_file"].as<std::string>();
  auto out_file = parsed_opts["out_file"].as<std::string>();
  auto spacing = parsed_opts["spacing"].as<float>();
  bool verbose = parsed_opts["verbose"].as<bool>();


  auto pdb_file = std::ifstream(in_file);
  if (!pdb_file.is_open()) {
    std::cerr << fmt::format("File '{}' does not exist", in_file);
    return -1;
  }

  auto [points, bounds] = parse_pdb(pdb_file);

  if (verbose) {
    std::cout << fmt::format("Parsed pdb file '{}': {} points\n", in_file,
                             points.size());
    std::cout << fmt::format(R"(
Bounds are: 
     MIN     MAX
x    {:.3f}  {:.3f}
y    {:.3f}  {:.3f}
z    {:.3f}  {:.3f}
)",
                             bounds.x.first, bounds.x.second, bounds.y.first,
                             bounds.y.second, bounds.z.first, bounds.z.second);
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
    std::cout << fmt::format("Generated grid to output '{}': {} grid points with spacing of {:.4f}\n",
                             out_file, num_grid_points, spacing);
  }
}
