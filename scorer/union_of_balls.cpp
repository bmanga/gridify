#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Union_of_balls_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_union_of_balls_3.h>
#include <list>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_self_intersections.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3                                  Point_3;

typedef CGAL::Skin_surface_traits_3<K>                      Traits;
typedef CGAL::Union_of_balls_3<Traits>                      Union_of_balls_3;
typedef Union_of_balls_3::Bare_point                        Bare_point;
typedef Union_of_balls_3::Weighted_point                    Weighted_point;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;

namespace PMP = CGAL::Polygon_mesh_processing;

Surface_mesh calc_surface_union_of_balls(const std::list<Weighted_point> &l)
{
  Union_of_balls_3 union_of_balls(l.begin(), l.end());
  Polyhedron p;
  CGAL::mesh_union_of_balls_3(union_of_balls, p);

  Surface_mesh mesh;
  CGAL::copy_face_graph(p, mesh);

  // CGAL bug workaround
  if (PMP::does_self_intersect(mesh)) {
    PMP::experimental::remove_self_intersections(mesh);
  }

  return mesh;
}

Surface_mesh calc_surface_union_of_balls(const std::vector<Point_3> &points, double radius)
{
  std::list<Weighted_point> l;
  for (const auto &p : points) {
    l.push_front(Weighted_point(Bare_point(p), radius));
  }

  return calc_surface_union_of_balls(l);
}

Surface_mesh calc_surface_union_of_balls(const std::vector<Point_3> &points, const std::vector<double> &radii)
{
  std::list<Weighted_point> l;
  for (int j = 0; j < points.size(); ++j) {
    const auto &p = points[j];
    double radius = radii[j];
    l.push_front(Weighted_point(Bare_point(p), radius));
  }

  return calc_surface_union_of_balls(l);
}