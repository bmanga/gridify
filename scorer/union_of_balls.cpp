#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Union_of_balls_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_union_of_balls_3.h>
#include <list>
#include <CGAL/Surface_mesh.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3                                  Point_3;

typedef CGAL::Skin_surface_traits_3<K>                      Traits;
typedef CGAL::Union_of_balls_3<Traits>                      Union_of_balls_3;
typedef Union_of_balls_3::Bare_point                        Bare_point;
typedef Union_of_balls_3::Weighted_point                    Weighted_point;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;

Surface_mesh calc_surface_union_of_balls(const std::vector<Point_3> &points, int radius)
{
  std::list<Weighted_point> l;
  for (const auto &p : points) {
    l.push_front(Weighted_point(Bare_point(p), radius));
  }

  Union_of_balls_3 union_of_balls(l.begin(), l.end());
  Polyhedron p;
  CGAL::mesh_union_of_balls_3(union_of_balls, p);

  Surface_mesh mesh;
  CGAL::copy_face_graph(p, mesh);

  std::ofstream f ("balls.off");
  f << mesh;
  f.close ();
  return mesh;
}