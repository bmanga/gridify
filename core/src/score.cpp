#include "core/score.h"
#include "core/radius.h"
#include "core/pca.h"

#include <vector>
#include <list>

#include <CGAL/Union_of_balls_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_union_of_balls_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_self_intersections.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

namespace {
typedef CGAL::Skin_surface_traits_3<K> Traits;
typedef CGAL::Union_of_balls_3<Traits> Union_of_balls_3;
typedef Union_of_balls_3::Bare_point Bare_point;
typedef Union_of_balls_3::Weighted_point Weighted_point;
typedef CGAL::Polyhedron_3<K> Polyhedron;

namespace PMP = CGAL::Polygon_mesh_processing;

Surface_Mesh calc_surface_union_of_balls(const std::list<Weighted_point> &l)
{
  Union_of_balls_3 union_of_balls(l.begin(), l.end());
  Polyhedron p;
  CGAL::mesh_union_of_balls_3(union_of_balls, p);

  Surface_Mesh mesh;
  CGAL::copy_face_graph(p, mesh);

  // CGAL bug workaround
  if (PMP::does_self_intersect(mesh)) {
    PMP::experimental::remove_self_intersections(mesh);
  }

  return mesh;
}

Surface_Mesh calc_surface_union_of_balls(const std::vector<Point_3> &points,
                                         double radius)
{
  std::list<Weighted_point> l;
  for (const auto &p : points) {
    l.push_front(Weighted_point(Bare_point(p), radius));
  }

  return calc_surface_union_of_balls(l);
}

Surface_Mesh calc_surface_union_of_balls(const std::vector<Point_3> &points,
                                         const std::vector<double> &radii)
{
  std::list<Weighted_point> l;
  for (int j = 0; j < points.size(); ++j) {
    const auto &p = points[j];
    double radius = radii[j];
    l.push_front(Weighted_point(Bare_point(p), radius));
  }

  return calc_surface_union_of_balls(l);
}
}


Surface_Mesh gen_ligand_geometry(const pdb_frame &f,
                         double scale_radius,
                         bool pca_align)
{
  std::vector<Point_3> points;
  std::vector<double> radii;
  const auto &radmatch = radius_matcher::get();
  for (const auto &a : f.atoms) {
    points.push_back(a.pos);
    radii.push_back(radmatch.radius(a) * scale_radius);
  }

  if (pca_align) {
    points = pca_aligned_points(points);
  }
  return calc_surface_union_of_balls(points, radii);
}
