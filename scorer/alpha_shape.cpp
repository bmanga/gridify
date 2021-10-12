#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/remove_outliers.h>

#include <CGAL/Nef_polyhedron_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>

#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>

#include <CGAL/Point_set_3.h>

#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>

#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/Point_set_3/IO.h>

#include <cassert>
#include <fstream>
#include <list>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Gt;
typedef CGAL::Alpha_shape_vertex_base_3<Gt>          Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>            Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>  Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>         Alpha_shape_3;
typedef Gt::Point_3                                  Point_3;
typedef Alpha_shape_3::Alpha_iterator                Alpha_iterator;
typedef Alpha_shape_3::Vertex_handle			    Vertex_handle;

typedef Gt::Vector_3 Vector_3;
typedef Gt::Sphere_3 Sphere_3;
typedef CGAL::Point_set_3<Point_3, Vector_3> Point_set;


using Nef_polyhedron_3 = CGAL::Nef_polyhedron_3<Gt>;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;


CGAL::Surface_mesh<Point_3> calc_alpha_shape_geometries(const std::vector<Point_3> &vertices)
{
  // compute alpha shape
  Alpha_shape_3 as(vertices.begin(), vertices.end());
  // find optimal alpha value
  Alpha_iterator opt = as.find_optimal_alpha(1);
  std::cout << "Optimal alpha value to get one connected component is " << *opt << std::endl;
  as.set_alpha(*opt);
  assert(as.number_of_solid_components() == 1);

  std::list<Vertex_handle> out_v;

  as.get_alpha_shape_vertices(std::back_inserter(out_v), Alpha_shape_3::REGULAR);

  Point_set points;

  for (const auto &o : vertices) {
    points.insert(o);
  }

  double spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag> (points, 6);

  typedef std::array<std::size_t, 3> Facet; // Triple of indices
  std::vector<Facet> facets;
  // The function is called using directly the points raw iterators
  CGAL::advancing_front_surface_reconstruction(points.points().begin(),
                                               points.points().end(),
                                               std::back_inserter(facets));

  // copy points for random access
  std::vector<Point_3> vertices2;
  vertices2.reserve (points.size());
  std::copy (points.points().begin(), points.points().end(), std::back_inserter (vertices2));
  CGAL::Surface_mesh<Point_3> output_mesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh (vertices2, facets, output_mesh);
  std::ofstream f ("out.off");
  f << output_mesh;
  f.close ();

  return output_mesh;
}
