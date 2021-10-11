#ifndef GRIDIFY_COMMON_H
#define GRIDIFY_COMMON_H

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;
using Surface_Mesh = CGAL::Surface_mesh<Point_3>;

using Primitive = CGAL::AABB_face_graph_triangle_primitive<Surface_Mesh>;
using Traits = CGAL::AABB_traits<K, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;
using Point_inside = CGAL::Side_of_triangle_mesh<Surface_Mesh, K>;

extern bool g_verbose;




#endif  // GRIDIFY_COMMON_H
