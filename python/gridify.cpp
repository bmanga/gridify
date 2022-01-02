#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fstream>

#include "core/grid.h"
#include "core/pdb.h"
#include "core/sitefind.h"
#include "core/score.h"
#include "core/pca.h"

namespace py = pybind11;
int add(int a, int b);

bool g_verbose = false;

/*
struct pdb_atom_entry {
  Point_3 pos;
  std::string residue;
  std::string atom;
  int atom_id;
  int residue_id;
  std::string chain;
  std::string kind;
};
*/
Surface_Mesh calc_surface_union_of_balls(const std::vector<Point_3> &points,
                                         double radius);
    PYBIND11_MODULE(gridify_python, m)
{
  m.doc() = "Python bindings for gridify";
  /******************* TYPES ********************/
  py::class_<Point_3>(m, "Point_3")
      .def("x", &Point_3::x)
      .def("y", &Point_3::y)
      .def("z", &Point_3::z);

  py::class_<pdb_atom_entry>(m, "pdb_atom_entry")
      .def_readonly("pos", &pdb_atom_entry::pos)
      .def_readonly("kind", &pdb_atom_entry::kind)
  .def_readonly("chain", &pdb_atom_entry::chain);

  py::class_<pdb_frame>(m, "pdb_frame")
      .def_readwrite("atoms", &pdb_frame::atoms);

  py::class_<binding_site>(m, "binding_site")
      .def_readwrite("atoms", &binding_site::atoms);

  /******************* PDB *******************/
  m.def("parse_pdb", [](const std::string &str) { return parse_pdb(str); });

  /******************* PCA *******************/
  m.def("pca_aligned_points", &pca_aligned_points);

  /******************* SITEFIND *******************/
  py::class_<ligand>(m, "ligand")
      .def_readwrite("chain", &ligand::chain)
      .def_readwrite("resid", &ligand::resid);

  m.def("discover_ligands", &discover_ligands);
  m.def("get_protein_residues_near_ligand",
        py::overload_cast<const pdb_frame &, const ligand &, double, bool>(
            &get_protein_residues_near_ligand));
  m.def("get_protein_residues_near_ligand",
        py::overload_cast<const std::vector<pdb_atom_entry> &,
                          const std::vector<pdb_atom_entry> &,
                          double, bool>(
            &get_protein_residues_near_ligand));

  /******************* GRIDIFY *******************/
  m.def("get_binding_site", &get_binding_site);
  m.def("gen_site_grid", &gen_site_grid);
  m.def("surface_only", &surface_only);

  /******************* SCORER *******************/
  py::class_<Surface_Mesh>(m, "Surface_Mesh");
  m.def(
      "gen_ligand_geometry", &gen_ligand_geometry);

  py::class_<site_ligand_stats>(m, "site_ligand_stats")
      .def_readwrite("site_volume", &site_ligand_stats::site_volume)
      .def_readwrite("ligand_volume", &site_ligand_stats::ligand_volume)
      .def_readwrite("intersection_volume",
                     &site_ligand_stats::intersection_volume)
      .def_readwrite("union_volume", &site_ligand_stats::union_volume);

  m.def("calc_grid_ligand_stats", &calc_grid_ligand_stats);
}