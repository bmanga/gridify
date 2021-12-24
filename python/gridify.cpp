#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fstream>

#include "core/grid.h"
#include "core/pdb.h"
#include "core/sitefind.h"

namespace py = pybind11;
int add(int a, int b);

bool g_verbose = true;

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

PYBIND11_MODULE(gridify_python, m)
{
  m.doc() = "Python bindings for gridify";
  /******************* TYPES ********************/
  py::class_<Point_3>(m, "Point_3")
      .def("x", &Point_3::x)
      .def("y", &Point_3::y)
      .def("z", &Point_3::z);

  py::class_<pdb_atom_entry>(m, "pdb_atom_entry")
      .def_readonly("pos", &pdb_atom_entry::pos);

  py::class_<pdb_frame>(m, "pdb_frame")
      .def_readwrite("atoms", &pdb_frame::atoms);

  /******************* PDB *******************/
  m.def("parse_pdb", [](const std::string &str) { return parse_pdb(str); });

  /******************* SITEFIND *******************/
  m.def("discover_ligands", &discover_ligands);
  m.def("get_protein_residues_near_ligand", &get_protein_residues_near_ligand);

  /******************* GRIDIFY *******************/
  m.def("generate_grid_points", &generate_grid_points);
}