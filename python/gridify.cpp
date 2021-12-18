#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "core/pdb.h"
#include <fstream>

namespace py = pybind11;
int add(int a, int b);

bool g_verbose = true;

PYBIND11_MODULE(gridify_python, m) {
    py::class_<pdb_frame>(m, "pdb_frame");
    m.doc() = "pybind11 example plugin"; // optional module docstring4
    m.def("parse_pdb", [](const std::string &str) { return parse_pdb(str); });
}