project(python-bindings)

find_package(pybind11 REQUIRED)

pybind11_add_module(gridify_python SHARED gridify.cpp)
target_link_libraries(gridify_python PUBLIC pybind11::pybind11 core)

