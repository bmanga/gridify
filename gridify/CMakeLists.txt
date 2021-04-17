
cmake_minimum_required (VERSION 3.12)
set(CMAKE_CXX_STANDARD 20)

find_package(CGAL REQUIRED)
find_package(fmt REQUIRED)
find_package(cxxopts REQUIRED)
add_executable (gridify "main.cpp" )
target_link_libraries(gridify PUBLIC CGAL::CGAL fmt::fmt cxxopts::cxxopts)
