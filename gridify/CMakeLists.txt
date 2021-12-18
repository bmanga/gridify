
cmake_minimum_required (VERSION 3.12)

find_package(CGAL REQUIRED)
find_package(fmt REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(concurrentqueue REQUIRED)

add_executable (gridify main.cpp common.hpp)
target_link_libraries(gridify PUBLIC core CGAL::CGAL fmt::fmt Eigen3::Eigen3
concurrentqueue::concurrentqueue)
