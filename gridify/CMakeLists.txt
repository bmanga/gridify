
cmake_minimum_required (VERSION 3.12)

find_package(CGAL REQUIRED)
find_package(fmt REQUIRED)
find_package(cxxopts REQUIRED)
find_package(cxxopts REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(concurrentqueue REQUIRED)

add_executable (gridify main.cpp common.hpp principal_lengths.cpp)
target_link_libraries(gridify PUBLIC common CGAL::CGAL fmt::fmt cxxopts::cxxopts nlohmann_json::nlohmann_json yaml-cpp::yaml-cpp Eigen3::Eigen3
concurrentqueue::concurrentqueue)

if (MSVC)
	add_definitions(/bigobj)
endif()
