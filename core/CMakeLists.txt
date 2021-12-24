add_library(core src/pdb.cpp src/pca.cpp src/grid.cpp "src/sitefind.cpp" "src/score.cpp" "include/core/sitefind.h" "include/core/score.h")
target_include_directories(core PUBLIC include)
find_package(CGAL REQUIRED)
find_package(concurrentqueue REQUIRED)
find_package(fmt REQUIRED)
find_package(LibCDS REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(cxxopts REQUIRED)




target_link_libraries(core PUBLIC CGAL::CGAL concurrentqueue::concurrentqueue fmt::fmt LibCDS::LibCDS nlohmann_json::nlohmann_json yaml-cpp::yaml-cpp cxxopts::cxxopts)