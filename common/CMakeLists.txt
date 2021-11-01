add_library(common src/pdb.cpp src/pca.cpp)
target_include_directories(common PUBLIC include)
find_package(CGAL REQUIRED)
find_package(concurrentqueue REQUIRED)
find_package(fmt REQUIRED)
find_package(LibCDS REQUIRED)
find_package(nlohmann_json REQUIRED)



target_link_libraries(common PUBLIC CGAL::CGAL concurrentqueue::concurrentqueue fmt::fmt LibCDS::LibCDS nlohmann_json::nlohmann_json)