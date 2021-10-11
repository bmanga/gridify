add_library(common src/pdb.cpp)
target_include_directories(common PUBLIC include)
find_package(CGAL REQUIRED)
find_package(concurrentqueue REQUIRED)
find_package(fmt REQUIRED)


target_link_libraries(common PUBLIC CGAL::CGAL concurrentqueue::concurrentqueue fmt::fmt)