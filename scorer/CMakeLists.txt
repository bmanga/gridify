add_executable(scorer scorer.cpp)
find_package(LibCDS REQUIRED)
find_package(cxxopts REQUIRED)
target_link_libraries(scorer PUBLIC core LibCDS::LibCDS cxxopts::cxxopts)