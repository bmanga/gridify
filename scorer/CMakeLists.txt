add_executable(scorer scorer.cpp alpha_shape.cpp union_of_balls.cpp)
find_package(LibCDS REQUIRED)

target_link_libraries(scorer PUBLIC common LibCDS::LibCDS)