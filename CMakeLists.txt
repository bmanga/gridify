# CMakeList.txt: file di progetto CMake di primo livello. Eseguire la configurazione globale
# e includere qui i sottoprogetti.
#
cmake_minimum_required (VERSION 3.12)

set(CMAKE_CXX_STANDARD 20)


project ("gridify")

list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})


# Includere i sottoprogetti.
add_subdirectory(common)
add_subdirectory ("gridify")
