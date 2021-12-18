# CMakeList.txt: file di progetto CMake di primo livello. Eseguire la configurazione globale
# e includere qui i sottoprogetti.
#
cmake_minimum_required (VERSION 3.12)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

project ("gridify")

list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})


if (MSVC)
	add_definitions(/bigobj)
endif()

# Includere i sottoprogetti.
add_subdirectory(core)
add_subdirectory ("gridify")
add_subdirectory(scorer)
add_subdirectory(sitefind)

