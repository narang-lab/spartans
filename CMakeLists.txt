project(TetrahedronTransport)
cmake_minimum_required(VERSION 3.15...3.19)

#set(CMAKE_BUILD_TYPE "Debug")

if(SKBUILD)

  #pybind11 cmake dir
  execute_process(
    COMMAND
      "${PYTHON_EXECUTABLE}" -c
      "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE _tmp_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
  list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
  
  #mpi4py include dir
  execute_process(
    COMMAND
      "${PYTHON_EXECUTABLE}" -c
      "import mpi4py as m; print(m.get_include())"
    OUTPUT_VARIABLE mpi4py_INCLUDE_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)

endif()

find_package(pybind11 CONFIG REQUIRED)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS} -O3 -std=c++0x -flto")
set(EXTERNAL_LIBRARIES ${HDF5_LIBRARIES} ${MPI_CXX_LIBRARIES})

set(CMAKE_AR  "gcc-ar")
set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
set(CMAKE_CXX_ARCHIVE_FINISH   true)

add_library(TetrahedronTransport src/Tetrahedron.cpp src/SubTet.cpp src/TriQuad.cpp src/MeshTransport.cpp src/OverlapInv.cpp src/TetLookup.cpp src/Util.cpp)
target_link_libraries(TetrahedronTransport ${EXTERNAL_LIBRARIES})

# Create standalone executable
add_executable(transport-standalone src/transport-standalone.cpp)
target_link_libraries(transport-standalone TetrahedronTransport ${EXTERNAL_LIBRARIES})

# Create py module
target_sources(TetrahedronTransport PRIVATE src/Transport.cpp)
pybind11_add_module(_core MODULE src/pyb11mpi.cpp)

target_link_libraries(_core PRIVATE TetrahedronTransport ${EXTERNAL_LIBRARIES})
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_include_directories(_core SYSTEM PRIVATE ${mpi4py_INCLUDE_DIRS})

install(TARGETS _core DESTINATION .)
