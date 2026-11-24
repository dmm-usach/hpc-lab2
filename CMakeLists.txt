cmake_minimum_required(VERSION 3.20)
project(hpc_lab2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED COMPONENTS CXX)

# BSPonMPI no se instala con un paquete del sistema: se busca la librería
# y el header, y se puede sobreescribir con -DBSPONMPI_ROOT=...
find_library(BSPONMPI_LIBRARY NAMES bsponmpi)
find_path(BSPONMPI_INCLUDE_DIR NAMES bsp.h)

if(NOT BSPONMPI_LIBRARY OR NOT BSPONMPI_INCLUDE_DIR)
    message(FATAL_ERROR
        "No se encontró BSPonMPI. Definí BSPONMPI_LIBRARY y BSPONMPI_INCLUDE_DIR "
        "(o mirá docker/Dockerfile para compilarlo).")
endif()

message(STATUS "BSPonMPI: ${BSPONMPI_LIBRARY} (headers en ${BSPONMPI_INCLUDE_DIR})")

function(hpc_target name source)
    add_executable(${name} ${source})
    target_include_directories(${name} PRIVATE ${BSPONMPI_INCLUDE_DIR})
    target_link_libraries(${name} PRIVATE
        ${BSPONMPI_LIBRARY}
        MPI::MPI_CXX
        OpenMP::OpenMP_CXX
    )
    # OpenMPI 4.x incluye bindings C++ obsoletos que chocan con el
    # extern "C" { #include <bsp.h> } de BSPonMPI.
    target_compile_definitions(${name} PRIVATE OMPI_SKIP_MPICXX)
endfunction()

hpc_target(gas src/gas.cpp)
hpc_target(hello tests/hello.cpp)
hpc_target(halo_test tests/halo.cpp)
