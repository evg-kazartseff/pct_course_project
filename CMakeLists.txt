cmake_minimum_required(VERSION 3.4)
project(mpi_jacobi_method)

set(CMAKE_C_STANDARD 11)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/bin)

add_subdirectory(serial/src)
add_subdirectory(parallel/src)
