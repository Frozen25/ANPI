cmake_minimum_required(VERSION 3.5)

project(tarea04 VERSION 0.0.0 LANGUAGES CXX)

## Ensure that all binaries are stored in the same directory (usually build/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake") 

## Options
option(ANPI_ENABLE_SIMD "Force the use of optimized code instead of generic" on)
option(ANPI_ENABLE_OpenMP "Force the use of OpenMP" off)

## All compiler options
include(CompilerFlags)

## Create libraries and executables
add_subdirectory (src)


CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/AnpiConfig.hpp.in ${CMAKE_SOURCE_DIR}/include/AnpiConfig.hpp)


## Also build unit tests ...
enable_testing ()
add_subdirectory (test)

## ... and benchmarks
add_subdirectory(benchmarks)
