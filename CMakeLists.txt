
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (C) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-01-30

#############################################################################
#                       THIS IS A TOP-LEVEL CMAKELISTS.txt                    #
#                                                                             #
#        It is intended to find all dependencies (required or optional)       #
#                    and set up corresponding variables                       #
###############################################################################

###############################################################################
##                            PRELIMINARIES                                  ##
###############################################################################

#Minimum CMAKE version.
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
cmake_policy(SET CMP0048 NEW)

# directly make an error if in-source build
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "In-source builds are not allowed.\n"
            "Please use ./config.sh and then ./clean_build.sh to build the code.")
endif ()

# Project Name and Version
project(exageostat-cpp VERSION 1.0.0)
# Set C++ standard to C++17.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
enable_language(CXX)

# Set Exageostat version.
set(EXAGEOSTAT_VERSION_MAJOR 1)
set(EXAGEOSTAT_VERSION_MINOR 0)
set(EXAGEOSTAT_VERSION_PATCH 0)
set(EXAGEOSTAT_VERSION ${EXAGEOSTAT_VERSION_MAJOR}.${EXAGEOSTAT_VERSION_MINOR}.${EXAGEOSTAT_VERSION_PATCH})
add_definitions("-DVERSION=${EXAGEOSTAT_VERSION}")

# CMAKE MODULES
# Add cmake files.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts/cmake)
# Include CMake to add ECRC cmake modules.
include(InitECRC)
# Include CMake to generate package config.
include(GenPkgConfig)

# Make sure pkg config is founded.
find_package(PkgConfig QUIET)
# Adding the headers directories needed for all projects.
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inst/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/prerequisites)

# Add src Directory to expose added libraries
add_subdirectory(src)

# Print installation path of Exageostat.
message(STATUS "Installation path : ${EXAGEOSTAT_INSTALL_PREFIX}")

# Enable/Disable Examples
option(EXAGEOSTAT_BUILD_EXAMPLES "Option to enable building example programs" OFF)
if (EXAGEOSTAT_BUILD_EXAMPLES)
    message(STATUS "Building Examples")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
endif ()

option(EXAGEOSTAT_BUILD_TESTS "Option to enable building tests" OFF)
if (EXAGEOSTAT_BUILD_TESTS)
    message(STATUS "Building Tests")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp-tests)
    include(CTest)
    add_test(NAME exageostat-tests
            COMMAND exageostat-tests
            )
endif ()

###############################################################################
##                            PRELIMINARIES                                  ##
###############################################################################


## Find R and Rcpp using FindR Module
#FIND_PACKAGE(R REQUIRED)
#if (${R_FOUND})
#    message(STATUS "Using R technology")
#    list(APPEND LIBS "R")
#    set(USE_R ON)
#    add_definitions(-DUSING_R)
#endif ()

