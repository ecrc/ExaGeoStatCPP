
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 2.0.0
# @brief CMake build script for the runtime-solvers library, which includes the RuntimeSolversMethods base class and the RuntimeSolversFactory
# @author Mahmoud ElKarargy
# @date 2024-11-04

# Include the concrete implementations of the RuntimeSolversMethods class
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/concrete)
# Define the sources for the library
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/RuntimeSolverFactory.cpp
        ${SOURCES}
        PARENT_SCOPE
        )
