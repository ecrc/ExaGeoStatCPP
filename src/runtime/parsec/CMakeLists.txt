
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @brief CMake build script for Parsec runtime
# @author Mahmoud ElKarargy
# @date 2024-03-10

# Define the sources for the library
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/ParsecFunctions.cpp
        ${SOURCES}
        PARENT_SCOPE
        )