
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief This file contains the CMake configuration for the kernels directory.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-04-11

# Add the Configurations.cpp file to the list of source files.
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/UnivariateMaternStationary.cpp
        ${SOURCES}
        PARENT_SCOPE
        )