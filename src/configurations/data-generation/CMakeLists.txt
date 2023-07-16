
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief This file contains the CMake configuration for the configurations directory.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-02-03

# Add the source files to the list of sources.
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/SyntheticDataConfigurations.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/DataConfigurations.cpp
        ${SOURCES}
        PARENT_SCOPE
        )