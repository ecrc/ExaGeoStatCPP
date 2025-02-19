
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief CMake configuration file for data-transformer directory
# @version 2.0.0
# @author Mahmoud ElKarargy
# @date 2024-10-15

# Add source files to the parent scope
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/DataTransformer.cpp
        ${SOURCES}
        PARENT_SCOPE
)
