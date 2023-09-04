
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief CMake configuration file for Data Generators module
# @version 1.0.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-02-14

# Add subdirectories for concrete implementations
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/concrete)

# Add source files to the parent scope
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/DataGenerator.cpp
        ${SOURCES}
        PARENT_SCOPE
        )
