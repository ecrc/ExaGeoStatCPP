
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief CMake configuration file for the Synthetic Generator module
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-02-14

# Add source files to the parent scope
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/SyntheticGenerator.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/CSVDataGenerator.cpp

        ${SOURCES}
        PARENT_SCOPE
        )