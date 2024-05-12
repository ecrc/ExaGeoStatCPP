
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief This file contains the CMake configuration for the data-units directory.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-01-31

# Add the Configurations.cpp file to the list of source files.
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/Configurations.cpp
        ${SOURCES}
        PARENT_SCOPE
        )