
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @brief CMake build script for the operators library, which includes the OperatorMethods base class and the OperatorFactory class for creating operators of different types.
# @author Mahmoud ElKarargy
# @date 2023-05-30

# Define the sources for the library
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/ExaGeoStat.cpp
        ${SOURCES}
        PARENT_SCOPE
        )