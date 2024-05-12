
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-08-15

# Include the concrete implementations of the ExaGeoStat Descriptor class
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/concrete)

set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/ExaGeoStatDescriptor.cpp
        ${SOURCES}
        PARENT_SCOPE
        )