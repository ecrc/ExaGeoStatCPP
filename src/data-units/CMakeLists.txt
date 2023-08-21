
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-02-27

# Include the concrete implementations of the ExaGeoStat Descriptor class
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/descriptor)

set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/Locations.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/DescriptorData.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ExaGeoStatData.cpp
        ${SOURCES}
        PARENT_SCOPE
        )