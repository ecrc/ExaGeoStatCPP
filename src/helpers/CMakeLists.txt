
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-06-08

set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/DistanceCalculationHelpers.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/CommunicatorMPI.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ByteHandler.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/BasselFunction.cpp
        ${SOURCES}
        PARENT_SCOPE
        )