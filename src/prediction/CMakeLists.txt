
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-06-08

set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/Prediction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/PredictionHelpers.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/PredictionAuxiliaryFunctions.cpp
        ${SOURCES}
        PARENT_SCOPE
        )