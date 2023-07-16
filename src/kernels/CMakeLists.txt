
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief This file contains the CMake configuration for the kernels directory.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-04-11

# Add the Configurations.cpp file to the list of source files.
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/Kernel.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternStationary.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternNonStationary.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/BivariateMaternFlexible.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/BivariateMaternParsimonious.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternNuggetsStationary.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateSpacetimeMaternStationary.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternDsigmaSquare.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternDnu.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternDbeta.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternDdsigmaSquare.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternDdsigmaSquareBeta.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternDdsigmaSquareNu.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternDdbetaBeta.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternDdbetaNu.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternDdnuNu.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/BivariateSpacetimeMaternStationary.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternNonGaussian.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateExpNonGaussian.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TrivariateMaternParsimonious.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/UnivariateMaternNonStat.cpp
        ${SOURCES}
        PARENT_SCOPE
        )