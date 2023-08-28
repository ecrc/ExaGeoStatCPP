
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @author Sameh Abdulah
# @author Mahmoud ElKarargy
# @date 2023-04-29
set(EXAGEOSTAT_TESTFILES

        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternStationary.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternNonStationary.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestBivariateMaternFlexible.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternDsigmaSquare.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestBivariateMaternParsimonious.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternNuggetsStationary.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateSpacetimeMaternStationary.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternDnu.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternDbeta.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternDdsigmaSquare.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternDdsigmaSquareBeta.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternDdsigmaSquareNu.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternDdbetaBeta.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternDdbetaNu.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternDdnuNu.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestBivariateSpacetimeMaternStationary.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternNonGaussian.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateExpNonGaussian.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestTrivariateMaternParsimonious.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestUnivariateMaternNonStat.cpp
        ${EXAGEOSTAT_TESTFILES}
        PARENT_SCOPE
        )
