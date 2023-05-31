
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-04-06

if (EXAGEOSTAT_USE_CHAMELEON)
    set(EXAGEOSTAT_TESTFILES
            ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestChameleonImplementationDense.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestChameleonImplementationDST.cpp
            ${EXAGEOSTAT_TESTFILES}
            PARENT_SCOPE
            )
endif()

if (EXAGEOSTAT_USE_HiCMA)
    set(EXAGEOSTAT_TESTFILES
            ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestHiCMAImplementationTLR.cpp
            ${EXAGEOSTAT_TESTFILES}
            PARENT_SCOPE
            )
endif()