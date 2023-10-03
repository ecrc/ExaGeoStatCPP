
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @author Mahmoud ElKarargy
# @date 2023-04-06

set(EXAGEOSTAT_TESTFILES
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestChameleonImplementationDense.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestChameleonImplementationDST.cpp
        ${EXAGEOSTAT_TESTFILES}
)

if (EXAGEOSTAT_USE_HICMA)
    list(APPEND EXAGEOSTAT_TESTFILES
            ${CMAKE_CURRENT_SOURCE_DIR}/concrete/TestHiCMAImplementationTLR.cpp
    )
endif ()

set(EXAGEOSTAT_TESTFILES ${EXAGEOSTAT_TESTFILES} PARENT_SCOPE)