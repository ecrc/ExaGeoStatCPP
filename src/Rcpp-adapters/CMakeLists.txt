
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2024-01-29

set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/FunctionsAdapter.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/RcppExports.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/RcppModules.cpp
        ${SOURCES}
        PARENT_SCOPE
        )