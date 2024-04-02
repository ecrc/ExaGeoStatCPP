
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @brief CMake build script for StarPu functions
# @author Mahmoud ElKarargy
# @date 2024-02-19

# Include the concrete implementations of the StarPu helpers classes and StarPu Helpes Factory
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/StarPuHelpersFactory.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/concrete/ChameleonStarPuHelpers.cpp
        ${SOURCES}
        )

if (USE_HICMA)
    list(APPEND SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/concrete/HicmaStarPuHelpers.cpp)
endif ()

set(SOURCES ${SOURCES} PARENT_SCOPE)
