
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief CMake configuration file for the data loader module
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2023-02-14

# Add source files to the parent scope
if (RUNTIME_TYPE STREQUAL "STARPU")
        set(SOURCES
                ${CMAKE_CURRENT_SOURCE_DIR}/CSVLoader.cpp
                ${SOURCES}
                PARENT_SCOPE
        )
elseif (RUNTIME_TYPE STREQUAL "PARSEC")
        set(SOURCES
                ${CMAKE_CURRENT_SOURCE_DIR}/ParsecLoader.cpp
                ${SOURCES}
                PARENT_SCOPE
        )
endif()
