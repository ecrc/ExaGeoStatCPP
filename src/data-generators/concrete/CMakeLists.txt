
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief CMake configuration file for the Synthetic Generator module
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2025-02-18

# Add source files to the parent scope
if (RUNTIME_TYPE STREQUAL "STARPU")
        set(SOURCES
                ${CMAKE_CURRENT_SOURCE_DIR}/runtime/StarpuGenerator.cpp
                ${SOURCES}
                )
elseif (RUNTIME_TYPE STREQUAL "PARSEC")
        set(SOURCES
                ${CMAKE_CURRENT_SOURCE_DIR}/runtime/ParsecGenerator.cpp
                ${SOURCES}
                )
endif()

set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/SyntheticGenerator.cpp
        ${SOURCES}
        PARENT_SCOPE
        )