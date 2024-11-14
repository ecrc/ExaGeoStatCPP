
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @brief CMake build script for StarPu functions
# @author Mahmoud ElKarargy
# @date 2024-03-10

# Include runtime directory,based on runtime flag.
if ("${RUNTIME_TYPE}" STREQUAL "PARSEC")
    add_subdirectory(parsec)
    set(JDF_GENERATED_SOURCES
            ${JDF_GENERATED_SOURCES}
            PARENT_SCOPE
            )
else ()
    #by default use StarPu runtime.
    add_subdirectory(starpu)
endif ()

# Define the sources for the library
set(SOURCES
        ${SOURCES}
        PARENT_SCOPE
        )
