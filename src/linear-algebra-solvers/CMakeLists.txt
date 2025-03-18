
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @brief CMake build script for the linear-algebra-solvers library, which includes the LinearAlgebraMethods base class and the LinearAlgebraFactory class for creating linear algebra solvers for different computations using HiCMA or Chameleon libraries.
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-20

# Include the concrete implementations of the LinearAlgebraMethods class
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/concrete)

# Define the sources for the library
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/LinearAlgebraFactory.cpp
        ${SOURCES}
        )

if (RUNTIME_TYPE STREQUAL "STARPU")
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/LinearAlgebraMethods.cpp
            ${SOURCES}
            )
endif ()

set(SOURCES ${SOURCES} PARENT_SCOPE)
