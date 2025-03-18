
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 2.0.0
# @brief CMake build script for the runtime-solvers library, which includes the concrete implementations of the
# RuntimeSolversMethods class based on the enabled runtime (StarPU or PaRSEC).
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2024-11-04

# Include the concrete implementations of the RuntimeSolversMethods class based on the enabled runtime (StarPU or PaRSEC)
if (RUNTIME_TYPE STREQUAL "STARPU")
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/StarpuRuntimeSolver.cpp
            ${SOURCES}
            )
elseif (RUNTIME_TYPE STREQUAL "PARSEC")
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/ParsecRuntimeSolver.cpp;
            ${SOURCES}
            )
endif ()
set(SOURCES ${SOURCES} PARENT_SCOPE)

