
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @brief CMake build script for the linear-algebra-solvers library, which includes the concrete implementations of the
# LinearAlgebraMethods class based on the enabled libraries (HiCMA or Chameleon).
# @author Sameh Abdulah
# @author Mahmoud ElKarargy
# @date 2023-03-20

# Include the concrete implementations of the LinearAlgebraMethods class based on the enabled libraries (HiCMA or Chameleon)
if (EXAGEOSTAT_USE_CHAMELEON)
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/dense/ChameleonImplementationDense.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/diagonal-super-tile/ChameleonImplementationDST.cpp
            ${SOURCES}
            PARENT_SCOPE
            )
endif ()
if (EXAGEOSTAT_USE_HICMA)
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/tile-low-rank/HicmaImplementation.cpp
            ${SOURCES}
            PARENT_SCOPE
            )
endif ()