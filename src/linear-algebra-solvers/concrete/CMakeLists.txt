
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @brief CMake build script for the linear-algebra-solvers library, which includes the concrete implementations of the
# LinearAlgebraMethods class based on the enabled libraries (HiCMA or Chameleon).
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-20

# Include the concrete implementations of the LinearAlgebraMethods class based on the enabled libraries (HiCMA or Chameleon)
if (RUNTIME_TYPE STREQUAL "STARPU")
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/chameleon/dense/ChameleonDense.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/chameleon/dst/ChameleonDST.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/chameleon/ChameleonImplementation.cpp
            ${SOURCES}
    )

    if (USE_HICMA)
        list(APPEND SOURCES
                ${CMAKE_CURRENT_SOURCE_DIR}/tlr/HicmaImplementation.cpp
        )
    endif ()
endif ()
set(SOURCES ${SOURCES} PARENT_SCOPE)

