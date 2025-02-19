
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @brief CMake build script for the descriptors library, which includes the concrete implementations of the
# Descriptors class based on the enabled libraries (HiCMA or Chameleon).
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-20

# Include the concrete implementations of the Descriptor class based on the enabled libraries (HiCMA or Chameleon)
if (RUNTIME_TYPE STREQUAL "STARPU")
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/ChameleonDescriptor.cpp
            ${SOURCES}
            )
    if (USE_HICMA)
        list(APPEND SOURCES
                ${CMAKE_CURRENT_SOURCE_DIR}/HicmaDescriptor.cpp
                )
    endif ()
else ()
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/ParsecDescriptor.cpp
            ${SOURCES}
            )
endif ()

set(SOURCES ${SOURCES} PARENT_SCOPE)
