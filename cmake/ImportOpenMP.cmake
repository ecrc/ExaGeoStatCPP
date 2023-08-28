
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportOpenMP.cmake
# @brief Find and include OpenMP library as a dependency.
# @version 1.0.0
# @author Sameh Abdulah
# @author Mahmoud ElKarargy
# @date 2023-03-13

# Add OpenMP if requested.
option(USE_OPENMP "Use OpenMP, if available" true)
if (NOT USE_OPENMP)
    message(STATUS "User has requested to NOT use OpenMP")
else ()
    find_package(OpenMP QUIET)
    IF (OPENMP_FOUND)
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(LIBS
                OpenMP::OpenMP_CXX
                ${LIBS}
                )
    ENDIF ()
endif ()
