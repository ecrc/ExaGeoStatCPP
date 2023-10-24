# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportBlas.cmake
# @brief This file searches for the BLAS library and includes it if not already included.
# @version 1.0.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-12

# search for BLAS library, if not already included
message("")
message("---------------------------------------- BLAS")
message(STATUS "Checking for BLAS")
include(macros/BuildDependency)

if (NOT TARGET BLAS)
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    find_package(BLAS QUIET)

    if (BLAS_FOUND)
        message("   Found BLAS: ${BLAS_LIBRARIES}")
    else ()
        message("   Can't find Blas, Installing it instead ..")
        # Set installation flags
        set(FLAGS "")
        set(ISCMAKE ON)
        set(ISGIT ON)
        set(AUTO_GEN OFF)
        set(build_tests "false")
        BuildDependency(BLAS "https://github.com/xianyi/OpenBLAS" "v0.3.21" ${FLAGS} ${ISCMAKE} ${ISGIT} ${AUTO_GEN})
        find_package(BLAS REQUIRED)
    endif ()

else ()
    message("   BLAS already included")
endif ()

message(STATUS "BLAS done")
