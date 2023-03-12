
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportBlas.cmake
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-12

# search for BLAS library, if not already included
message("")
message("---------------------------------------- BLAS")
message(STATUS "Checking for BLAS")

include(macros/BuildDependency)

if (NOT TARGET BLAS)

    find_package(BLAS QUIET)

    if (BLAS_FOUND)
        message("   Found BLAS: ${BLAS_LIBRARIES}")
    else ()
        set(build_tests_save "${build_tests}")
        set(build_tests "false")
        BuildDependency(blas "https://github.com/xianyi/OpenBLAS" "v0.3.21")
        set(build_tests "${build_tests_save}")
    endif ()

else ()
    message("   BLAS already included")
endif ()

message(STATUS "BLAS done")
