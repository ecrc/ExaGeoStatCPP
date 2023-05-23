
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportLapack.cmake
# @brief Find and include LAPACK library as a dependency.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-12

# search for LAPACK library, if not already included
message("")
message("---------------------------------------- LAPACK")
message(STATUS "Checking for LAPACK")

include(macros/BuildD)

if (NOT TARGET LAPACK)
    find_package(LAPACK REQUIRED)
    if (LAPACK_FOUND)
        message("   Found LAPACK: ${LAPACK_LIBRARIES}")
    else ()
        set(build_tests_save "${build_tests}")
        set(build_tests "false")
        BuildDependency(blas "https://github.com/xianyi/OpenBLAS" "v0.3.21")
        set(build_tests "${build_tests_save}")
        find_package(LAPACK REQUIRED)
    endif ()
else ()
    message("   LAPACK already included")
endif ()

message(STATUS "LAPACK done")
