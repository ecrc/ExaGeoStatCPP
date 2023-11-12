
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportLapack.cmake
# @brief Find and include LAPACK library as a dependency.
# @version 1.0.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-12

# search for LAPACK library, if not already included
message("")
message("---------------------------------------- LAPACK")
message(STATUS "Checking for LAPACK")

include(macros/BuildDependency)

if (NOT TARGET LAPACK)
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    find_package(LAPACK QUIET)

    if (LAPACK_FOUND)
        message("   Found LAPACK: ${LAPACK_LIBRARIES}")
    else ()
        message("   Can't find Blas, Installing it instead ..")
        # Set installation flags
        set(FLAGS "")
        set(ISCMAKE ON)
        set(ISGIT ON)
        set(AUTO_GEN OFF)
        set(build_tests "false")
        BuildDependency(LAPACK "https://github.com/xianyi/OpenBLAS" "v0.3.21" ${FLAGS} ${ISCMAKE} ${ISGIT} ${AUTO_GEN})
        find_package(LAPACK REQUIRED)
    endif ()
else ()
    message("   LAPACK already included")
endif ()

message(STATUS "LAPACK done")
