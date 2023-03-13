
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportBlasPP.cmake
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-12

# search for BLAS library, if not already included
message("")
message("---------------------------------------- BLAS++")
message(STATUS "Checking for BLAS++")

if (NOT TARGET blaspp)

    include(ImportBlas)

    find_package(blaspp QUIET)

    message(${blaspp_FOUND})
    if (blaspp_FOUND)
        message("Found BLAS++: ${blaspp_DIR}")
    elseif (EXISTS "${CMAKE_SOURCE_DIR}/blaspp/CMakeLists.txt")
        set(build_tests_save "${build_tests}")
        set(build_tests "false")
        add_subdirectory("blaspp")

        set(build_tests "${build_tests_save}")
        set(blaspp_DIR "${CMAKE_BINARY_DIR}/blaspp")
    else ()
        set(build_tests_save "${build_tests}")
        set(build_tests "false")
        set(url "https://bitbucket.org/icl/blaspp")
        set(tag "2021.04.01")
        message(STATUS "Fetching BLAS++ ${tag} from ${url}")
        include(FetchContent)
        FetchContent_Declare(
                blaspp GIT_REPOSITORY "${url}" GIT_TAG "${tag}")
        FetchContent_MakeAvailable(blaspp)
        set(build_tests "${build_tests_save}")
    endif ()
else ()
    message("   BLAS++ already included")
endif ()

set(LIBS
        blaspp
        ${LIBS}
        )
message(STATUS "BLAS++ done")
