# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportBlasPP.cmake
# @brief This file searches for the BLAS++ library and includes it if not already included.
# @version 1.0.0
# @author Sameh Abdulah
# @author Mahmoud ElKarargy
# @date 2023-03-12

# search for BLAS library, if not already included
message("")
message("---------------------------------------- BLAS++")
message(STATUS "Checking for BLAS++")

if (NOT TARGET blaspp)

    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    include(ImportBlas)

    find_package(blaspp QUIET)

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
        set(url "https://github.com/icl-utk-edu/blaspp")
        set(tag "v2023.01.00")
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
