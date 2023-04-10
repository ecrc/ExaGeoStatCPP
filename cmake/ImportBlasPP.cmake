
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportBlasPP.cmake
# @brief This file searches for the BLAS++ library and includes it if not already included.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-12

# search for BLAS library, if not already included
message("")
message("---------------------------------------- BLAS++")
message(STATUS "Checking for BLAS++")

if (NOT TARGET blaspp)

    include(ImportBlas)

    # Find the BLAS++ library
    find_package(blaspp QUIET)

    message(${blaspp_FOUND})
    # If BLAS++ is found, include it
    if (blaspp_FOUND)
        message("Found BLAS++: ${blaspp_DIR}")
        # If BLAS++ is not found, but its CMakeLists.txt file exists, add it as a subdirectory and include it
    elseif (EXISTS "${CMAKE_SOURCE_DIR}/blaspp/CMakeLists.txt")
        set(build_tests_save "${build_tests}")
        set(build_tests "false")
        add_subdirectory("blaspp")

        set(build_tests "${build_tests_save}")
        set(blaspp_DIR "${CMAKE_BINARY_DIR}/blaspp")
        # If BLAS++ is not found and its CMakeLists.txt file does not exist, fetch it from the specified URL and include it
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

# Add BLAS++ to the list of libraries
set(LIBS
        blaspp
        ${LIBS}
        )

message(STATUS "BLAS++ done")