
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportLapackPP.cmake
# @brief Find and include LAPACK++ library as a dependency.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-12

# search for LAPACK library, if not already included
message("")
message("---------------------------------------- LAPACK++")
message(STATUS "Checking for LAPACK++")
if (NOT TARGET lapackpp)

    # Include the ImportLapack script to ensure that LAPACK is included.
    include(ImportLapack)

    find_package(lapackpp QUIET)

    # If LAPACK++ is found, print its location.
    if (lapackpp_FOUND)
        message("   Found LAPACK++: ${lapackpp_DIR}")
        # If not found, install it.
    elseif (EXISTS "${CMAKE_SOURCE_DIR}/lapackpp/CMakeLists.txt")
        # Save the current value of build_tests and set it to false.
        set(build_tests_save "${build_tests}")
        set(build_tests "false")

        # Build LAPACK++ from source.
        add_subdirectory("lapackpp")

        # Restore the value of build_tests.
        set(build_tests "${build_tests_save}")

        # Set the location of LAPACK++.
        set(lapackpp_DIR "${CMAKE_BINARY_DIR}/lapackpp")
    else()
        # Save the current value of build_tests and set it to false.
        set(build_tests_save "${build_tests}")
        set(build_tests "false")

        # Fetch LAPACK++ source from the repository.
        set(url "https://github.com/icl-utk-edu/lapackpp")
        set(tag "master")
        message(STATUS "Fetching LAPACK++ ${tag} from ${url}")
        include(FetchContent)
        FetchContent_Declare(lapackpp GIT_REPOSITORY "${url}" GIT_TAG "${tag}")
        FetchContent_MakeAvailable(lapackpp)

        # Restore the value of build_tests.
        set(build_tests "${build_tests_save}")
    endif()
else()
    message("   LAPACK++ already included")
endif()

# Add LAPACK++ to the linking libraries.
set(LIBS
        lapackpp
        ${LIBS}
        )

# Add definition indicating version.
if ("${lapackpp_defines}" MATCHES "LAPACK_ILP64")
    set(COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS} -DHCORE_HAVE_LAPACK_WITH_ILP64")
endif()

message(STATUS "LAPACK++ done")