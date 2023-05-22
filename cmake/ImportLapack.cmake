
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
    # Try to find LAPACK.
    find_package(LAPACK REQUIRED)

    # If LAPACK is found, print its location.
    if (LAPACK_FOUND)
        message("   Found LAPACK: ${LAPACK_LIBRARIES}")
        # If not found, install it.
    else ()
        # Save the current value of build_tests and set it to false.
        set(build_tests_save "${build_tests}")
        set(build_tests "false")

        # Build OpenBLAS from source.
        BuildDependency(blas "https://github.com/xianyi/OpenBLAS" "v0.3.21")

        # Restore the value of build_tests.
        set(build_tests "${build_tests_save}")

        # Find LAPACK after installation.
        find_package(LAPACK REQUIRED)
    endif ()
else ()
    message("   LAPACK already included")
endif ()

# Include LAPACK libraries in the project.
list(APPEND LIBS  ${LAPACK_LIBRARIES})
link_directories(${LAPACK_LIBRARY_DIRS})

# Include LAPACK headers in the project.
include_directories(${LAPACK_INCLUDE_DIRS})

message(STATUS "LAPACK done")