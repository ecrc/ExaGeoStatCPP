
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportHwloc.cmake
# @brief Find and include Hwloc library as a dependency.
# @version 1.0.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-15


message("")
message("---------------------------------------- Hwloc")
message(STATUS "Checking for Hwloc")

include(macros/BuildDependency)

# If Hwloc library is not already included as a target, try to find it.
if (NOT TARGET HWLOC)
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    find_package(HWLOC 1.11.5 QUIET)

    # If Hwloc is found, print its location.
    if (HWLOC_FOUND)
        message("   Found HWLOC: ${HWLOC_LIBRARIES}")
        # If not found, install it.
    else ()
        message("   Can't find Hwloc, Installing it instead ..")

        # Set the flags to be passed to the build command.
        set(FLAGS --prefix=${PROJECT_SOURCE_DIR}/installdir/_deps/HWLOC/)
        set(ISCMAKE OFF)
        set(ISGIT ON)
        set(AUTO_GEN ON)

        # Build Hwloc from source.
        set(HWLOC_DIR  ${PROJECT_SOURCE_DIR}/installdir/_deps/HWLOC/)
        BuildDependency(HWLOC "https://github.com/open-mpi/hwloc" "hwloc-2.4.0" ${FLAGS} ${ISCMAKE} ${ISGIT} ${AUTO_GEN})

        # Clear the flags.
        set(FLAGS "")

        # Find Hwloc after installation.
        find_package(HWLOC 1.11.5 REQUIRED)
    endif ()
else ()
    message("   HWLOC already included")
endif ()

# Include Hwloc libraries in the project.
list(APPEND LIBS  ${HWLOC_LIBRARIES})
link_directories(${HWLOC_LIBRARY_DIRS_DEP})

# Include Hwloc headers in the project.
include_directories(${HWLOC_INCLUDE_DIRS})

message(STATUS "HWLOC done")