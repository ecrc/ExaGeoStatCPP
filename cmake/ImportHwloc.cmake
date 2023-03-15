
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportHwloc.cmake
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-15


message("")
message("---------------------------------------- Hwloc")
message(STATUS "Checking for Hwloc")

include(macros/BuildHwloc)

if (NOT TARGET HWLOC)
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)

    find_package(HWLOC 1.11.5 QUIET)

    if (HWLOC_FOUND)
        message("   Found HWLOC: ${HWLOC_LIBRARIES}")
    else ()
        set(hwloc_installpath ${CMAKE_BINARY_DIR}/_deps/hwloc-install)
        set(HWLOC_DIR "${hwloc_installpath}/")
        BuildHwloc(HWLOC "https://github.com/open-mpi/hwloc" "hwloc-2.4.0")
        find_package(HWLOC 1.11.5 REQUIRED)
    endif ()
else ()
    message("   HWLOC already included")
endif ()

list(APPEND LIBS  ${HWLOC_LIBRARIES})
link_directories(${HWLOC_LIBRARY_DIRS_DEP})
include_directories(${HWLOC_INCLUDE_DIRS})

message(STATUS "HWLOC done")
