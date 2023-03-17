
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

include(macros/BuildDependency)

if (NOT TARGET HWLOC)
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    find_package(HWLOC 1.11.5 QUIET)

    if (HWLOC_FOUND)
        message("   Found HWLOC: ${HWLOC_LIBRARIES}")
    else ()
        message("   Can't find Hwloc, Installing it instead ..")
        set(FLAGS --prefix=${PROJECT_SOURCE_DIR}/installdir/_deps)
        set(ISCMAKE OFF)
        set(ISGIT ON)
        BuildDependency(HWLOC "https://github.com/open-mpi/hwloc" "hwloc-2.4.0" ${FLAGS} ${ISCMAKE} ${ISGIT})
        find_package(HWLOC 1.11.5 REQUIRED)
    endif ()
else ()
    message("   HWLOC already included")
endif ()

list(APPEND LIBS  ${HWLOC_LIBRARIES})
link_directories(${HWLOC_LIBRARY_DIRS_DEP})
include_directories(${HWLOC_INCLUDE_DIRS})

message(STATUS "HWLOC done")
