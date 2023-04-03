
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportHcore.cmake
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-15

message("")
message("---------------------------------------- Hcore")
message(STATUS "Checking for Hcore")
include(macros/BuildDependency)

if (NOT TARGET HCORE_FOUND)
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    find_package(HCORE QUIET)

    if (HCORE_FOUND)
        message("   Found Hcore: ${HCORE_LIBDIR}")
    else ()
        message("   Can't find Hcore, Installing it instead ..")
        set(FLAGS -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/installdir/_deps/HCORE)
        set(ISCMAKE ON)
        set(ISGIT ON)
        BuildDependency(HCORE "https://github.com/ecrc/hcore.git" "v0.1.3" ${FLAGS} ${ISCMAKE} ${ISGIT})
        set(FLAGS "")
        find_package(HCORE REQUIRED)
    endif ()
else()
    message("   HCORE already included")
endif()

message(STATUS "Hcore done")

list(APPEND LIBS  ${HCORE_LIBRARIES})
link_directories(${HCORE_LIBRARY_DIRS_DEP})
include_directories(${HCORE_INCLUDE_DIRS})