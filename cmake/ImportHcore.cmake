
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
include(macros/BuildHcore)

if (NOT TARGET HCORE_FOUND)
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)

    find_package(HCORE QUIET)

    if (HCORE_FOUND)
        message("   Found Chameleon: ${CHAMELEON_LIBDIR}")
    else ()
        set(hcore_installpath ${CMAKE_BINARY_DIR}/_deps/hcore-install)
        set(HCORE_DIR "${chameleon_installpath}/")
        BuildHcore(HCORE "https://github.com/ecrc/hcore.git" "v0.1.3")
        find_package(HCORE REQUIRED)
    endif ()
else()
    message("   HCORE already included")
endif()

message(STATUS "Hcore done")