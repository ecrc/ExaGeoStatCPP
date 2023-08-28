
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportChameleon.cmake
# @brief This script checks for Chameleon and includes it in the project if it is not already a target.
# @version 1.0.0
# @author Sameh Abdulah
# @author Mahmoud ElKarargy
# @date 2023-03-13

message("")
message("---------------------------------------- Catch2")
message(STATUS "Checking for Catch2")
include(macros/BuildDependency)

IF (NOT TARGET Catch2)
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    find_package(Catch2 QUIET)

    # If Catch2 is not found, fetch and build it
    if (Catch2_FOUND)
        message("   Found Catch2")
    else ()
        message("   Can't find catch2, Installing it instead ..")
        include(FetchContent)
        set(FETCHCONTENT_QUIET OFF)
        FetchContent_Declare(
                Catch2
                GIT_REPOSITORY https://github.com/catchorg/Catch2.git
                GIT_TAG v3.3.2 # Replace with the version of Catch2 you want to use for v3
                GIT_SHALLOW TRUE
        )
        FetchContent_MakeAvailable(Catch2)
    endif ()
else ()
    message(STATUS "Catch2 already included")
endif ()
