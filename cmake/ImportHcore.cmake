
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportHcore.cmake
# @brief Checks for the Hcore library and includes it in the project if it is not already present.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-15

message("")
message("---------------------------------------- Hcore")
message(STATUS "Checking for Hcore")

include(macros/BuildDependency)

# Check if the Hcore library is already included in the project.
if (NOT TARGET HCORE_FOUND)

    # If not, attempt to find it with PkgConfig and CMake's find_package function.
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    find_package(HCORE QUIET)

    # If the Hcore library is found, add it to the project's libraries.
    if (HCORE_FOUND)
        message("   Found Hcore: ${HCORE_LIBDIR}")
    else()

        # If the Hcore library is not found, install it and add it to the project's libraries.
        message("   Can't find Hcore, Installing it instead ..")
        set(FLAGS -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/installdir/_deps/HCORE)
        set(ISCMAKE ON)
        set(ISGIT ON)
        BuildDependency(HCORE "https://github.com/ecrc/hcore.git" "v0.1.3" ${FLAGS} ${ISCMAKE} ${ISGIT})
        set(FLAGS "")
        find_package(HCORE REQUIRED)
    endif()
else()
    message("   HCORE already included")
endif()

# Add the Hcore library to the project's list of libraries.
list(APPEND LIBS ${HCORE_LIBRARIES})
link_directories(${HCORE_LIBRARY_DIRS_DEP})
include_directories(${HCORE_INCLUDE_DIRS})

message(STATUS "Hcore done")