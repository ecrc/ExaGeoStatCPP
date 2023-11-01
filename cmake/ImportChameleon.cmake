
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportChameleon.cmake
# @brief This script checks for Chameleon and includes it in the project if it is not already a target.
# @version 1.0.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-13

message("")
message("---------------------------------------- Chameleon")
message(STATUS "Checking for Chameleon")
include(macros/BuildDependency)

if (NOT TARGET CHAMELEON_FOUND)
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    find_package(CHAMELEON QUIET)

    # If Chameleon is found, include it
    if (CHAMELEON_FOUND)
        message("   Found Chameleon: ${CHAMELEON_LIBDIR}")
    # If Chameleon is not found, install it
    else ()
        message("   Can't find Chameleon, Installing it instead ..")
        # Set installation flags
        set(FLAGS -DCHAMELEON_USE_CUDA=${USE_CUDA} \-DBLA_VENDOR=Intel10_64lp \-DCHAMELEON_USE_MPI=${USE_MPI} \-DCHAMELEON_SCHED_STARPU=ON \-DCHAMELEON_ENABLE_TESTING=OFF)
        set(ISCMAKE ON)
        set(ISGIT ON)
        set(AUTO_GEN OFF)
        # Install Chameleon
        BuildDependency(CHAMELEON "https://gitlab.inria.fr/solverstack/chameleon.git" "v1.1.0" ${FLAGS} ${ISCMAKE} ${ISGIT} ${AUTO_GEN})
        # Reset flags
        set(FLAGS "")
        find_package(CHAMELEON REQUIRED)
    endif ()
else()
    message("   CHAMELEON already included")
endif()

link_directories(${CHAMELEON_LIBRARY_DIRS_DEP})

include_directories(AFTER ${CHAMELEON_INCLUDE_DIRS_DEP})
include_directories(AFTER ${CHAMELEON_DIR_FOUND}/include/coreblas)
include_directories(${CHAMELEON_DIR_FOUND}/chameleon-src)

if (CHAMELEON_LINKER_FLAGS)
    list(APPEND CMAKE_EXE_LINKER_FLAGS "${CHAMELEON_LINKER_FLAGS} ")
endif ()
if (CHAMELEON_LIBRARY_DIRS)
    # the RPATH to be used when installing
    list(APPEND CMAKE_INSTALL_RPATH "${CHAMELEON_LIBRARY_DIRS}")
endif ()
if (CHAMELEON_LIBRARIES)
    if (CHAMELEON_LIBRARIES_DEP)
        list(APPEND LIBS ${CHAMELEON_LIBRARIES_DEP})
    else ()
        list(APPEND LIBS ${CHAMELEON_LIBRARIES})
    endif ()
endif ()

message(STATUS "Chameleon done")
