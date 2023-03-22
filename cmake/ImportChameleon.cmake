
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportChameleon.cmake
# @version 1.0.0
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

    if (CHAMELEON_FOUND)
        message("   Found Chameleon: ${CHAMELEON_LIBDIR}")
    else ()
        message("   Can't find Chameleon, Installing it instead ..")
        set(FLAGS -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/installdir/_deps/CHAMELEON/ " " -DCHAMELEON_USE_CUDA=${USE_CUDA} " " -DCHAMELEON_USE_MPI=${USE_MPI} " " -DCHAMELEON_SCHED_STARPU=ON)
        set(ISCMAKE ON)
        set(ISGIT ON)
        set(CHAMELEON_DIR  ${PROJECT_SOURCE_DIR}/installdir/_deps/CHAMELEON/)
        BuildDependency(CHAMELEON "https://gitlab.inria.fr/solverstack/chameleon.git" "v1.1.0" ${FLAGS} ${ISCMAKE} ${ISGIT})
        find_package(CHAMELEON REQUIRED)
    endif ()
else()
    message("   CHAMELEON already included")
endif()

link_directories(${CHAMELEON_LIBRARY_DIRS_DEP})

include_directories(AFTER ${CHAMELEON_INCLUDE_DIRS_DEP})
include_directories(AFTER ${CHAMELEON_DIR_FOUND}/include/coreblas)
if (CHAMELEON_LINKER_FLAGS)
    list(APPEND CMAKE_EXE_LINKER_FLAGS "${CHAMELEON_LINKER_FLAGS}")
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