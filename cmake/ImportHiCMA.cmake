
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportHiCMA.cmake
# @brief Find and include HiCMA library as a dependency.
# @version 1.0.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-13

message("")
message("---------------------------------------- Hicma")
message(STATUS "Checking for HiCMA")
include(macros/BuildDependency)

if (NOT TARGET HICMA_FOUND)
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    find_package(HICMA QUIET)

    # If HiCMA is found, print its location.
    if (HICMA_FOUND)
        message("   Found HiCMA: ${HICMA_LIBDIR}")
        # If not found, install it.
    else ()
        message("   Can't find HiCMA, Installing it instead ..")

        # Set the flags to be passed to the build command.
        set(FLAGS -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/installdir/_deps/HICMA/ \-DHICMA_USE_MPI=${USE_MPI})
        set(ISCMAKE ON)
        set(ISGIT ON)
        set(AUTO_GEN OFF)
        # Build HiCMA from source.
        set(HICMA_DIR ${PROJECT_SOURCE_DIR}/installdir/_deps/HICMA/)
        BuildDependency(HiCMA "https://github.com/ecrc/hicma.git" "v1.0.0" ${FLAGS} ${ISCMAKE} ${ISGIT} ${AUTO_GEN})

        # Clear the flags.
        set(FLAGS "")

        # Find HiCMA after installation.
        find_package(HICMA REQUIRED)
    endif ()
else()
    message("   HiCMA already included")
endif()

# Include HiCMA headers in the project.
include_directories(${HICMA_INCLUDE_DIRS_DEP})
# TODO: Fix install control headers
include_directories(${PROJECT_SOURCE_DIR}/installdir/_deps/HICMA/hicma-src/hicma_ext/)

# Include HiCMA libraries in the project.
if (HICMA_LINKER_FLAGS)
    list(APPEND CMAKE_EXE_LINKER_FLAGS "${HICMA_LINKER_FLAGS}")
endif ()
if (HICMA_LIBRARY_DIRS)
    list(APPEND CMAKE_INSTALL_RPATH "${HICMA_LIBRARY_DIRS}")
    link_directories(${HICMA_LIBRARY_DIRS})
endif ()

# Add HiCMA libraries to the dependencies of the project.
if (HICMA_LIBRARIES_DEP)
    list(APPEND LIBS ${HICMA_LIBRARIES_DEP})
else ()
    list(APPEND LIBS ${HICMA_LIBRARIES})
endif()

message(STATUS "HiCMA done")