
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportSTARPU.cmake
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-13

message("")
message("---------------------------------------- StarPU")
message(STATUS "Checking for StarPU")

include(macros/BuildSTARPU)

if (NOT TARGET STARPU)
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)

    find_package(STARPU QUIET)

    if (STARPU_FOUND)
        message("   Found StarPU: ${STARPU_LIBRARIES}")
    else ()
        set(starpu_installpath ${CMAKE_BINARY_DIR}/_deps/starpu-install)
        set(STARPU_DIR "${starpu_installpath}/")
        BuildStarPU(STARPU "https://gitlab.inria.fr/starpu/starpu.git" "starpu-1.3.9")
        find_package(STARPU REQUIRED)
    endif ()
else ()
    message("   STARPU already included")
endif ()

list(APPEND LIBS  ${STARPU_LIBRARIES})
link_directories(${STARPU_LIBRARY_DIRS_DEP})
include_directories(${STARPU_INCLUDE_DIRS})

message(STATUS "starpu done")
