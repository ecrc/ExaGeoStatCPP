
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

    find_package(STARPU 1.3.9 QUIET COMPONENTS ${STARPU_COMPONENT_LIST})

    if (STARPU_FOUND)
        message("   Found StarPU: ${STARPU_LIBRARIES}")
    else ()
        set(STARPU_DIR  ${PROJECT_SOURCE_DIR}/installdir/_deps/STARPU/)
        BuildStarPU(STARPU "https://gitlab.inria.fr/starpu/starpu.git" "starpu-1.3.9")
        find_package(STARPU 1.3.9 REQUIRED COMPONENTS ${STARPU_COMPONENT_LIST})
    endif ()
else ()
    message("   STARPU already included")
endif ()

list(APPEND LIBS  ${STARPU_LIBRARIES})
link_directories(${STARPU_LIBRARY_DIRS_DEP})
include_directories(${STARPU_INCLUDE_DIRS})
include_directories(${STARPU_INCLUDE_DIRS}/runtime/starpu)

include_directories(${STARPU_INCLUDE_DIRS_DEP})
if (STARPU_LINKER_FLAGS)
    list(APPEND CMAKE_EXE_LINKER_FLAGS "${STARPU_LINKER_FLAGS}")
endif ()
set(CMAKE_REQUIRED_INCLUDES "${STARPU_INCLUDE_DIRS_DEP}")
foreach (libdir ${STARPU_LIBRARY_DIRS_DEP})
    list(APPEND CMAKE_REQUIRED_FLAGS "-L${libdir}")
endforeach ()
set(CMAKE_REQUIRED_LIBRARIES "${STARPU_LIBRARIES_DEP}")

message(STATUS "starpu done")
