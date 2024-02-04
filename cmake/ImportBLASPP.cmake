# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportBLASPP.cmake
# @brief This file searches for the BLAS++ library and includes it if not already included.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2024-02-04

include(ImportBLAS)

#Configurations
set(name blaspp)
string(TOUPPER ${name} capital_name)
set(tag "v2023.01.00")
set(url "https://github.com/icl-utk-edu/blaspp")
set(${name}_DIR "${CMAKE_INSTALL_PREFIX}/${capital_name}/${name}-build/")

message("")
message("---------------------------------------- ${capital_name}")
message(STATUS "Checking for ${capital_name} with Version ${version}")

# Check if the target is already included
IF (NOT TARGET ${name})
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    find_package(${name} ${version} QUIET COMPONENTS ${components})

    # If the package is found, print a message
    if (${name}_FOUND)
        message("   Found ${capital_name}; ${${name}_DIR} ${${name}_LIBRARIES}")
    else ()
        # If the package is not found, install it using BuildDependency
        message("   Can't find ${capital_name}, Installing it instead ..")
        include(FetchContent)
        set(FETCHCONTENT_BASE_DIR ${CMAKE_INSTALL_PREFIX}/${capital_name})
        FetchContent_Declare(${name} GIT_REPOSITORY "${url}" GIT_TAG "${tag}")
        FetchContent_MakeAvailable(${name})

    endif ()
else ()
    message(STATUS "${capital_name} already included")
endif ()

set(LIBS ${name} ${LIBS})
message(STATUS "${name} done")