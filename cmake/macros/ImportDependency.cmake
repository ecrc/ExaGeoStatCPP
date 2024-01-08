# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportDependency.cmake
# @brief CMake script for importing and building external dependencies.
# @version 1.0.1
# @author Mahmoud ElKarargy
# @date 2023-12-28

# ImportDependency Macro:
# This macro is designed to check for the presence of an external dependency and install it if not found.
# It takes the following parameters:
# - raw_name: The name of the dependency.
# - tag: The tag of the dependency to fetch.
# - version: The version of the dependency.
# - url: The URL of the repository or source tarball.
# - flag: Additional flags to pass to the configure/make commands.
# - is_cmake: A boolean flag indicating whether the dependency uses CMake as its build system.
# - is_git: A boolean flag indicating whether the dependency is hosted on a git repository.
# - auto_gen: A boolean flag indicating whether to use autogen scripts or not.

# The macro checks whether the dependency is already included. If not, it attempts to find the package.
# If the package is found, it prints a message. If not, it calls the BuildDependency macro to fetch,
# configure, build, and install the dependency. Finally, it attempts to find the package again to validate the installation.

macro(ImportDependency name tag version url flag is_cmake is_git auto_gen)

    # First, Check if no path is set for installation.
    if (CMAKE_INSTALL_PREFIX MATCHES "/usr/")
        message(WARNING "Installation path not specified. Please set the installation path using -DCMAKE_INSTALL_PREFIX=path/to/install or execute ./config.sh. Otherwise, please note that administrative privileges may be required to install in system paths.")
    endif ()

    # Convert name to uppercase for consistency
    string(TOUPPER ${name} capital_name)
    message("")
    message("---------------------------------------- ${capital_name}")
    message(STATUS "Checking for ${capital_name} with Version ${version}")

    # Include the BuildDependency macro
    include(macros/BuildDependency)

    # Check if the target is already included
    IF (NOT TARGET ${name})
        include(FindPkgConfig)
        find_package(PkgConfig QUIET)
        find_package(${name} ${version} QUIET)

        # If the package is found, print a message
        if (${name}_FOUND)
            message("   Found ${capital_name}; ${${name}_DIR} ${${name}_LIBRARIES}")
        else ()
            # If the package is not found, install it using BuildDependency
            message("   Can't find ${capital_name}, Installing it instead ..")
            BuildDependency(${name} ${url} ${tag} "${flag}" ${is_cmake} ${is_git} ${auto_gen})
            find_package(${name} ${version} REQUIRED)
        endif ()
    else ()
        message(STATUS "${capital_name} already included")
    endif ()

    # Include and link to the installation directory of the dependency
    link_directories(${${name}_LIBRARY_DIRS_DEP})
    link_directories(${${name}_LIBRARY_DIRS})
    include_directories(${${name}_INCLUDE_DIRS})
    include_directories(AFTER ${${name}_INCLUDE_DIRS_DEP})
    list(APPEND LIBS ${${name}_LIBRARIES})
    list(APPEND LIBS ${${name}_LIBRARIES_DEP})

endmacro()