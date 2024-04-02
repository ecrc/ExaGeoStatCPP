# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportDependency.cmake
# @brief CMake script for importing and building external dependencies.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Amr Nasr
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

# Define a macro named ImportDependency for handling external dependencies. The macro checks for the dependency's presence and installs it if missing.
macro(ImportDependency name tag version url flag components is_cmake is_git auto_gen)

    # Check if the installation prefix is set to a system path (like /usr/) and warn the user about potential need for administrative privileges.
    if (CMAKE_INSTALL_PREFIX MATCHES "/usr/")
        message(WARNING "Installation path not specified. Please set the installation path using -DCMAKE_INSTALL_PREFIX=path/to/install or execute ./config.sh. Otherwise, please note that administrative privileges may be required to install in system paths.")
    endif ()

    # Convert the dependency name to uppercase for consistent messaging.
    string(TOUPPER ${name} capital_name)
    # Begin a section in the output to visually separate the handling of this dependency.
    message("")
    message("---------------------------------------- ${capital_name}")
    # Log the attempt to check for the specified version of the dependency.
    message(STATUS "Checking for ${capital_name} with Version ${version}")

    # Include the previously defined BuildDependency macro script for potential use.
    include(macros/BuildDependency)

    # If the dependency has not already been targeted for building in the current CMake process, proceed to check its presence.
    IF (NOT TARGET ${name})
        # Use the FindPkgConfig module to potentially use pkg-config for finding installed libraries.
        include(FindPkgConfig)
        find_package(PkgConfig QUIET)
        # Attempt to find the specified version of the package quietly, without generating much output.
        find_package(${name} ${version} QUIET COMPONENTS ${components})

        # If the package is found, output a message detailing the found configuration.
        if (${name}_FOUND)
            message("   Found ${capital_name}; ${${name}_DIR} ${${name}_LIBRARIES}")
        else ()
            # If the package is not found, notify and invoke BuildDependency to install it.
            message("   Can't find ${capital_name}, Installing it instead ..")
            BuildDependency(${name} ${url} ${tag} "${flag}" ${is_cmake} ${is_git} ${auto_gen})
            # After attempting installation, forcibly attempt to find the package again, this time requiring its presence.
            find_package(${name} ${version} REQUIRED COMPONENTS ${components})
        endif ()
    else ()
        # If the dependency target already exists, log that it's already been included.
        message(STATUS "${capital_name} already included")
    endif ()

    # Setup link and include directories based on the found or installed package configuration.
    # Add the dependency's library directories to the link directories for the current CMake target.
    link_directories(${${name}_LIBRARY_DIRS_DEP})
    link_directories(${${name}_LIBRARY_DIRS})
    # Add the dependency's include directories to the include path.
    include_directories(${${name}_INCLUDE_DIRS})
    include_directories(AFTER ${${name}_INCLUDE_DIRS_DEP})

    # If the dependency is not GSL, append its libraries to the list of libraries to be linked against.
    if(NOT ${name} STREQUAL "GSL")
        list(APPEND LIBS ${${name}_LIBRARIES})
    endif()
    # Append any additional dependency libraries to the list of libraries.
    list(APPEND LIBS ${${name}_LIBRARIES_DEP})

endmacro()
