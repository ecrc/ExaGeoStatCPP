
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportDoxygen.cmake
# @brief This script checks for Doxygen and installs it if it is not already installed using the ImportDependency method.
# @version 1.0.1
# @author Mahmoud ElKarargy
# @date 2023-01-22


# Set configurations for Doxygen
set(name "Doxygen")
set(tag "Release_1_9_1")  # Use the tag corresponding to the desired version
set(version "1.9.1")
set(flag -G "Unix Makefiles")
set(url "https://github.com/doxygen/doxygen.git")
set(is_cmake ON)  # Doxygen uses CMake
set(is_git ON)    # We will be cloning from a git repository
set(auto_gen OFF) # Auto generation of project files is not needed

if(EXISTS "${CMAKE_INSTALL_PREFIX}/DOXYGEN/bin/doxygen")
    message("here")
    set(DOXYGEN_EXECUTABLE "${CMAKE_INSTALL_PREFIX}/DOXYGEN/bin/doxygen")
endif()

# Include the ImportDependency script
include(macros/ImportDependency)
# Call ImportDependency to handle the fetching and building of Doxygen
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

message(STATUS "${name} done")
