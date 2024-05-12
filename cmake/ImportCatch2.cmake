
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportChameleon.cmake
# @brief This script checks for Chameleon and includes it in the project if it is not already a target.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2023-03-13

# Sets the 'name' variable to "Catch2", identifying the dependency being imported.
set(name "Catch2")
# Specifies the 'tag' variable as "v3.3.2", which likely corresponds to a specific release or version of Catch2.
set(tag "v3.3.2")
# Sets the 'version' variable to "3.3.2", defining the expected version of Catch2 to be used in the project.
set(version "3.3.2")
# Initializes the 'flag' variable as an empty string, which could be used for additional configuration options during the build or installation process.
set(flag "")
# Sets the 'is_cmake' flag to ON, indicating that Catch2 uses CMake as its build system.
set(is_cmake ON)
# Sets the 'is_git' flag to ON, suggesting that Catch2's source code is hosted in a git repository.
set(is_git ON)
# Initializes the 'auto_gen' flag to OFF, implying that there's no need for auto-generation scripts (like autogen.sh) in the building process of Catch2.
set(auto_gen OFF)
# Defines the 'url' variable with the GitHub repository URL of Catch2, specifying the source location from which Catch2 will be fetched.
set(url "https://github.com/catchorg/Catch2.git")

# Includes the 'ImportDependency' macro script located in the 'macros' directory. This script is responsible for checking the presence of the specified dependency and importing it if necessary.
include(macros/ImportDependency)
# Calls the 'ImportDependency' macro with the previously set variables as arguments. This macro will check for Catch2's presence, and if it's not found, it will attempt to fetch and install it using the provided details.
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})
# Logs a message indicating the completion of the Catch2 import process.
message(STATUS "${name} done")
