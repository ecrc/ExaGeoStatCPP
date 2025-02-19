
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportHwloc.cmake
# @brief Find and include Hwloc library as a dependency.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-15

# Configuration settings for integrating the HWLOC library
# 'name' sets the identifier for the HWLOC library within this script to "HWLOC".
set(name "Hwloc")
# 'tag' specifies "hwloc-2.10.0" as the version tag, identifying a specific release of HWLOC to be used.
set(tag "hwloc-2.10.0")
# 'version' defines "2.10.0" as the version of HWLOC, ensuring it meets project compatibility requirements.
set(version "2.10.0")
# 'flag' is available for additional build configuration options but remains empty for HWLOC.
set(flag "")
# 'is_cmake' indicates whether HWLOC uses CMake for its build system. It's set to OFF, suggesting an alternative build system is used.
set(is_cmake OFF)
# 'is_git' denotes that HWLOC's source code is maintained in a Git repository, set to ON.
set(is_git ON)
# 'auto_gen' signals the need for running autogen scripts as part of the build process for HWLOC, set to ON.
set(auto_gen ON)
# 'url' provides the GitHub repository URL for HWLOC, indicating the source code's location.
set(url "https://github.com/open-mpi/hwloc")

# Includes the 'ImportDependency' macro script, located in the 'macros' directory, responsible for handling the import and setup of dependencies.
include(macros/ImportDependency)
# The 'ImportDependency' macro is called with the configuration parameters above to manage the detection, fetching, and setup of HWLOC.
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

# A message is logged to indicate the successful integration of the HWLOC library into the project.
message(STATUS "${name} done")

