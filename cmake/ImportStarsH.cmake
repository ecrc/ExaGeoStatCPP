
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief Find and include STARSH library as a dependency.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-13

# Configuration parameters for integrating the STARSH library
# 'name' is set to "STARSH" to identify the STARSH library within this script.
set(name "STARSH")
# 'tag' specifies "v0.3.1" as the version tag for STARSH, denoting the exact release to be used.
set(tag "v0.3.1")
# 'version' sets "0.3.1" as the version of the STARSH library, ensuring it aligns with project requirements.
set(version "0.3.1")
# 'flag' is used for additional build configuration options, specifically disabling StarPU and optionally enabling MPI.
set(flag \-DSTARPU=OFF \-DMPI=${USE_MPI})
# 'is_cmake' indicates that STARSH uses CMake as its build system, set to ON.
set(is_cmake ON)
# 'is_git' denotes that the source code for STARSH is hosted on a Git repository, set to ON.
set(is_git ON)
# 'auto_gen' signals whether autogen scripts are needed for the build process; it is set to OFF for STARSH.
set(auto_gen OFF)
# 'url' provides the GitHub repository URL for STARSH, specifying the source code's location.
set(url "https://github.com/ecrc/stars-h.git")

# The 'ImportDependency' macro, located in the 'macros' directory, is included to manage the import and setup of the STARSH library.
include(macros/ImportDependency)
# The 'ImportDependency' macro is called with the configuration parameters set above to manage the detection, fetching, and setup of STARSH.
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

# A message is output to indicate the successful integration of the STARSH library into the project.
message(STATUS "${name} done")
