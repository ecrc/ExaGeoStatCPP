
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportLapack.cmake
# @brief Find and include LAPACK library as a dependency.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-12

# Configuration settings for the integration of the LAPACK library
# 'name' is designated as "LAPACK" to identify the LAPACK library within the scope of this script.
set(name "LAPACK")
# 'tag' is set to "v0.3.21", specifying the particular version tag of LAPACK to be utilized.
set(tag "v0.3.21")
# 'version' denotes the LAPACK library version as "0.3.21", aligned with the tag for consistency in versioning.
set(version "0.3.21")
# 'flag' is intended for any additional flags needed for configuration or building, but is left blank in this case.
set(flag "")
# 'is_cmake' is a boolean flag indicating that LAPACK uses CMake for its build process, set to ON.
set(is_cmake ON)
# 'is_git' signifies that the source code for LAPACK is available in a Git repository, set to ON.
set(is_git ON)
# 'auto_gen' indicates whether autogen scripts are necessary for the configuration process; it is set to OFF for LAPACK.
set(auto_gen OFF)
# 'url' provides the repository URL for LAPACK, pointing to the location where the source code can be accessed.
set(url "https://github.com/xianyi/OpenBLAS")

# The 'ImportDependency' macro script from the 'macros' directory is included, which facilitates the import and setup of dependencies.
include(macros/ImportDependency)
# The 'ImportDependency' macro is executed with the above-defined parameters to manage the detection, retrieval, and configuration of LAPACK.
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

# A status message is output to indicate the successful completion of the LAPACK setup process.
message(STATUS "${name} done")

