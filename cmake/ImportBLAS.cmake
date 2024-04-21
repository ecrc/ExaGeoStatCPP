# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportBLAS.cmake
# @brief This file searches for the BLAS library and includes it if not already included.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2023-03-12

# Set basic configuration variables for the BLAS library.
# 'name' is set to "BLAS", which is the identifier used for the dependency throughout the script.
set(name "BLAS")
# 'tag' specifies the version tag of the BLAS library to fetch, indicating a specific state of the source code in the repository.
set(tag "v0.3.21")
# 'version' specifies the version of the BLAS library. This may be used to ensure compatibility or meet specific requirements.
set(version "0.3.21")
# 'flag' can be used to pass additional flags to the configure/make commands when building the dependency, but it's empty here.
set(flag "")
# 'is_cmake' is a boolean flag indicating whether the BLAS library uses CMake for its build system. It's set to ON, meaning it does.
set(is_cmake ON)
# 'is_git' is a boolean flag indicating whether the BLAS library's source code is hosted in a git repository. It's set to ON.
set(is_git ON)
# 'auto_gen' is a boolean flag indicating whether to use autogen scripts for the configuration process. It's set to OFF here.
set(auto_gen OFF)
# 'url' specifies the location of the BLAS library's source code repository.
set(url "https://github.com/xianyi/OpenBLAS")

# Include the 'ImportDependency' macro script located in the 'macros' directory.
# This macro is responsible for importing and possibly installing the dependency.
include(macros/ImportDependency)
# Call the 'ImportDependency' macro with the previously set configuration parameters.
# This macro checks if BLAS is already available; if not, it proceeds to fetch, configure, build, and install it.
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})
# Print a message indicating the completion of the BLAS setup process.
message(STATUS "${name} done")
