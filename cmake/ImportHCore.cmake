
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportHCore.cmake
# @brief Checks for the Hcore library and includes it in the project if it is not already present.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-15

# Configurations for the HCORE library
# 'name' is set to "HCORE", serving as the identifier for the dependency throughout this script.
set(name "HCORE")
# 'tag' defines the version tag of the HCORE library to be fetched, represents a specific snapshot of the source code in the repository.
set(tag "v0.1.3")
# 'version' indicates the version of the HCORE library, used to ensure compatibility or fulfill certain requirements.
set(version "0.1.3")
# 'flag' is intended for passing additional flags to the configure/make commands during the building of the dependency, but is left empty here.
set(flag "")
# 'is_cmake' is a boolean flag that denotes whether the HCORE library utilizes CMake for its build process. Here, it's set to ON.
set(is_cmake ON)
# 'is_git' is a boolean flag that signifies if the source code for HCORE is maintained in a git repository. This is set to ON.
set(is_git ON)
# 'auto_gen' is a boolean flag that indicates whether autogen scripts are required for the configuration process. It is set to OFF in this context.
set(auto_gen OFF)
# 'url' provides the location of the HCORE library's source code repository.
set(url "https://github.com/ecrc/hcore.git")

# The 'ImportDependency' macro script, located in the 'macros' directory, is included. This macro is crucial for importing and potentially installing the dependency.
include(macros/ImportDependency)
# The 'ImportDependency' macro is invoked with the configuration parameters set above to handle the detection, fetching, and setup of HCORE.
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})
# A message is logged to indicate the successful completion of the HCORE setup process.
message(STATUS "${name} done")

