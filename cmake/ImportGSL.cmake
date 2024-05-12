
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportGSL.cmake
# @brief Checks for the GSL library and includes it in the project if it is not already present.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-16

# Configurations for integrating the GSL (GNU Scientific Library)
# 'name' is assigned "GSL" to identify the GNU Scientific Library within this script.
set(name "GSL")
# 'tag' specifies the version tag "v2.7.1" for the GSL library, indicating the exact version to be used.
set(tag "v2.7.1")
# 'version' sets "2.7.1" as the version of the GSL library, ensuring compatibility with project requirements.
set(version "2.7.1")
# 'flag' is available for additional configuration options during build or installation, but remains empty here.
set(flag "")
# 'is_cmake' indicates whether GSL uses CMake for building. It is set to OFF, implying an alternative build system is used.
set(is_cmake OFF)
# 'is_git' denotes if GSL's source code is hosted in a Git repository. It is set to OFF, suggesting the source is obtained from a different location.
set(is_git OFF)
# 'auto_gen' signifies the need for autogen scripts in the build process. Here, it is set to OFF, indicating they are not needed.
set(auto_gen OFF)
# 'url' provides the download location for the GSL source code, pointing to the GNU archive.
set(url "https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz")

# Include the 'ImportDependency' macro, responsible for managing the GSL library's import and setup process.
include(macros/ImportDependency)
# Execute the 'ImportDependency' macro with the previously established parameters to handle the detection, downloading, and integration of GSL.
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

# Add GSL to the project's list of linked libraries, making its functionality accessible within the project.
list(APPEND LIBS gsl)

# Output a message signaling the successful integration of the GSL library into the project.
message(STATUS "${name} done")

