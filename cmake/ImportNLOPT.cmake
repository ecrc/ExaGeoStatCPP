
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportNLOPT.cmake
# @brief Find and include NLOPT library as a dependency.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-26

# Configuration settings for the integration of the NLOPT library
# 'name' is assigned to "NLOPT", serving as the identifier for this library within the script.
set(name "NLOPT")
# 'tag' defines "v2.7.1" as the version tag of NLOPT, indicating the specific release to be utilized.
set(tag "v2.7.1")
# 'version' specifies "2.7.1" as the version of the NLOPT library, ensuring compatibility with the project's requirements.
set(version "2.7.1")
# 'flag' is intended for additional configuration options during the build process. Disable all language bindings to avoid Python compatibility issues.
set(flag -DNLOPT_PYTHON=OFF \-DNLOPT_SWIG=OFF \-DNLOPT_OCTAVE=OFF \-DNLOPT_MATLAB=OFF \-DNLOPT_GUILE=OFF)
# 'is_cmake' indicates that NLOPT uses CMake for its build system, which is set to ON.
set(is_cmake ON)
# 'is_git' denotes that the NLOPT source code is hosted in a Git repository, which is set to ON.
set(is_git ON)
# 'auto_gen' signals whether autogen scripts are required for the build process, which is set to OFF for NLOPT.
set(auto_gen OFF)
# 'url' provides the location of the NLOPT source code repository on GitHub.
set(url "https://github.com/stevengj/nlopt")

# The 'ImportDependency' macro script, located in the 'macros' directory, is included for managing the import and setup of the NLOPT library.
include(macros/ImportDependency)
# The 'ImportDependency' macro is invoked with the above-defined parameters to handle the detection, fetching, and integration of NLOPT into the project.
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

# A status message is outputted to indicate the successful integration of the NLOPT library into the project.
message(STATUS "${name} done")

