
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportSTARPU.cmake
# @brief Find and include STARPU library as a dependency.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-13

# Configuration settings for integrating the STARPU library into the project
# 'name' is set to "STARPU" to identify this specific library within the script.
set(name "STARPU")
# 'tag' specifies "starpu-1.3.10" as the version tag, indicating the exact version of STARPU to be used.
set(tag "starpu-1.3.10")
# 'version' sets "1.3.10" as the version of the STARPU library, ensuring project compatibility.
set(version "1.3.10")

# Conditional setting of 'flag' based on project configurations for CUDA and MPI.
if (USE_CUDA AND USE_MPI)
    # Sets flags for enabling CUDA and MPI, and disables OpenCL, documentation build, and export dynamic when both CUDA and MPI are used.
    set(flag \--enable-cuda  \--disable-opencl  \--enable-shared  \--disable-build-doc  \--disable-export-dynamic  \--enable-mpi)
elseif(USE_CUDA)
    # Sets flags for enabling CUDA and shared libraries, and disables OpenCL, documentation build, export dynamic, and MPI when only CUDA is used.
    set(flag  \--enable-cuda  \--disable-opencl  \--enable-shared  \--disable-build-doc  \--disable-export-dynamic  \--disable-mpi)
elseif(USE_MPI)
    # Sets flags for enabling MPI and shared libraries, and disables CUDA, OpenCL, documentation build, and export dynamic when only MPI is used.
    set(flag  \--disable-cuda  \--disable-opencl  \--enable-shared  \--disable-build-doc  \--disable-export-dynamic  \--enable-mpi)
else()
    # Sets flags for enabling shared libraries, and disables CUDA, OpenCL, documentation build, export dynamic, and MPI when neither CUDA nor MPI is used.
    set(flag \--disable-cuda  \--disable-opencl  \--enable-shared  \--disable-build-doc  \--disable-export-dynamic  \--disable-mpi)
endif()

# 'is_cmake' is set to OFF indicating STARPU does not use CMake as its primary build system.
set(is_cmake OFF)
# 'is_git' is set to ON, denoting that STARPU's source code is maintained in a Git repository.
set(is_git ON)
# 'auto_gen' is set to ON, signaling the need for autogen scripts to be run as part of the build process.
set(auto_gen ON)
# 'url' provides the location of the STARPU source code repository.
set(url "https://gitlab.inria.fr/starpu/starpu.git")

# Include the 'ImportDependency' macro script, responsible for handling the import and setup of dependencies.
include(macros/ImportDependency)
# The 'ImportDependency' macro is invoked with the configuration parameters to manage the detection, fetching, and integration of STARPU.
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "${STARPU_COMPONENT_LIST}" ${is_cmake} ${is_git} ${auto_gen})

# A message is logged to indicate the successful integration of the STARPU library into the project.
message(STATUS "${name} done")

