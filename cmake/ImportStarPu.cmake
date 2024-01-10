
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportSTARPU.cmake
# @brief Find and include STARPU library as a dependency.
# @version 1.0.1
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-13

#Configurations
set(name "STARPU")
set(tag "starpu-1.3.9")
set(version "1.3.9")

if (USE_CUDA AND USE_MPI)
    set(flag \--enable-cuda  \--disable-opencl  \--enable-shared  \--disable-build-doc  \--disable-export-dynamic  \--enable-mpi)
elseif(USE_CUDA)
    set(flag  \--enable-cuda  \--disable-opencl  \--enable-shared  \--disable-build-doc  \--disable-export-dynamic  \--disable-mpi)
elseif(USE_MPI)
    set(flag  \--disable-cuda  \--disable-opencl  \--enable-shared  \--disable-build-doc  \--disable-export-dynamic  \--enable-mpi)
else()
    set(flag \--disable-cuda  \--disable-opencl  \--enable-shared  \--disable-build-doc  \--disable-export-dynamic  \--disable-mpi)
endif()

set(is_cmake OFF)
set(is_git ON)
set(auto_gen ON)
set(url "https://gitlab.inria.fr/starpu/starpu.git")

include(macros/ImportDependency)
message(" ${STARPU_COMPONENT_LIST}")
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" ${STARPU_COMPONENT_LIST} ${is_cmake} ${is_git} ${auto_gen})

message(STATUS "${name} done")
