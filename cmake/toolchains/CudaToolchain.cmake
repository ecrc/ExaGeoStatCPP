
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CudaToolchain.cmake
# @brief This file is used to set up the CUDA toolchain for compilation.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-12

# Set CUDA compilation options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the CUDA architectures to be targeted
set(CUDA_ARCHITECTURES "35;50;72")

# Find the CUDA toolkit
find_package(CUDAToolkit REQUIRED)
set(ENV{LDFLAGS} "-L$ENV{CUDA_DIR}/lib64")
