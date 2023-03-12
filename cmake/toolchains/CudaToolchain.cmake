
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CudaToolchain.cmake
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-12

# Set cuda compilation options.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CUDA_ARCHITECTURES "35;50;72")

find_package(CUDAToolkit REQUIRED)
