
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportCuSolver.cmake
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-13

set(cudart_lib CUDA::cudart)
set(cublas_lib CUDA::cublas)
set(LIBS
        CUDA::cusolver
        CUDA::cublas
        CUDA::cublasLt
        ${LIBS}
        )
