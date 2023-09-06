
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportCuSolver.cmake
# @brief This script sets CUDA libraries and adds CuSolver, CuBlas, and CuBlasLt to the list of libraries.
# @version 1.0.0
# @author Mahmoud ElKarargy
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
