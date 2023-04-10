
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportGFortran.cmake
# @brief Defines the gfortran library for use in the project.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-14

# Export the list of libraries, including gfortran, for use in the project.
set(LIBS
        gfortran
        ${LIBS}
        )
