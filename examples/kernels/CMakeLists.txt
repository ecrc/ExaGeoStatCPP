
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief CMake file for building the Examples on kernels executable.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-04-12

# Define the target executable
add_executable(Example_kernels Kernels_main.cpp)

# Link the target executable with the project and any additional libraries
target_link_libraries(Example_kernels PRIVATE ${PROJECT_NAME} ${LIBS})