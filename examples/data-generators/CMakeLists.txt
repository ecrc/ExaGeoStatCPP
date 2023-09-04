
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief CMake file for building the Example_Synthetic_Data_Generation executable.
# @version 1.0.0
# @author Mahmoud ElKarargy
# @date 2023-03-04

# Define the target executable
add_executable(Example_Synthetic_Data_Generation SyntheticLocationsGeneration.cpp)

# Link the target executable with the project and any additional libraries
target_link_libraries(Example_Synthetic_Data_Generation PRIVATE ${PROJECT_NAME} ${LIBS} ${PROJECT_NAME}_INTERFACE)
