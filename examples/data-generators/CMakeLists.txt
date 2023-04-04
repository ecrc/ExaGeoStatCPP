
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-04

add_executable(Example_Synthetic_Data_Generation SyntheticData_main.cpp)
target_link_libraries(Example_Synthetic_Data_Generation ${PROJECT_NAME} ${LIBS})