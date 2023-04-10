
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief Defines an executable and links it with the ExaGeoStat library and other libraries.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-01-31

# Define an executable named "Example_Synthetic_Data_Configurations".
add_executable(Example_Synthetic_Data_Configurations ${CMAKE_CURRENT_SOURCE_DIR}/data-generation/main_synthetic.cpp)

# Link the executable with the ExaGeoStat library and other libraries.
target_link_libraries(Example_Synthetic_Data_Configurations ${PROJECT_NAME} ${LIBS})