
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief Defines an executables and links them with the ExaGeoStat library and other libraries.
# @version 2.0.0
# @author Mahmoud ElKarargy
# @date 2024-09-23

# Define the target executable
add_executable(Example_Stage_Zero ${CMAKE_CURRENT_SOURCE_DIR}/StageZeroDataGenerator.cpp)

# Link the target executable with the project and any additional libraries
target_link_libraries(Example_Stage_Zero PUBLIC ${PROJECT_NAME}_INTERFACE)
