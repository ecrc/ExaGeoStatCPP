
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief This file defines the build instructions for the Example_ExaGeoStat_Descriptor executable.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-07-17

# Define the executable target
add_executable(Example_ExaGeoStat_Descriptors ExaGeoStatDescriptorMain.cpp)

# Link the executable target to the ExaGeoStat library and any required external libraries
target_link_libraries(Example_ExaGeoStat_Descriptors ${PROJECT_NAME} ${LIBS} ${PROJECT_NAME}_INTERFACE)
