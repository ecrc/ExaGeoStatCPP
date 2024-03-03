
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2023-01-31

# Define the target executable
add_executable(Example_Data_Generation ${CMAKE_CURRENT_SOURCE_DIR}/DataGeneration.cpp)
add_executable(Example_Data_Modeling ${CMAKE_CURRENT_SOURCE_DIR}/DataModeling.cpp)
add_executable(Example_Data_Generation_and_Modeling ${CMAKE_CURRENT_SOURCE_DIR}/DataGenerationAndModeling.cpp)
add_executable(Example_Data_Generation_Modeling_and_Prediction ${CMAKE_CURRENT_SOURCE_DIR}/DataGenerationModelingAndPrediction.cpp)
add_executable(Example_Data_Prediction ${CMAKE_CURRENT_SOURCE_DIR}/DataPrediction.cpp)
add_executable(Example_Data_Generation_and_Prediction ${CMAKE_CURRENT_SOURCE_DIR}/DataGenerationAndPrediction.cpp)

# Link the target executable with the project and any additional libraries
target_link_libraries(Example_Data_Generation PUBLIC ${PROJECT_NAME}_INTERFACE)
target_link_libraries(Example_Data_Modeling PUBLIC ${PROJECT_NAME}_INTERFACE)
target_link_libraries(Example_Data_Generation_and_Modeling PUBLIC ${PROJECT_NAME}_INTERFACE)
target_link_libraries(Example_Data_Generation_Modeling_and_Prediction PUBLIC ${PROJECT_NAME}_INTERFACE)
target_link_libraries(Example_Data_Prediction PUBLIC ${PROJECT_NAME}_INTERFACE)
target_link_libraries(Example_Data_Generation_and_Prediction PUBLIC ${PROJECT_NAME}_INTERFACE)
