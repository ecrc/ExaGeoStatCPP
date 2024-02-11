
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief Defines an executable and links it with the ExaGeoStat library and other libraries.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2023-01-31

# Define an executable named "Example_Configurations".
add_executable(Example_Configurations_Setup ${CMAKE_CURRENT_SOURCE_DIR}/SetupConfigurations.cpp)
add_executable(Example_Different_Configurations_Running ${CMAKE_CURRENT_SOURCE_DIR}/RunningWithDifferentConfigurations.cpp)

# Link the executable with the ExaGeoStat library and other libraries.
target_link_libraries(Example_Configurations_Setup ${PROJECT_NAME}_INTERFACE)
target_link_libraries(Example_Different_Configurations_Running ${PROJECT_NAME}_INTERFACE)