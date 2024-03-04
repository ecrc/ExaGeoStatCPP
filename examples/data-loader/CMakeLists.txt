
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2024-02-14

# Define the target executable
add_executable(Example_CSV_Loader ${CMAKE_CURRENT_SOURCE_DIR}/CSVLoader.cpp)

# Link the target executable with the project and any additional libraries
target_link_libraries(Example_CSV_Loader PRIVATE ${PROJECT_NAME}_INTERFACE)
