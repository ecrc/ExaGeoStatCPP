
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.1
# @author Mahmoud ElKarargy
# @date 2024-02-14

# Define the target executable
add_executable(Example_Hardware ${CMAKE_CURRENT_SOURCE_DIR}/ExaGeoStatHardware.cpp)

# Link the target executable with the project and any additional libraries
target_link_libraries(Example_Hardware PRIVATE ${PROJECT_NAME}_INTERFACE)

