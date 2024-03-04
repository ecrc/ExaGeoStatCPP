
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2024-02-14

# Define the target executable
add_executable(Example_Chameleon_Descriptor ChameleonDescriptor.cpp)

if (USE_HICMA)
    add_executable(Example_Hicma_Descriptor HicmaDescriptor.cpp)
    add_executable(Example_Chameleon_Hicma_Converter ChameleonToHicmaConverter.cpp)

    target_link_libraries(Example_Hicma_Descriptor PRIVATE ${PROJECT_NAME}_INTERFACE)
    target_link_libraries(Example_Chameleon_Hicma_Converter PRIVATE ${PROJECT_NAME}_INTERFACE)
endif ()

# Link the target executable with the project and any additional libraries
target_link_libraries(Example_Chameleon_Descriptor PRIVATE ${PROJECT_NAME}_INTERFACE)