
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief Includes subdirectories for different modules of the ExaGeoStat software package.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2024-02-24

# Include subdirectories for end-to-end module, configurations module and data-generators module.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/configurations)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/data-generators)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/data-loader)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/descriptors)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/end-to-end)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/hardware)
