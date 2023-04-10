
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief Configures and generates documentation using Doxygen.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-12

# Configure the Doxygen configuration file from the template.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.in"
        "${CMAKE_CURRENT_BINARY_DIR}/config")

# Generate documentation using Doxygen.
add_custom_target(docs COMMAND
        "${DOXYGEN_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/config")