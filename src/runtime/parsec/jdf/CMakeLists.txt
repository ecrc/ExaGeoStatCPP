
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @brief CMake build script for jdf directory
# @author Mahmoud ElKarargy
# @date 2024-10-17

# Add .jdf files
set(JDF_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/ReadCSVTimeSlot.jdf
        ${CMAKE_CURRENT_SOURCE_DIR}/ReadCSVComplex.jdf
        ${CMAKE_CURRENT_SOURCE_DIR}/ReadCSVToComplex.jdf
        ${CMAKE_CURRENT_SOURCE_DIR}/ReadCSVToComplexTimeSlot.jdf
        ${CMAKE_CURRENT_SOURCE_DIR}/ReadCSV.jdf
        ${CMAKE_CURRENT_SOURCE_DIR}/ForwardSHT.jdf
        ${CMAKE_CURRENT_SOURCE_DIR}/MeanSquaredError.jdf
        ${CMAKE_CURRENT_SOURCE_DIR}/DifferenceDouble.jdf
        ${CMAKE_CURRENT_SOURCE_DIR}/ForwardSHTReshape.jdf
        ${CMAKE_CURRENT_SOURCE_DIR}/GetMatrixNorm.jdf
        ${CMAKE_CURRENT_SOURCE_DIR}/InverseSHT.jdf
        ${CMAKE_CURRENT_SOURCE_DIR}/MatrixCompress.jdf
)

# Output directory for generated sources
set(JDF_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated_jdf")

# Make sure the output directory exists
file(MAKE_DIRECTORY ${JDF_OUTPUT_DIR})

foreach (jdf_file ${JDF_FILES})
    get_filename_component(jdf_name ${jdf_file} NAME_WE)
    execute_process(
            COMMAND ${HICMA_X_BIN_DIR}/parsec-ptgpp -i ${jdf_file} -o ${jdf_name}
            ERROR_QUIET
    )
    execute_process(
            COMMAND mv ${jdf_name}.c ${CMAKE_CURRENT_BINARY_DIR}/generated_jdf/${jdf_name}.c
            COMMAND mv ${jdf_name}.h ${CMAKE_CURRENT_BINARY_DIR}/generated_jdf/${jdf_name}.h
    )
    list(APPEND JDF_GENERATED_SOURCES ${JDF_OUTPUT_DIR}/${jdf_name}.c)
endforeach ()

set(JDF_GENERATED_SOURCES
        ${JDF_GENERATED_SOURCES}
        PARENT_SCOPE
)
