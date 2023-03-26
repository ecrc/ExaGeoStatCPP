
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-20

if(EXAGEOSTAT_USE_CHAMELEON)
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/dense/ChameleonImplementation.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/diagonal-super-tile/ChameleonImplementation.cpp
            ${SOURCES}
            PARENT_SCOPE
    )
endif()
if(EXAGEOSTAT_USE_HiCMA)
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/tile-low-rank/HicmaImplementation.cpp
            ${SOURCES}
            PARENT_SCOPE
            )
endif()