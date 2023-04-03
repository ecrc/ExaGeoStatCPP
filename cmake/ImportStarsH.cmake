
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-13

message("")
message("---------------------------------------- Stars-H")
message(STATUS "Checking for STARSH")
include(macros/BuildDependency)

if (NOT TARGET STARSH_FOUND)
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    find_package(STARSH QUIET)

    if (STARSH_FOUND)
        message("   Found STARSH: ${STARSH_LIBRARIES}")
    else()
        message("   Can't find STARSH, Installing it instead ..")
        set(FLAGS -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/installdir/_deps/STARSH/ \-DSTARPU=OFF \-DMPI=${USE_MPI})
        set(ISCMAKE ON)
        set(ISGIT ON)
        BuildDependency(STARSH "https://github.com/ecrc/stars-h.git" "v0.3.1" ${FLAGS} ${ISCMAKE} ${ISGIT})
        set(FLAGS "")
        find_package(STARSH REQUIRED)
    endif()
else()
    message("   STARSH already included")
endif()

include_directories(${STARSH_INCLUDE_DIRS_DEP})
if (STARSH_LINKER_FLAGS)
    list(APPEND CMAKE_EXE_LINKER_FLAGS "${STARSH_LINKER_FLAGS}")
endif ()
if (STARSH_LIBRARY_DIRS)
    # the RPATH to be used when installing
    list(APPEND CMAKE_INSTALL_RPATH "${STARSH_LIBRARY_DIRS}")
endif ()


if (STARSH_LIBRARIES)
    # look for gsl
    find_library(_STARSH_LIB NAME starsh PATHS ${STARSH_LIBRARY_DIRS})
    if (_STARSH_LIB AND NOT "${STARSH_LIBRARIES_DEP}" MATCHES "gsl")
        execute_process(COMMAND nm ${_STARSH_LIB} COMMAND grep gsl RESULT_VARIABLE GSL_IN_STARSH)
        if (${GSL_IN_STARSH} EQUAL 0)
            message(STATUS "STARSH depends on gsl. Adding it to dependency list")
            find_package(GSL REQUIRED)
            if (GSL_FOUND)
                if (STARSH_LIBRARIES_DEP)
                    list(APPEND STARSH_LIBRARIES_DEP ${GSL_LIBRARIES})
                else ()
                    list(APPEND STARSH_LIBRARIES ${GSL_LIBRARIES})
                endif ()
            endif ()
        endif ()
    endif ()
    # insert to dependencies
    if (STARSH_LIBRARIES_DEP)
        list(APPEND LIBS ${STARSH_LIBRARIES_DEP})
        link_directories(${STARSH_LIBRARY_DIRS_DEP})
        link_directories(${STARSH_LIBRARIES_DEP})
    else ()
        list(APPEND LIBS ${STARSH_LIBRARIES})
        link_directories(${STARSH_LIBRARIES})
    endif ()
    list(APPEND LIBS ${STARSH_LIBRARIES} )
    link_directories(${STARSH_LIBRARY_DIRS_DEP})
    include_directories(${STARSH_INCLUDE_DIRS})

endif()

message(STATUS "StarsH Done")
