
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief Find and include STARSH library as a dependency.
# @version 1.0.1
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-13

#Configurations
set(name "STARSH")
set(tag "v0.3.1")
set(version "0.3.1")
set(flag \-DSTARPU=OFF \-DMPI=${USE_MPI})
set(is_cmake ON)
set(is_git ON)
set(auto_gen OFF)
set(url "https://github.com/ecrc/stars-h.git")

include(macros/ImportDependency)
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

message(STATUS "${name} done")


message("")
message("---------------------------------------- Stars-H AGAIN WITH OLD CMAKE")
find_package(STARSH)
if (NOT STARSH_FOUND)

    message("")
    message(STATUS "Can't find STARSH using PkgConfig, Installing it through HiCMA")
    include(ImportHiCMA)
    if(NOT HiCMA_INSTALLED)
        message(FATAL_ERROR "HiCMA Installation failed")
    endif()

    message(STATUS "Now trying to find STARSH")
endif()
find_package(STARSH REQUIRED)
if (STARSH_FOUND)
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
                if (STARSH_LIBRARIES_DEP)
                    list(APPEND STARSH_LIBRARIES_DEP "gsl")
                else ()
                    list(APPEND STARSH_LIBRARIES "gsl")
                endif ()
            endif ()
        endif ()
        # insert to dependencies
        if (STARSH_LIBRARIES_DEP)
            list(INSERT EXAGEOSTAT_DEP 0 ${STARSH_LIBRARIES_DEP})
        else ()
            list(INSERT EXAGEOSTAT_DEP 0 ${STARSH_LIBRARIES})
        endif ()
    endif ()
endif()
message(STATUS "StarsH Done")