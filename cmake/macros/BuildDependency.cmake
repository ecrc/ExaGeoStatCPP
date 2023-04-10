# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file BuildDependency.cmake
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-12

macro(BuildDependency raw_name url tag ${FLAGS} ${ISCMAKE} ${ISGIT})
    # Set the a lower case and upper case name of the dependency.
    string(TOLOWER ${raw_name} name)
    string(TOUPPER ${raw_name} capital_name)
    message(STATUS "Fetching ${name} ${tag} from ${url}")
    include(FetchContent)
    set(FETCHCONTENT_BASE_DIR ${PROJECT_SOURCE_DIR}/installdir/_deps/${capital_name}/)
    # Fetch the dependency, depending on which it's a git repo or no.
    if (ISGIT)
        FetchContent_Declare(${name} GIT_REPOSITORY "${url}" GIT_TAG "${tag}")
    else()
        FetchContent_Declare(${name} URL "${url}")
    endif ()
    FetchContent_Populate(${name})
    # Installation of the source files will be in bin/_deps/${name}-src
    set(${name}_srcpath ${PROJECT_SOURCE_DIR}/installdir/_deps/${capital_name}/${name}-src)
    # The bin directory where the code will get build is also in bin/_deps/${name}_srcpath/bin
    set(${name}_binpath ${${name}_srcpath}/bin)
    # The installation will be installdir/capital_name/ .To avoid deleting it when building software multiple time
    set(${name}_installpath ${PROJECT_SOURCE_DIR}/installdir/_deps/${capital_name}/)
    file(MAKE_DIRECTORY ${${name}_binpath})

    # Configure subproject into <subproject-build-dir>
    if (ISCMAKE)
        execute_process(COMMAND ${CMAKE_COMMAND} ${FLAGS}
                ${${name}_srcpath}
                WORKING_DIRECTORY
                ${${name}_binpath})
    else()
        execute_process(COMMAND ./configure ${FLAGS}
                WORKING_DIRECTORY ${${name}_srcpath}
                COMMAND_ERROR_IS_FATAL ANY)
    endif ()

    # Build and install subproject
    include(ProcessorCount)
    ProcessorCount(N)

    if (ISCMAKE)
        execute_process(COMMAND make -j ${N}
                WORKING_DIRECTORY ${${name}_binpath}
                COMMAND_ERROR_IS_FATAL ANY)

        execute_process(COMMAND make install -j ${N}
                WORKING_DIRECTORY ${${name}_binpath}
                COMMAND_ERROR_IS_FATAL ANY)
    else()
        execute_process(COMMAND make -j ${N}
                WORKING_DIRECTORY ${${name}_srcpath}
                COMMAND_ERROR_IS_FATAL ANY)

        execute_process(COMMAND make install -j ${N}
                WORKING_DIRECTORY ${${name}_srcpath}
                COMMAND_ERROR_IS_FATAL ANY)
    endif ()

    set(ENV{LD_LIBRARY_PATH} "${${name}_installpath}/lib:${${name}_installpath}/lib64:$ENV{LD_LIBRARY_PATH}")
    set(ENV{LIBRARY_PATH} "${${name}_installpath}/lib:${${name}_installpath}/lib64:$ENV{LIBRARY_PATH}")
    set(ENV{CPATH} "${${name}_installpath}/include:$ENV{CPATH}")
    set(ENV{PKG_CONFIG_PATH} "${${name}_installpath}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    set(${capital_name}_DIR "${${name}_installpath}/lib")
    include_directories(${${name}_installpath}/include)
    link_directories(${${name}_installpath}/lib)
    install(
            DIRECTORY
            "${${name}_installpath}/lib"
            DESTINATION
            .
    )
    install(
            DIRECTORY
            "${${name}_installpath}/include"
            DESTINATION
            .
    )
    install(
            DIRECTORY
            "${${name}_installpath}/share"
            DESTINATION
            .
    )
endmacro()