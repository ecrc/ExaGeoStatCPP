
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file BuildNLOPT.cmake
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-03-26

macro(BuildNLOPT raw_name url tag)
string(TOLOWER ${raw_name} name)
string(TOUPPER ${raw_name} capital_name)
message(STATUS "Fetching ${name} ${tag} from ${url}")
include(FetchContent)
set(FETCHCONTENT_BASE_DIR ${PROJECT_SOURCE_DIR}/installdir/_deps/NLOPT/)
FetchContent_Declare(${name} GIT_REPOSITORY "${url}" GIT_TAG "${tag}")
FetchContent_Populate(${name})
set(${name}_srcpath ${PROJECT_SOURCE_DIR}/installdir/_deps/NLOPT/${name}-src)
set(${name}_installpath ${PROJECT_SOURCE_DIR}/installdir/_deps/NLOPT/)

# Configure subproject into <subproject-build-dir>
include(ProcessorCount)
ProcessorCount(N)
execute_process(COMMAND sh autogen.sh
WORKING_DIRECTORY ${${name}_srcpath}
COMMAND_ERROR_IS_FATAL ANY)

execute_process(COMMAND ./configure --prefix=${${name}_installpath}
WORKING_DIRECTORY ${${name}_srcpath}
COMMAND_ERROR_IS_FATAL ANY)

execute_process(COMMAND make -j ${N}
WORKING_DIRECTORY ${${name}_srcpath}
COMMAND_ERROR_IS_FATAL ANY)

execute_process(COMMAND make install -j ${N}
WORKING_DIRECTORY ${${name}_srcpath}
COMMAND_ERROR_IS_FATAL ANY)

set(ENV{LD_LIBRARY_PATH} "${${name}_installpath}/lib:${${name}_installpath}/lib64:$ENV{LD_LIBRARY_PATH}")
set(ENV{LIBRARY_PATH} "${${name}_installpath}/lib:${${name}_installpath}/lib64:$ENV{LIBRARY_PATH}")
set(ENV{CPATH} "${${name}_installpath}/include:$ENV{CPATH}")
set(ENV{PKG_CONFIG_PATH} "${${name}_installpath}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
set(${capital_name}_DIR "${${name}_installpath}")
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
