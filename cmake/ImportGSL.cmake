
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportGSL.cmake
# @brief Checks for the GSL library and includes it in the project if it is not already present.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-16

#Configurations
set(name "GSL")
set(tag "v2.7.1")
set(version "2.7.1")
set(flag "")
set(is_cmake OFF)
set(is_git OFF)
set(auto_gen OFF)
set(url "https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz")

include(macros/ImportDependency)
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

# Add the GSL library to the project's list of libraries.
list(APPEND LIBS gsl)

message(STATUS "${name} done")
