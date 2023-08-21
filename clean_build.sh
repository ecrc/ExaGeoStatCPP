#!/bin/bash
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file clean_build.sh
# @brief This script cleans and builds a software package called ExaGeoStat.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-01-30

# Define variables.
verbose=""
num_proc="-j $(nproc)"  # Use the number of available processors by default.

# Parse command-line arguments.
while getopts "vj:h" opt; do
  case $opt in
    v)
      verbose="VERBOSE=1"
      echo "Using verbose mode"
      ;;
    j)
      num_proc="-j $OPTARG"
      echo "Using $OPTARG threads to build"
      ;;
    h)
      # Print help information and exit.
      echo "Usage: $(basename "$0") [-v] [-j <thread_number>] [-h]"
      echo "Clean and build the ExaGeoStat software package."
      echo ""
      echo "Options:"
      echo "  -v                 Use verbose output."
      echo "  -j <thread_number> Build with a specific number of threads."
      echo "  -h                 Show this help message."
      exit 0
      ;;
    *)
      # Print an error message and exit.
      echo "Invalid flag. Use the -h flag for help."
      exit 1
      ;;
  esac
done

# Change to the bin directory, or exit if it doesn't exist.
cd bin/ || {
  echo "Error: bin directory not found."
  exit 1
}

# Clean the directory and build the code with the specified options.
make clean
make all $num_proc $verbose