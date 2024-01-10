#!/bin/bash
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file config.sh
# @version 1.0.1
# @author Mahmoud ElKarargy
# @date 2023-01-30

# Set variables and default values
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m'

INSTALL_PREFIX=$PWD/installdir/_deps
PROJECT_SOURCE_DIR=$(dirname "$0")
BUILDING_TESTS="OFF"
BUILDING_HEAVY_TESTS="OFF"
BUILDING_EXAMPLES="OFF"
USING_HiCMA="OFF"
VERBOSE="OFF"
USE_CUDA="OFF"
USE_MPI="OFF"
BLAS_VENDOR=""
PACKAGE="OFF"

# Parse command line options
while getopts ":tevhHi:cmspT" opt; do
  case $opt in
    i) ##### Define installation path  #####
       echo -e "${YELLOW}Installation path set to $OPTARG.${NC}"
       INSTALL_PREFIX=$OPTARG
       ;;
    t) ##### Building tests enabled #####
      echo -e "${GREEN}Building tests enabled.${NC}"
      BUILDING_TESTS="ON"
      ;;
    T) ##### Building heavy tests enabled #####
      echo -e "${GREEN}Building heavy tests enabled.${NC}"
      BUILDING_HEAVY_TESTS="ON"
      ;;
    e) ##### Building examples enabled #####
      echo -e "${GREEN}Building examples enabled.${NC}"
      BUILDING_EXAMPLES="ON"
      ;;
    H) ##### Using HiCMA #####
      echo -e "${GREEN}Using HiCMA.${NC}"
      USING_HiCMA="ON"
      ;;
    c)##### Using cuda enabled #####
        echo -e "${GREEN}Cuda enabled ${NC}"
        USE_CUDA=ON
        ;;
    m)##### Using MPI enabled #####
        echo -e "${GREEN}MPI enabled ${NC}"
        USE_MPI=ON
        ;;
    v) ##### printing full output of make #####
      echo -e "${GREEN}printing make with details.${NC}"
      VERBOSE=ON
      ;;
    s) ##### Passing BLA vendor with mkl #####
      echo -e "${GREEN}MKL as a BLA vendor${NC}"
      BLAS_VENDOR="Intel10_64lp"
      ;;
    p) ##### Enabling packaging system for distribution #####
      echo -e "${GREEN}CPACK enabled${NC}"
      PACKAGE=ON
      ;;
    \?) ##### Error unknown option #####
      echo "Option $OPTARG parameter is unknown, please -h for help"
      exit 1
      ;;
    :) ##### Error in an option #####
      echo "Option $OPTARG requires parameter(s)"
      exit 0
      ;;
    h) ##### Prints the help #####
      echo "Usage of $(basename "$0"):"
      echo ""
      printf "%20s %s\n" "-i [path] :" "specify installation path, default = ${PWD}/installdir/_deps/"
      printf "%20s %s\n" "-t :" "to enable building tests."
      printf "%20s %s\n" "-T :" "to enable building heavy tests."
      printf "%20s %s\n" "-e :" "to enable building examples."
      printf "%20s %s\n" "-H :" "to enable using HiCMA."
      printf "%20s %s\n" "-c :" "to enable using CUDA."
      printf "%20s %s\n" "-m :" "to enable using MPI."
      printf "%20s %s\n" "-v :" "to enable verbose printings."
      printf "%20s %s\n" "-d :" "to enable debug mode."
      printf "%20s %s\n" "-s :" "to manually pass MKL as your bla vendor."
      printf "%20s %s\n" "-p :" "to enable a packaging system for distribution."
      printf "%20s %s\n" "-h :" "Help."
      echo ""
      exit 1
      ;;
    esac
done

if [ -z "$BUILDING_TESTS" ]; then
  BUILDING_TESTS="OFF"
  echo -e "${RED}Building tests disabled.${NC}"
fi

if [ -z "$BUILDING_EXAMPLES" ]; then
  BUILDING_EXAMPLES="OFF"
  echo -e "${RED}Building examples disabled.${NC}"
fi

echo -e "${BLUE}Installation path set to $INSTALL_PREFIX.${NC}"

if [ -z "$USING_HiCMA" ]; then
   echo -e "${RED}Using HiCMA is disabled.${NC}"
fi

if [ -z "$USE_CUDA" ]; then
  USE_CUDA="OFF"
  echo -e "${RED}Using CUDA disabled${NC}"
fi

if [ -z "$USE_MPI" ]; then
  USE_MPI="OFF"
  echo -e "${RED}Using MPI disabled${NC}"
fi

echo ""
echo -e "${YELLOW}Use -h to print the usages of exageostat-cpp flags.${NC}"
echo ""
rm -rf bin/
mkdir -p bin/installdir

cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
  -DBUILD_TESTS="${BUILDING_TESTS}" \
  -DBUILD_HEAVY_TESTS="${BUILDING_HEAVY_TESTS}" \
  -DBUILD_EXAMPLES="${BUILDING_EXAMPLES}" \
  -DUSE_HICMA="${USING_HiCMA}" \
  -DCMAKE_VERBOSE_MAKEFILE:BOOL=${VERBOSE} \
  -DUSE_CUDA="${USE_CUDA}" \
  -DUSE_MPI="${USE_MPI}" \
  -DBLA_VENDOR="${BLAS_VENDOR}" \
  -DCREATE_PACKAGE="${PACKAGE}" \
  -H"${PROJECT_SOURCE_DIR}" \
  -B"${PROJECT_SOURCE_DIR}/bin" \
  -G "Unix Makefiles"
