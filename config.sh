#!/bin/bash
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (C) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file config.sh
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-01-30

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m'

INSTALL_PREFIX=$PWD/bin/installdir
PROJECT_SOURCE_DIR=$(dirname "$0")

while getopts ":tevhHCi:dc" opt; do
  case $opt in
    i) ##### Define installation path  #####
       echo -e "${YELLOW}Installation path set to $OPTARG.${NC}"
       INSTALL_PREFIX=$OPTARG
       ;;
    t) ##### Building tests enabled #####
      echo -e "${GREEN}Building tests enabled.${NC}"
      BUILDING_TESTS="ON"
      ;;
    e) ##### Building examples enabled #####
      echo -e "${GREEN}Building examples enabled.${NC}"
      BUILDING_EXAMPLES="ON"
      ;;
    H) ##### Using HiCMA #####
      echo -e "${GREEN}Using HiCMA.${NC}"
      USING_HiCMA="ON"
      ;;
    C) ##### Using Chameleon #####
        echo -e "${GREEN}Using Chameleon.${NC}"
        USING_CHAMELEON="ON"
        ;;
    c)##### Using cuda enabled #####
        echo -e "${GREEN}Cuda enabled ${NC}"
        USE_CUDA=ON
        ;;
    v) ##### printing full output of make #####
      echo -e "${YELLOW}printing make with details.${NC}"
      VERBOSE=ON
      ;;
    d)##### Using debug mode to build #####
      echo -e "${RED}Debug mode enabled ${NC}"
      BUILD_TYPE="DEBUG"
      ;;
    \?) ##### using default settings #####
      BUILDING_TESTS="OFF"
      BUILDING_EXAMPLES="OFF"
      USING_HiCMA="OFF"
      USING_CHAMELEON="OFF"
      INSTALL_PREFIX=$PWD/bin/installdir
      VERBOSE=OFF
      BUILD_TYPE="RELEASE"
      USE_CUDA="OFF"

      echo -e "${RED}Building tests disabled.${NC}"
      echo -e "${RED}Building examples disabled.${NC}"
      echo -e "${BLUE}Installation path set to $INSTALL_PREFIX.${NC}"
      echo -e "${BLUE}Using HiCMA is disabled.${NC}"
      echo -e "${BLUE}Using Chameleon is disabled.${NC}"
      ;;
    :) ##### Error in an option #####
      echo "Option $OPTARG requires parameter(s)"
      exit 0
      ;;
    h) ##### Prints the help #####
      echo "Usage of $(basename "$0"):"
      echo ""
      printf "%20s %s\n" "-t :" "to enable building tests."
      echo ""
      printf "%20s %s\n" "-e :" "to enable building examples."
      echo ""
      printf "%20s %s\n" "-i [path] :" "specify installation path."
      printf "%20s %s\n" "" "default = /exageostat-cpp/bin/installdir"
      echo ""
      printf "%20s %s\n" "-H :" "to enable using HiCMA."
      echo ""
      printf "%20s %s\n" "-C :" "to enable using chameleon."
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

if [ -z "$USING_CHAMELEON" ]; then
   echo -e "${RED}Using Chameleon is disabled.${NC}"
fi

if [ -z "$BUILD_TYPE" ]; then
  BUILD_TYPE="RELEASE"
  echo -e "${GREEN}Building in release mode${NC}"
fi
if [ -z "$USE_CUDA" ]; then
  USE_CUDA="OFF"
  echo -e "${RED}Using CUDA disabled${NC}"
fi


echo ""
echo -e "${YELLOW}Use -h to print the usages of exageostat-cpp flags.${NC}"
echo ""
rm -rf bin/
mkdir -p bin/installdir

cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DEXAGEOSTAT_INSTALL_PREFIX="$INSTALL_PREFIX" \
  -DEXAGEOSTAT_BUILD_TESTS="$BUILDING_TESTS" \
  -DEXAGEOSTAT_BUILD_EXAMPLES="$BUILDING_EXAMPLES" \
  -DEXAGEOSTAT_USE_HICMA="$USING_HiCMA" \
  -DCMAKE_VERBOSE_MAKEFILE:BOOL=$VERBOSE \
  -DUSE_CUDA="$USE_CUDA"\
  -H"${PROJECT_SOURCE_DIR}" \
  -B"${PROJECT_SOURCE_DIR}/bin"
