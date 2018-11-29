#!/bin/bash
#B136013
#
#define variables
DATA_FOLDER=./data
BUILD_FOLDER=./build
TIMING_HEADER="Input File\tProcesses Number\tAverage Iteration Time (ms)"

#make required folders
mkdir -p $DATA_FOLDER
mkdir -p $BUILD_FOLDER

#make required subfolders
mkdir -p $DATA_FOLDER/average_pixel_prod/

#clean subfolders
rm -rf $DATA_FOLDER/average_pixel_prod/*

#clean build files
make clean_build

# make production mode executable
make imagenew

# running production mode executable
echo "running production mode executable..."
echo -e $TIMING_HEADER > $DATA_FOLDER/time_results_prod.tsv
for RESOURCE in ./resources/*.pgm
do
  mpirun -n 16 ./imagenew $RESOURCE
done

python test.py
