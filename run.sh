#!/bin/bash
#B136013
#
#define variables
DATA_FOLDER=./data
GRAPHS_FOLDER=./graphs
OUTPUT_FOLDER=./output
BUILD_FOLDER=./build
TIMING_HEADER="Input File\tProcesses Number\tAverage Iteration Time (ms)"
RESOURCES_LIST="./resources/edgenew192x128.pgm_./resources/edgenew256x192.pgm_./resources/edgenew512x384.pgm_./resources/edgenew768x768.pgm_./resources/edgenew1600x1200.pgm"
RESOURCES_BIG_LIST="./resources/edgenew768x768.pgm_./resources/edgenew1600x1200.pgm"
# RESERVATION="-q R380254"

#make required folders
mkdir -p $DATA_FOLDER
mkdir -p $GRAPHS_FOLDER
mkdir -p $BUILD_FOLDER

#clean folders
rm -rf $DATA_FOLDER/*
rm -rf $GRAPHS_FOLDER/*
rm -rf $OUTPUT_FOLDER/*

#make required subfolders
mkdir -p $DATA_FOLDER/average_pixel_prod/
mkdir -p $DATA_FOLDER/average_pixel_test/
mkdir -p $DATA_FOLDER/average_pixel_intervals_test/
mkdir -p $GRAPHS_FOLDER/average_pixel/

#clean subfolders
rm -rf $DATA_FOLDER/average_pixel_prod/*
rm -rf $DATA_FOLDER/average_pixel_test/*
rm -rf $DATA_FOLDER/average_pixel_intervals_test/*
rm -rf $GRAPHS_FOLDER/average_pixel/*

#clean executables and build files
make clean

# check average pixel with intervals in each iteration
make clean_build
make imagenew_average_pixel
echo "checking average pixel with intervals in each iteration..."
NODES_NUM=select=1:ncpus=36
PROC_LIST="16"
echo -e $TIMING_HEADER > $DATA_FOLDER/time_results_average_pixel_test.tsv
qsub $RESERVATION -v PROC_LIST="${PROC_LIST}",RESOURCES_LIST="${RESOURCES_LIST}",MPIPROG="imagenew_average_pixel" -l $NODES_NUM imagenew.pbs

# check average iteration time
make clean_build
make imagenew_timing
echo "checking average iteration time..."
PROC_LIST="1_2_3_4_5_7_8_11_12_16_20_24_28_32_36"
echo -e $TIMING_HEADER > $DATA_FOLDER/time_results_timing_test.tsv
qsub $RESERVATION -v PROC_LIST="${PROC_LIST}",RESOURCES_LIST="${RESOURCES_LIST}",MPIPROG="imagenew_timing" -l $NODES_NUM imagenew.pbs

# check average iteration time for big input
echo "checking average iteration time for big input..."
NODES_NUM=select=4:ncpus=36
PROC_LIST="40_48_52_60_64_72_80_88_96_104_112_116_128_136"
qsub $RESERVATION -v PROC_LIST="${PROC_LIST}",RESOURCES_LIST="${RESOURCES_BIG_LIST}",MPIPROG="imagenew_timing" -l $NODES_NUM imagenew.pbs

# check average iteration time with intervals in each iteration
make clean_build
make imagenew_timing_intervals
echo "checking average iteration time with intervals in each iteration..."
echo -e $TIMING_HEADER > $DATA_FOLDER/time_results_timing_intervals_test.tsv
NODES_NUM=select=1:ncpus=36
PROC_LIST="1_2_3_4_5_7_8_11_12_16"
qsub $RESERVATION -v PROC_LIST="${PROC_LIST}",RESOURCES_LIST="./resources/edgenew768x768.pgm",MPIPROG="imagenew_timing_intervals" -l $NODES_NUM imagenew.pbs

# check correctness in production
make clean_build
make imagenew
echo "checking correctness in production..."
echo -e $TIMING_HEADER > $DATA_FOLDER/time_results_prod.tsv
PROC_LIST="16"
qsub $RESERVATION -v PROC_LIST="${PROC_LIST}",RESOURCES_LIST="${RESOURCES_LIST}",MPIPROG="imagenew" -l $NODES_NUM imagenew.pbs

echo "wait for the experiments to finish, then run \$python data_analyzer.py to create the appropriate graphs"
