#!/bin/sh

#Setup

#Global variables
in_dir='/Volumes/MARTIN_LOZA/Projects/HK_Interactions/Data/HiC/PolyII_dynamics_NG/corrected/';
out_dir='/Users/martin/Documents/Projects/HK_Interactions/Analysis/2023_10/2023_10_02/';
date=230928

# Analysis

#set resolution to use (microC data)
resolution=2k

#set input and output files
file_dir="${in_dir}res${resolution}.h5"
file_out="${out_dir}tracks_${date}.ini"

#echo "${out_dir}tracks_230928.ini"
Echo "Finished"


