#!/bin/bash

# Move to main directory
cd /Users/martin/Documents/Projects/HK_Interactions/Analysis/2023_12/2023_12_26/Results/

# Define resolutions to use
resolutions=('5000' '10000' '25000' '50000' '100000' '250000' '500000' '1000000' '2500000' '5000000')

# For each resolution
for resolution in "${resolutions[@]}"
do
    #print test text
    # echo ${resolution}
    # Cooler load command
    cooler load -f bg2 \
     bins/bins_${resolution}.bed \
     average_bedpe/mean_values_${resolution}.bed cooler/${resolution}.cool \
     --count-as-float
    
# Create single mcool file for all resolutions

# Concatenate "cooler/" with each resolution and ".cool" at the end
resolutions_with_prefix=("${resolutions[@]/#/cooler/}")
resolutions_with_suffix=("${resolutions_with_prefix[@]/%/.cool}")

# Print all resolutions with prefix and suffix
# echo "${resolutions_with_suffix[@]}"

# Use hicConvertFormat to convert all resolutions to mcool
hicConvertFormat \
    --matrices ${resolutions_with_suffix[@]} \
    --outFileName cooler/average.mcool \
    --inputFormat cool \
    --outputFormat mcool \
    # --resolutions ${resolutions[@]} \
