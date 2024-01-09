#!/bin/bash

# Main directory
INPUT_DIR=/mnt/c/Users/Marti/Documents/Projects/HK_Interactions/Data/HiC/4DN_portal/average
# Directory of bins
BIN_DIR=${INPUT_DIR}/02_bins
# Directory of bedpe files
BEDPE_DIR=${INPUT_DIR}/01_bedpe
# Output directory
COOL_DIR=${INPUT_DIR}/03_cool
# MCOOL directory
MCOOL_DIR=${INPUT_DIR}/04_mcool

# Define resolutions to use
resolutions=('2000' '5000' '10000' '25000' '50000' '100000' '250000' '500000' '1000000')
# resolutions=('5000' '10000' '25000' '50000' '100000' '250000' '500000' '1000000' '2500000' '5000000')
#test resolutions
# resolutions=('5000000' '10000000')

# For each resolution
# for resolution in "${resolutions[@]}"
# do
#     #print test text
#     # echo ${resolution}
#     # Cooler load command
#     cooler load -f bg2 \
#      $BIN_DIR/bins_${resolution}.bed \
#      $BEDPE_DIR/average_matrix_${resolution}.bed ${COOL_DIR}/${resolution}.cool \
#      --count-as-float
# done

# Create single mcool file for all resolutions

# Concatenate "cooler/" with each resolution and ".cool" at the end
resolutions_with_prefix=("${resolutions[@]/#/$COOL_DIR/}")
resolutions_with_suffix=("${resolutions_with_prefix[@]/%/.cool}")

# Print all resolutions with prefix and suffix
echo "${resolutions_with_suffix[@]}"

# Use hicConvertFormat to convert all resolutions to mcool
hicConvertFormat \
    --matrices ${resolutions_with_suffix[@]} \
    --outFileName $MCOOL_DIR/average_3.mcool \
    --inputFormat cool \
    --outputFormat mcool \
    # --resolutions ${resolutions[@]} \
    # --load_raw_values