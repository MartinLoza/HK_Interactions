#!/bin/sh

# Define the directory paths
data_dir="/Volumes/MARTIN_LOZA/Projects/HK_Interactions/Data/"
in_dir="/Users/martin/Documents/Projects/HK_Interactions/Analysis/2023_10/2023_10_17/Results_deeptools/bed/"
out_dir="/Users/martin/Documents/Projects/HK_Interactions/Analysis/2023_10/2023_10_24/Results_deeptools/"

# echo $data_dir

cd ${out_dir}

# computeMatrix scale-regions \
#  -S ${data_dir}ChIP_TF/ENCODE/CTCF_K562.bigWig \
#     ${data_dir}ChIP_TF/ENCODE/CTCFL_K562.bigWig \
#     ${data_dir}ChIP_TF/ENCODE/YY1_K562_ENCFF913JYP.bigWig \
#     ${data_dir}ChIP_TF/ENCODE/SP1_K562_ENCFF475BKW.bigWig \
#  -R ${in_dir}hub_hk.bed \
#     ${in_dir}nonhub_hk.bed \
#     ${in_dir}cts_cp.bed \
#  -o ${out_dir}matrices/t_matrix_231019.mat.gz \
#  --skipZeros
 
plotHeatmap \
 -m ${out_dir}matrices/t_matrix_231019.mat.gz \
 -out ${out_dir}plots/TFs_hub_nonHub_cts_231019.png \
 --colorMap Blues \
 --zMin 0 --zMax 3  
#for scaling the heatmap, the plot in the top is not scaled
 #--zMin -3 --zMax 3  
 