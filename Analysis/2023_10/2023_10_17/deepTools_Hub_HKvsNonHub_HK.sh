#!/bin/sh

cd ~/Documents/Projects/HK_Interactions/Analysis/2023_10/2023_10_17/

in_dir="~/Documents/Projects/HK_Interactions/Analysis/2023_10/2023_10_17/Results_deeptools/bed/"

# computeMatrix scale-regions \
#  -S /Volumes/MARTIN_LOZA/Projects/HK_Interactions/Data/ChIP_Histone/ENCODE/H3K27ac_K562.bigWig \
#     /Volumes/MARTIN_LOZA/Projects/HK_Interactions/Data/ChIP_Histone/ENCODE/H3K4me3_K562.bigWig \
#     /Volumes/MARTIN_LOZA/Projects/HK_Interactions/Data/ChIP_Histone/ENCODE/H3K4me1_K562.bigWig \
#     /Volumes/MARTIN_LOZA/Projects/HK_Interactions/Data/ChIP_TF/ENCODE/CTCF_K562.bigWig \
#     /Volumes/MARTIN_LOZA/Projects/HK_Interactions/Data/ChIP_TF/ENCODE/CTCFL_K562.bigWig \
#  -R Results_deeptools/bed/hub_hk.bed \
#     Results_deeptools/bed/nonhub_hk.bed \
#     Results_deeptools/bed/cts_cp.bed \
#  --skipZeros \
#  -o Results_deeptools/Hub_HKvsNonHub_HK_matrix_histone_231018.mat.gz 

plotHeatmap \
 -m Results_deeptools/Hub_HKvsNonHub_HK_matrix_histone_231018.mat.gz \
 -out Results_deeptools/CTCF_histone_Hub_HKvsNonHub_HK_vsCTS_CP_231018.png \
 --colorMap Blues \
 --zMin 0 --zMax 3 
#for scaling the heatmap, the plot in the top is not scaled
 #--zMin -3 --zMax 3  
 

