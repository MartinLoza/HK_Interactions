#!/bin/sh

cd ~/Documents/Projects/HK_Interactions/Analysis/2023_10/2023_10_17/

# computeMatrix scale-regions \
#  -S /Volumes/MARTIN_LOZA/Projects/HK_Interactions/Data/ChIP_TF/ENCODE/CTCF_K562.bigWig \
#     /Volumes/MARTIN_LOZA/Projects/HK_Interactions/Data/ChIP_TF/ENCODE/CTCFL_K562.bigWig \
#  -R Results_deeptools/bed/hk_cp.bed \
#     Results_deeptools/bed/cts_cp.bed \
#  --skipZeros \
#  -o Results_deeptools/HKvsCTS_matrix_231018.mat.gz 

plotHeatmap \
 -m Results_deeptools/HKvsCTS_matrix_231018.mat.gz \
 -out Results_deeptools/CTCF_CTCFL_CTSvsHK_231018.png \
 --colorMap Blues \
 --zMin 0 --zMax 3  
#for scaling the heatmap, the plot in the top is not scaled
 #--zMin -3 --zMax 3  
 

