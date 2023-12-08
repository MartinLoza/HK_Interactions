#!/bin/sh

cd ~/Documents/Projects/HK_Interactions/Analysis/2023_10/2023_10_17/

# computeMatrix scale-regions \
#  -S /Volumes/MARTIN_LOZA/Projects/HK_Interactions/Data/ChIP_TF/ENCODE/CTCF_K562.bigWig \
#     /Volumes/MARTIN_LOZA/Projects/HK_Interactions/Data/ChIP_TF/ENCODE/CTCFL_K562.bigWig \
#  -R tmp_results/hub_ppi_top_1000_3k.bed \
#  --skipZeros \
#  -o Results_deeptools/CTCF_CTCFL_matrix_231018.mat.gz \
#  --beforeRegionStartLength 100000 \
#  --afterRegionStartLength 100000

plotHeatmap \
 -m Results_deeptools/CTCF_CTCFL_matrix_231018.mat.gz \
 -out Results_deeptools/CTCF_CTCFL_HubPPI_Top1k_3kb.png \
 --colorMap Blues \
 --zMin 0 --zMax 3 
#for scaling the heatmap, the plot in the top is not scaled
 #--zMin -3 --zMax 3  
 

