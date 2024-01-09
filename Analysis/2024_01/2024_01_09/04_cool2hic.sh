#!/bin/bash

# Set the path to the input .mcool file
input_mcool="/mnt/c/Users/Marti/Documents/Projects/HK_Interactions/Data/HiC/4DN_portal/average/04_mcool/average.mcool"

# Set the path to the bedpe file
input_bedpe="/mnt/c/Users/Marti/Documents/Projects/HK_Interactions/Data/HiC/4DN_portal/average/05_hic/bedpe/"

# Set the path to the output .hic file
output_hic="/mnt/c/Users/Marti/Documents/Projects/HK_Interactions/Data/HiC/4DN_portal/average/05_hic//average_100k.hic"

# Set the path to the chrom.sizes file
chrom_sizes="/mnt/c/Users/Marti/Documents/Projects/HK_Interactions/Data/HiC/4DN_portal/average/05_hic/4DN.hg38.chrom.sizes"

# Set the path to the juicer_tools jar file
juicer_tools_jar="/mnt/c/Users/Marti/Documents/Projects/HK_Interactions/Data/HiC/4DN_portal/average/05_hic/juicer_tools_1.22.01.jar"

#set resolution
res=100000

# Bedpe file
output_bedpe=${input_bedpe}average_matrix_${res}.bed
echo $output_bedpe 

echo "Converting to short format"
# Convert the ginteractions file to short format with score using awk
awk -F "\t" '{print 0, $1, $2, 0, 0, $4, $5, 1, $7}' ${output_bedpe} > ${output_bedpe}.short

# Sort the short format with score file
echo "Sorting"
sort -k2,2d -k6,6d ${output_bedpe}.short > ${output_bedpe}.short.sorted

# Convert the short format with score file to .hic using juicer pre
echo "Converting to .hic"
java -Xms512m -Xmx2048m -jar $juicer_tools_jar pre -r 5000,10000,25000,50000,100000,250000,500000,1000000 ${output_bedpe}.short.sorted $output_hic $chrom_sizes