#!/bin/bash

# Set the path to the input .mcool file
input_mcool="./average.mcool"

# Set the path to the output .hic file
output_hic="./average_50k.hic"

# Set the path to the chrom.sizes file
chrom_sizes="./hg38.chrom.sizes.txt"

# Set the path to the juicer_tools jar file
juicer_tools_jar="./juicer_tools_1.22.01.jar"

#set resolution
res=50000
# Get the resolutions stored in the .mcool file
resolutions=$(h5ls -r $input_mcool | grep -Eo 'resolutions/[0-9]+' | cut -d '/' -f 2 | sort -n | uniq)
echo $resolutions
# highest_res=$(echo $resolutions | tr ' ' '\n' | head -n 1)
# echo "highest resolution: $highest_res"

# Use Cooler to write the .mcool matrix as interactions in bedpe format
# output_bedpe=$(echo $input_mcool | sed "s/.mcool/.${highest_res}.bedpe/")
# echo -e "cooler dump --join -r $highest_res $input_mcool::/resolutions/$highest_res"
# cooler dump --join $input_mcool::/resolutions/$highest_res > $output_bedpe

output_bedpe=$(echo $input_mcool | sed "s/.mcool/.${res}.bedpe/")
echo -e "cooler dump --join -r $res $input_mcool::/resolutions/$res"
cooler dump --join $input_mcool::/resolutions/$res > $output_bedpe

# Convert the ginteractions file to short format with score using awk
awk -F "\t" '{print 0, $1, $2, 0, 0, $4, $5, 1, $7}' ${output_bedpe} > ${output_bedpe}.short

# Sort the short format with score file
sort -k2,2d -k6,6d ${output_bedpe}.short > ${output_bedpe}.short.sorted

# Convert the short format with score file to .hic using juicer pre
java -Xms512m -Xmx2048m -jar $juicer_tools_jar pre -r 50000,100000,250000,500000,1000000 ${output_bedpe}.short.sorted $output_hic $chrom_sizes