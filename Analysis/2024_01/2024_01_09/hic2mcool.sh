# Input data folder
INPUT_DIR=/mnt/c/Users/Marti/Documents/Projects/HK_Interactions/Data/HiC/ENCODE/raw
# Output data folder
OUTPUT_DIR=/mnt/c/Users/Marti/Documents/Projects/HK_Interactions/Data/HiC/ENCODE/raw

# Load the metadata file
METADATA=$INPUT_DIR/metadata.tsv

# Our columns of interest are the 1st and 11th, File accession and Biosample term name
# test. print the first 5 lines of the 1st and 11th columns 
# cut -f 1,11 $METADATA | head -n 5

# For each file in input directory, convert to mcool format
# for FILE in $INPUT_DIR/*.hic
# do
#     .
#     Get the file accession number
#     FILE_ACCESSION=$(basename $FILE .hic)

#     # This commented part is for saving the cool files using their biosample name.
#     #However, as some names are too long, I will use their accession number instead
#     # # Get the biosample term name
#     # BIOSAMPLE_TERM_NAME=$(grep $FILE_ACCESSION $METADATA | cut -f 11)
    
#     # Convert to cool format (hicexplorer cant convert hic to mcool directly)
#     Convert to mcool format
#     hicConvertFormat -m $FILE -o $OUTPUT_DIR/$FILE_ACCESSION.cool \
#     --inputFormat hic \
#     --outputFormat cool \
#     --resolutions 1000000 5000000
# done



#test for using hicconvert to convert .hic to .mcool
hicConvertFormat -m $INPUT_DIR/ENCFF247KEK.hic -o $OUTPUT_DIR/test.cool \
--inputFormat hic \
--outputFormat cool \
--resolutions 10000 \
--load_raw_values