# Init libraries
import pandas as pd
import cooler

# Global variables
in_dir = "/mnt/c/Users/Marti/Documents/Projects/HK_Interactions/Data/HiC/4DN_portal/average/01_bedpe"
out_dir = "/mnt/c/Users/Marti/Documents/Projects/HK_Interactions/Data/HiC/4DN_portal/average/02_bins"
raw_cool_dir ="/mnt/c/Users/Marti/Documents/Projects/HK_Interactions/Data/HiC/4DN_portal/raw/"
date = "240106"

# In previous test I verify that we can create cool files using the cooler function load.
# For this process we need to provide the bins information and the bedpe files with the mean values. 
# The function need to be run in the command line,
# but I will setup the neccesary files in this notebook.

#define the resolutions to be used
resolutions = ['2000', '5000','10000','25000','50000','100000','250000','500000','1000000','2500000','5000000']

#test resolutions
# resolutions = ['100000','250000']
# resolutions = ['5000000']

#for each resolution, we need to save the bin information and the average matrix
for res in resolutions:
   
    #### Save the bins information ####

    # Load the mcool file for the resolution. We can use any cell type, since the bin information is the same for all of them.
    cool_data = cooler.Cooler(raw_cool_dir + "4DNFII4IPDTR.mcool::resolutions/" + str(res))
    #get the bins table
    bins = cool_data.bins()[:]
    #remove the columns that are not needed
    bins = bins[['chrom', 'start', 'end']]
    #save the bins table as a bed file
    bins.to_csv(f'{out_dir}/bins_{res}.bed', sep='\t', index=False, header=False)

     #### Remove duplicated rows in the average data if exist ####
    # Load the average matrix for the resolution.
    average_matrix = pd.read_csv(f'{in_dir}/average_matrix_{res}.bed', sep='\t', header=None)

    #check for duplicated rows
    # Get the index of duplicated rows
    dup_idx = average_matrix[average_matrix.duplicated].index   
    # if dup_idx is not empty, drop the duplicated rows
    if dup_idx.empty == False:
        average_matrix.drop(dup_idx, inplace=True)
        #save the new average matrix
        average_matrix.to_csv(f'{in_dir}/average_matrix_{res}.bed', sep='\t', index=False, header=False)
    