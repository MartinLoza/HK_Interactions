# ---
# title: "Housekeeping P-P interactions"
# author: "Martin Loza"
# date: "2023/09/11
# ---

# I would like to explore the existence of housekeeping promoter-promoter interactions. In other words, promoter-promoter interactions found in a large number of cell types.

# Setup
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(GenomicRanges)


## Global variables
in_dir  <- "/Volumes/MARTIN_LOZA/Projects/Enhancer_Grammar/Data/"
out_dir <- "~/Documents/Projects/Enhancer_Grammar/Analysis/Weekly_analyses/2023_09/Results"

date = 230911
red = "#CC6677"
blue = "#6699CC"
yellow = "#DDCC77"
dpi = 300
text_size = 14
fig_size = 4
point_size = 1
alpha = 0.3

#Functions
source(file = "~/Documents/Projects/Enhancer_Grammar/Analysis/Functions.R")
source(file = "~/Documents/Projects/Enhancer_Grammar/Analysis/Functions_Visualizations.r")


print("Load data")

# load cres with their predicted targets
cres <- readRDS(file = "/Volumes/MARTIN_LOZA/Projects/Enhancer_Grammar/Data/v4/annotated/all_regions_annotated_230317.rds")
cres %>% filter(!duplicated(region_id)) %>% pull(n_celltypes) %>% table()

# load the original regions containing the whole set of cell types
regions_all <- readRDS(file = "/Volumes/MARTIN_LOZA/Projects/Enhancer_Grammar/Data/v4/raw_enhancers/merged/05_raw_enhancers_selected_50cts_221121.rds")


print("Set up")
# Filter only the cres from selected cell types
regions_all <- regions_all %>% filter(CellType %in% cres$celltypes)

#select important columns 
cres <- cres %>% select(chr, start, end, region_id, center, length, celltypes, n_celltypes, test_n_celltypes, gene, tss, distance, n_genes, nearest_gene, nearest_transcript_id, nearest_transcript_tss, dist_to_nearest_transcript, ann_v1)

# I would like to select interactions only for CP
cres <- cres %>% filter(ann_v1 == "CP")
cres$ann_v1 <- droplevels(cres$ann_v1)
cres %>% filter(!duplicated(region_id)) %>% pull(n_celltypes) %>% table()

print("Analysis")

# Let's rename CREs to interactions. 
interactions_cp <- cres 
rm(cres)
gc()
#add new identifier with region_id and target gene
interactions_cp <- interactions_cp %>% mutate(id = paste0(region_id,"_", gene))
#test 
any(duplicated(interactions_cp$id))

# To remove any bias added by a core-promoter predicted to interact with itself, let's remove CREs annotated that are predicted to interact with a target within 500 bps
#remove cres_gene interactions of CP targetting itself withing 500 bps
rmv_df <- interactions_cp %>% 
  filter(gene == nearest_gene) %>% 
  filter(distance <= 500)

#remove identified self-target from interactions df
interactions_cp <- interactions_cp %>% filter(!id %in% rmv_df$id)

rm(rmv_df)

# I would like count, for each interaction in how many cell types is found
# For each unique cres, we would like to count the number of celltypes that a given target appears

#get unique cres
u_cres <- interactions_cp %>% filter(!duplicated(region_id)) 
#for test
# u_cres <- u_cres %>% filter(test_n_celltypes == 50)
print("Loop")
# n_ct_interaction <- lapply(X = u_cres$region_id, FUN = function(c_region_id){
  #test
  n_ct_interaction <- lapply(X = u_cres$region_id[1:10], FUN = function(c_region_id){
  c_region <- u_cres %>% filter(region_id == c_region_id)
  
  #get overlapping interactions in the original data set
  ovl_idx <- WhichRegionsOverlap(query_enhancer = c_region,
                                 subject_enhancer = regions_all)
    #get overlap regions from the original raw data
  all_overlaps <- regions_all[subjectHits(ovl_idx),]
  #count the number of times a target gene appear, e.g. how many cell types it appears as target
  c_n_cts <- all_overlaps %>% group_by(TargetGene) %>% count() %>%
    mutate(id = paste0(c_region$region_id,"_", TargetGene)) %>% ungroup() %>% select(id, n)
  
  return(c_n_cts)
})

n_ct_interaction <- Reduce(f = rbind, x = n_ct_interaction)

#save tmp results
saveRDS(object = n_ct_interaction, file = paste0(out_dir,"/tmp/n_ct_interactions_",date,".rds"))
#read tmp results
# n_ct_interaction <- readRDS(file = paste0(out_dir,"/tmp/n_ct_interactions_",date,".rds"))








