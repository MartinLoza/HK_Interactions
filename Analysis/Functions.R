# Functions used on this project

## Function name: GetRegionsCenter
## output:    Get the center of a given region. It needs the start and the length of the regions
GetRegionsCenter <- function(.data){
  #if length feature doesn't exist
  if(!"length" %in% colnames(.data)){
    .data <- .data %>% mutate("length" = end - start)
  }
  #get the center
  center = .data$start + floor(.data$length/2)
  .data <- .data %>% mutate("center" = center)
  return(.data)
}

## Function name: RemoveDuplicatedRows
## output:    Remove duplicated rows in a df
RemoveDuplicatedRows <- function(.data){
  is_duplicated <- duplicated(x = .data)
  return(.data[!is_duplicated,])
}

## Function name: GetUniqueRegions
## output:    Df containing only unique regions
GetUniqueRegions <- function(.data = NULL){
  .data <- .data %>% filter(!duplicated(region_id))
  return(.data)
}

# Which raw and merged regions overlaps
WhichOverlapEnhancer <- function(query_enhancer = NULL, subject_enhancer = NULL, min_overlap = 0){
  #get the query and subject genomic ranges objects
  gr_query <- GRanges(seqnames = query_enhancer$chr,
                      ranges = IRanges(start = query_enhancer$start,
                                       end = query_enhancer$end))
  gr_subject <- GRanges(seqnames = subject_enhancer$chr,
                        ranges = IRanges(start = subject_enhancer$start,
                                         end = subject_enhancer$end))
  #Find the overlaps
  overlaps <- findOverlaps(query = gr_query, subject = gr_subject,
                           minoverlap = min_overlap)
  #Get the overlapping region's ids
  overlaps_ids <- unique(subjectHits(overlaps))
  #return the ids
  return(overlaps_ids)
}

#Check enhancers data frame
CheckEnhancers <- function(enhancers = NULL){
  #check length
  cat("\nLength OK: \n", 
      sum(enhancers$length == (enhancers %>% mutate(t = end - start) %>% pull(t))) == nrow(enhancers))
  #check center
  cat("\nCenter OK: \n", 
      sum(enhancers$center == (enhancers %>% mutate(t = (start + floor(length/2))) %>% pull(t))) == nrow(enhancers))
  #check chromosome
  cat("\nChromosome OK: \n", 
      sum(enhancers$chr == factor(x = enhancers$chr, levels = paste0("chr", c(1:22,"X")))) == nrow(enhancers))
  #check distance
  cat("\nDistance OK: \n", 
      sum(enhancers$distance == (enhancers %>% mutate(t = abs(tss - center)) %>% pull(t))) == nrow(enhancers))
  #check region ids
  cat("\nRegion ids OK:\n",
      sum(enhancers %>% 
            mutate(t = paste0(chr,":",start,"-",end)) %>%
            mutate (tt = (t == region_id)) %>% 
            pull(tt)) == nrow(enhancers))
  #check unique regions ids
  cat("\nUnique region ids OK:\n",
      enhancers %>% 
        pull(region_id) %>% 
        unique() %>% length() == nrow(enhancers))
}

SetupRegions <- function(regions = NULL){
  
  # Get the length (sometimes we miss one base)
  regions <- regions %>% mutate(length = end - start) 
  # Get center of region
  regions <- regions %>% mutate(center = start + floor(length/2))
  # Create region id for regions
  regions <- regions %>% mutate(region_id = paste0(chr,":",start,"-",end))
  # Set chromosome information as factor
  chr_levels <- paste0("chr", c(1:22,"X"))
  regions <- regions %>% mutate(chr = factor(x = chr, levels = chr_levels))
  
  return(regions)
}

# Get regions overlapping
WhichRegionsOverlap <- function(chr = NULL, query_enhancer = NULL, subject_enhancer = NULL, min_overlap = 0){
  
  #get the query and subject genomic ranges objects
  if(is.null(chr)){
    gr_query <- GRanges(seqnames = query_enhancer[,"chr"],
                        ranges = IRanges(start = query_enhancer[,"start"],
                                         end = query_enhancer[,"end"]))
    gr_subject <- GRanges(seqnames = subject_enhancer[,"chr"],
                          ranges = IRanges(start = subject_enhancer[,"start"],
                                           end = subject_enhancer[,"end"]))
  }else{
    gr_query <- GRanges(seqnames = chr,
                        ranges = IRanges(start = query_enhancer[,"start"],
                                         end = query_enhancer[,"end"]))
    gr_subject <- GRanges(seqnames = chr,
                          ranges = IRanges(start = subject_enhancer[,"start"],
                                           end = subject_enhancer[,"end"]))
  }
  
  #Find the overlaps
  overlaps <- findOverlaps(query = gr_query, subject = gr_subject,
                           minoverlap = min_overlap, )
  #return the ids
  return(overlaps)
}


#Create region ids using the chr, start and end
GetRegionIds <- function(.data){
  .data <- .data %>% mutate(region_id = paste0(chr,":",start,"-",end))
  return(.data)
}

#Filter matrix (similar to dplry one)
FilterMatrix <- function(m = NULL, col = NULL, wt = NULL){
  idx <- which(m[,col] == wt)
  return(m[idx, , drop = FALSE])
}

#clustering with louvain algorithm
ClusterLouvain <- function(x, k = 10, resolution = 0.5) {
  
  g <- bluster::makeSNNGraph(x, k = k)
  res <- igraph::cluster_louvain(g, resolution = resolution)
  
  memberships <- igraph::membership(res)
  
  return(memberships)
}

#Rename columns
RenameColumn <- function(.data = NULL, old_column_name = NULL, new_column_name = NULL){
  #get column index
  col_idx <- which(colnames(.data) == old_column_name)
  #replace the name
  colnames(.data)[col_idx] <- new_column_name
  return(.data)
}

RemoveColumn <- function(.data = NULL, column = NULL){
  #ut. Columns in colnames
  if(sum(!column %in% colnames(.data)) == length(column)){
    stop(call. = TRUE, "ERROR: Columns not found in the data.")
  }
  idx <- which(colnames(.data) %in% column)
  .data <- .data[,-idx]
  return(.data)
}

# MapGeneLabels <- function(genes_labels = NULL, from_type = NULL, to_type = NULL){
#   mapping <- AnnotationDbi::select(x = org.Hs.eg.db, 
#                                    keys = genes_labels,
#                                    keytype=from_type,
#                                    columns = to_type)
#   colnames(mapping) = c("original_name", "mapped_name")
#   #remove unmapped labels
#   mapping <- mapping %>% filter(!is.na(mapped_name))
#   return(mapping)
# }

GetGenesAnnotations <- function(keys = NULL, keys_type = "SYMBOL", annotations = c("SYMBOL","ENSEMBL","ENTREZID"), organism = "human"){
  if(organism == "human"){
    gene_annotations <- AnnotationDbi::select(x = org.Hs.eg.db, 
                                              keys = keys,
                                              columns = annotations,
                                              keytype = keys_type)
  }
  return(gene_annotations)
}
 

