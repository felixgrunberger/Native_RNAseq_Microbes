###########################################################################
###########################################################################
###
### TRANSCRIPTIONAL SUB UNIT DETECTION AND ANNOTATION
###
###########################################################################
###########################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES AND PLOTTING FUNCTION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(here)
source(here("Rscripts/load_libraries.R"))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#.................................split
splitAt <- function(x, pos) {
  unname(split(x, findInterval(x, pos)))
}

#.................................split full operons based on trnascript coverage on 3´end
divide_operons_to_singles <- function(filtering_table, collapsed_table, wanted_strand){
  
  if (wanted_strand == "+"){
    coverage_table <- pfu_coverage_forward
  }else if(wanted_strand == "-"){
    coverage_table <- pfu_coverage_reverse
  }
  stepsize <- 20 # upstream and downstream reads of TSS || TTS
  what <- filtering_table %>%
    dplyr::filter(strand == wanted_strand) %>%
    dplyr::select(start,end, width, strand, identity, gene, median_utr5, median_utr3) %>%
    mutate(tss = as.integer(median_utr5), tts = as.integer(median_utr3)) %>%
    dplyr::filter(!is.na(tss), !is.na(tts)) %>%
    rowwise() %>%
    mutate(coverage_tss_downstream = ifelse(strand == "+", mean(coverage_table$depth[(tss):(tss+stepsize)]),
                                            mean(coverage_table$depth[(tss-stepsize):(tss)])),
           coverage_tss_upstream   = ifelse(strand == "+", mean(coverage_table$depth[(tss-stepsize):(tss)]),
                                            mean(coverage_table$depth[(tss):(tss+stepsize)])),
           coverage_tts_downstream = ifelse(strand == "+", mean(coverage_table$depth[(tts-stepsize):(tts)]),
                                            mean(coverage_table$depth[(tts):(tts+stepsize)])),
           coverage_tts_upstream   = ifelse(strand == "+", mean(coverage_table$depth[(tts):(tts+stepsize)]),
                                            mean(coverage_table$depth[(tts-stepsize):(tts)])),
           drop_factor_tss = coverage_tss_downstream/coverage_tss_upstream,
           drop_factor_tts = coverage_tts_downstream/coverage_tts_upstream) %>%
    distinct(tss, .keep_all = TRUE) 
  
  searched_interval <- NA
  
  collapsed_table <- collapsed_table %>%
    dplyr::filter(strand_operon == wanted_strand)
  
  for (i in seq_along(collapsed_table$seqnames)){
    searched_interval[i] <- list(collapsed_table$start_operon[i]:collapsed_table$end_operon[i])
  }
  
  list_of_tts <- sort(what$tts[what$drop_factor_tts > 1.5])
  list_of_tss <- what$tss[what$drop_factor_tss > 1.5]
  list_of_genes <- what$gene[what$drop_factor_tts > 1.5]
  
  final_list <- NA
  for (i in seq_along(collapsed_table$seqnames)){
    final_list <- c(final_list,list(splitAt(unlist(searched_interval[i]),list_of_tts+1)))
  }
  
  start_positions <- NA
  end_positions   <- NA
  
  for (i in seq_along(final_list)){
    for(j in seq_along(final_list[[i]])){
      start_positions <- c(start_positions,min(final_list[[i]][[j]]))
      end_positions   <- c(end_positions, max(final_list[[i]][[j]]))
    }
  }
  
  split_table <- matrix(nrow = length(start_positions), ncol = 2) %>%
    as_tibble() %>%
    mutate(start = start_positions, 
           end = end_positions,
           width = abs(end - start),
           chr = "CP023154",
           strand = wanted_strand) %>%
    dplyr::filter(!is.na(start), !is.na(end))
  
  collapsed_table_all_sub <- collapse_ids_subsets(split_table) 
  return(collapsed_table_all_sub)
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SUB-TU DETECTION (ONLY SHOWN FOR PFU)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#....................................PYROCOCCUS
#...................................load read coverage
#.................................forward
pfu_coverage_forward <- fread(here("data/coverage_data/pfu_tex_plus_depth.txt.gz")) %>%
  dplyr::rename(position = 2, depth = 3) %>%
  mutate(depth = depth + 1) # give every position at least one read!

#.................................reverse
pfu_coverage_reverse <- fread(here("data/coverage_data/pfu_tex_minus_depth.txt.gz")) %>%
  dplyr::rename(position = 2, depth = 3) %>%
  mutate(depth = depth + 1) # give every position at least one read!

#.................................split operons based on 3´end coverage
pfu_sub_operons_plus  <- divide_operons_to_singles(pfu_filtered_ids, pfu_collapsed_ids, wanted_strand = "+")
pfu_sub_operons_minus <- divide_operons_to_singles(pfu_filtered_ids, pfu_collapsed_ids, wanted_strand = "-")

pfu_sub_operons <- rbindlist(list(pfu_sub_operons_plus, pfu_sub_operons_minus)) %>%
  distinct(start_operon, .keep_all = TRUE) %>%
  arrange(start_operon)

#.................................annotate with gene names using granges objects
#...............................gff
gff_cds <- read.gff(here("data/genome_data/pfu.gff")) %>%
  as_tibble() %>%
  dplyr::filter(type %in% "CDS") %>%
  mutate(gene = str_split_fixed(str_split_fixed(attributes, "ID=",2)[,2],".p01;",2)[,1],
         gene_number = 1:n()) %>%
  dplyr::select(start, end, strand, gene, gene_number)

#...............................Make the Granges object
regions <- makeGRangesFromDataFrame(pfu_sub_operons %>% mutate(start = start_operon, end = end_operon, strand = strand_operon))

#...............................Make new metadata column called "feature"
mcols(regions)$feature <- ""

#...............................Make Granges object for the feature of interest and annotate with gene name
genes <- makeGRangesFromDataFrame(gff_cds %>% mutate(seqnames = "CP023154"))
mcols(genes)$feature <- gff_cds$gene

#...............................Find overlaps and assign feature to regions
hits <- findOverlaps(query = regions, subject = genes, ignore.strand = FALSE)
mcols(hits)$feature <- gff_cds$gene[subjectHits(hits)]
mcols(regions[queryHits(hits)])$feature <- hits@elementMetadata$feature

#...............................combine all
pfu_sub_operons_annotated <- as.data.table(hits) %>%
  left_join(as.data.table(regions) %>% mutate(queryHits = 1:n()), by = "queryHits") %>%
  group_by(start) %>%
  arrange(start) %>%
  mutate(genes_in_operon = paste(feature.x, collapse = ","),
         size_operon = 1 + str_count(genes_in_operon, ",")) %>%
  rowid_to_column("operon_id") %>%
  dplyr::rename(start_operon = start, 
                end_operon = end, 
                strand_operon = strand,
                width_operon = width) %>%
  dplyr::select(seqnames, start_operon, end_operon, width_operon, strand_operon, operon_id, genes_in_operon,size_operon) %>%
  distinct(start_operon, .keep_all = TRUE)

#...............................write to table
write_tsv(pfu_sub_operons_annotated, here("data/operon_data/pfu_tex_suboperons.tsv"))

