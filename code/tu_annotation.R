###########################################################################
###########################################################################
###
### TRANSCRIPTIONAL UNIT DETECTION AND ANNOTATION
###
###########################################################################
###########################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
packages <- c("tidyverse", "here", "ggthemes", "data.table", 
              "ggExtra", "Rsamtools", "GenomicAlignments", "seqTools", 
              "Rsubread", "ape", "DT", "ggpubr", "ggridges", "ggsci", 
              "GenomicRanges","IRanges")
invisible(lapply(packages, require, character.only = TRUE))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................publication_theme_white
theme_Publication_white <- function(base_size=10) {
  (theme_foundation(base_size=base_size, base_family="Helvetica")
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(3), hjust = 0.5),
           text = element_text(color = "black"),
           panel.background = element_rect(colour = NA, fill = "white"),
           plot.background = element_rect(colour = NA, fill = "white"),
           panel.border = element_rect(colour = "black"),
           axis.title = element_text(face = "bold",size = rel(1.5)),
           axis.title.y = element_text(angle=90,vjust =2, size = rel(1.2)),
           axis.title.x = element_text(vjust = -0.2, size = rel(1.2)),
           axis.text = element_text(size = rel(1.2)), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="grey80"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(color = NA, fill = "white"),
           legend.position = "bottom",
           legend.background= element_rect(color = NA, fill = "white"),
           legend.direction = "horizontal",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing = unit(0, "cm"),
           legend.text = element_text(color = "black", size = rel(1.2)),
           legend.title = element_text(face="italic", color = "black"),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="grey90",fill="grey70"),
           strip.text = element_text(face="bold")
   ))
}

#...................................merge overlapping "regions"
mergeOverlapping <- function(input, minfrac=0) {
  data <- makeGRangesFromDataFrame(input)
  hits <- findOverlaps(data)
  x <- data[queryHits(hits)]
  y <- data[subjectHits(hits)]
  relative_overlap <- width(pintersect(x, y)) / pmin(width(x), width(y))
  hits <- hits[relative_overlap >= minfrac]
  gr1 <- mergeConnectedRanges(data, hits)
  
  list_of_id <- as.data.frame(gr1$revmap) %>%
    dplyr::select(group, value) %>%
    rowwise() %>%
    mutate(corresponding_gene = input$gene[value]) %>%
    distinct(corresponding_gene, group) %>%
    group_by(group) %>%
    arrange(group) %>%
    mutate(genes_in_operon = paste(corresponding_gene, collapse = ","),
           size_operon = 1 + str_count(genes_in_operon, ",")) %>%
    distinct(group, genes_in_operon, size_operon)
  
  new_operon_table <- gr1 %>%
    as.data.frame() %>%
    mutate(operon_id = 1:n()) %>%
    dplyr::left_join(list_of_id, by = c("operon_id" = "group")) %>%
    dplyr::rename(start_operon = start,
                  end_operon = end,
                  width_operon = width,
                  strand_operon = strand) %>%
    dplyr::select(-revmap)
  
  return(new_operon_table)
}

mergeOverlapping_subsets <- function(input, minfrac=0) {
  
  data <- makeGRangesFromDataFrame(input)
  hits <- findOverlaps(data)
  x <- data[queryHits(hits)]
  y <- data[subjectHits(hits)]
  relative_overlap <- width(pintersect(x, y)) / pmin(width(x), width(y))
  hits <- hits[relative_overlap >= minfrac]
  gr1 <- mergeConnectedRanges(data, hits)
  
  # > create bed.like object
  new_operon_table <- gr1 %>%
    as.data.frame() %>%
    mutate(operon_id = 1:n()) %>%
    dplyr::rename(start_operon = start,
                  end_operon = end,
                  width_operon = width,
                  strand_operon = strand) %>%
    dplyr::select(-revmap) %>%
    dplyr::filter(width_operon > 100)
  
  return(new_operon_table)
}

#...................................filter id table for annotation of CDS 
filter_id_table <- function(input_table, amount_of_reads, blacklist,range_5utr_downstream, range_5utr_upstream, range_3utr_downstream, range_3utr_upstream, identity_cutoff){
  
  #.................filter id table for CDS mapping or tRNA mapping reads only! (filter out blacklist regions, eg. rDNA loci)
  working_table <- input_table %>%
    ungroup() %>%
    dplyr::filter(mapped_type == "CDS" | mapped_type == "tRNA", identity > identity_cutoff) %>%
    mutate(length_mapped = abs(start - end)) %>%
    group_by(minion_read_name) %>%
    dplyr::filter(sum((start:end %not in% blacklist) == F) == 0 & length_mapped < 2 * aligned_reads)
  
  #.................filter out all 5´UTR regions that are shorter than -20 long (IN the gene)
  utr5_table <- working_table %>%
    group_by(gene) %>%
    dplyr::mutate(real_utr5 = ifelse(strand == "+" & start < (start_gene + 20), T,
                                     ifelse(strand == "-" & end > (end_gene - 20), T, F))) %>% 
    dplyr::filter(real_utr5 == T) %>%
    mutate(median_utr5 = ifelse(strand == "+", median(start), median(end))) %>%
    dplyr::select(median_utr5, gene) %>%
    distinct()
  
  #.................filter out all 3´UTR regions that are shorter than -20 long (IN the gene)
  utr3_table <- working_table %>%
    group_by(gene) %>%
    dplyr::mutate(real_utr3 = ifelse(strand == "+" & end > (end_gene - 20), T,
                                     ifelse(strand == "-" & start < (start_gene + 20), T, F))) %>% 
    dplyr::filter(real_utr3 == T) %>%
    mutate(median_utr3 = ifelse(strand == "+", median(end), median(start))) %>%
    dplyr::select(median_utr3, gene) %>%
    distinct()
  
  #.................add TSS and TTS information to table with all other information, check sequencing depth after adding TTS and TSS informations
  full_id_table_filtered <- working_table %>%
    left_join(utr5_table, by = "gene") %>%
    left_join(utr3_table, by = "gene") %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(n_mapped_after_filtering = n()) %>%
    ungroup() %>%
    dplyr::filter(n_mapped_after_filtering >= amount_of_reads) %>%
    dplyr::mutate(length_gene = abs(start_gene - end_gene),
                  true_utr3 = ifelse(strand == "+" & end < (median_utr3 + range_3utr_downstream) & end > (median_utr3 - range_3utr_upstream), T, 
                                     ifelse(strand == "-" & start < (median_utr3 + range_3utr_downstream) & start > (median_utr3 - range_3utr_upstream), T, F)),
                  true_utr5 = ifelse(strand == "+" & start < (median_utr5 + range_5utr_upstream) & start > (median_utr5 - range_5utr_downstream), T, 
                                     ifelse(strand == "-" & end < (median_utr5 + range_5utr_downstream) & end > (median_utr5 -range_5utr_upstream), T, F))) %>%
    dplyr::filter((true_utr3 == T & aligned_reads > 0.5 * length_gene) | true_utr5 == T & aligned_reads > 0.5 * length_gene)
  return(full_id_table_filtered)
}

#...................................function to collapse reads that were filtered before
collapse_ids <- function(filtered_id_table){
  
  # > plus strand
  In.df_plus <- filtered_id_table %>%
    ungroup %>%
    dplyr::filter(strand == "+") 
  
  # > minus strand
  In.df_minus <- filtered_id_table %>%
    ungroup %>%
    dplyr::filter(strand == "-") 
  
  # > merge overlapping data sets and combine them in table
  merged_table_plus <- mergeOverlapping(In.df_plus)
  merged_table_minus <- mergeOverlapping(In.df_minus) 
  merged_table_all <- rbind(as.data.table(merged_table_plus), as.data.table(merged_table_minus))
  
  return(merged_table_all)
}

#...................................read in single-read tables, filter by strand, 
collapse_ids_subsets <- function(filtered_id_table){
  
  In.df_plus <- filtered_id_table %>%
    ungroup %>%
    dplyr::filter(strand == "+") 
  In.df_minus <- filtered_id_table %>%
    ungroup %>%
    dplyr::filter(strand == "-") 
  
  if(length(In.df_plus$start) > 0){
    merged_table_plus <- mergeOverlapping_subsets(In.df_plus)
  }
  if(length(In.df_minus$start) > 0){
    merged_table_minus <- mergeOverlapping_subsets(In.df_minus) 
  }
  if(length(In.df_plus$start) > 0 & length(In.df_minus$start) > 0){
    merged_table_all <- rbindlist(list(as.data.table(merged_table_plus), as.data.table(merged_table_minus)))
    return(merged_table_all)
  }else if(length(In.df_plus$start) == 0 & length(In.df_minus$start) > 0){
    merged_table_all <- as.data.table(merged_table_minus)
    return(merged_table_all)
  }else if(length(In.df_plus$start) > 0 & length(In.df_minus$start) == 0){
    merged_table_all <- as.data.table(merged_table_plus)
    return(merged_table_all)
  }
}








#...................................negative function for %in%
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

#...................................Merge ranges that are "connected" (directly or indirectly) via a hit (or several hits) in 'hits'.
mergeConnectedRanges <- function(x, hits){
  stopifnot(is(x, "GenomicRanges"))
  stopifnot(is(hits, "Hits"))
  stopifnot(queryLength(hits) == subjectLength(hits))
  stopifnot(queryLength(hits) == length(x))
  clusters <- extractClustersFromSelfHits(hits)
  ans <- range(extractList(x, clusters))
  if (any(elementNROWS(ans) != 1L))
    stop(wmsg("some connected ranges are not on the same ",
              "chromosome and strand, and thus cannot be ",
              "merged"))
  ans <- unlist(ans)
  mcols(ans)$revmap <- clusters
  ans
}

#...................................extract clusters from Hits object
#...from https://support.bioconductor.org/p/68021/ 
extractClustersFromSelfHits <- function(hits){
  stopifnot(is(hits, "Hits"))
  stopifnot(queryLength(hits) == subjectLength(hits))
  hits <- union(hits, t(hits))
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  cid <- seq_len(queryLength(hits))  # cluster ids
  while (TRUE) {
    h <- Hits(qh, cid[sh],
              queryLength(hits), subjectLength(hits))
    cid2 <- pmin(cid, selectHits(h, "first"))
    if (identical(cid2, cid))
      break
    cid <- cid2
  }
  unname(splitAsList(seq_len(queryLength(hits)), cid))
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load gff genome annotation and blacklist area around rDNA loci

#.................................escherichia coli
gff_table_ecoli <- read.gff(here("data/genome_data/ecoli.gff")) %>%
  as_tibble() %>%
  dplyr::filter(type %in% "rRNA") %>%
  mutate(start = start - 500,
         end = end + 500) %>%
  dplyr::select(start, end)

blacklist_rRNA_ecoli <- list()
for (i in 1:length(gff_table_ecoli$start)){
  blacklist_rRNA_ecoli <- unlist(c(blacklist_rRNA_ecoli, c(gff_table_ecoli$start[i]:gff_table_ecoli$end[i])))
}

#.................................pyrococcus
gff_table_pfu <- read.gff(here("data/genome_data/pfu.gff")) %>%
  as_tibble() %>%
  dplyr::filter(type %in% "rRNA") %>%
  mutate(start = start - 500,
         end = end + 500) %>%
  dplyr::select(start, end)

blacklist_rRNA <- list()
for (i in 1:length(gff_table_pfu$start)){
  blacklist_rRNA <- unlist(c(blacklist_rRNA, c(gff_table_pfu$start[i]:gff_table_pfu$end[i])))
}

#.................................haloferax
gff_table_hvo <- read.gff(here("data/genome_data/hvo.gff")) %>%
  as_tibble() %>%
  dplyr::filter(type %in% "rRNA") %>%
  mutate(start = start - 500,
         end = end + 500) %>%
  dplyr::select(start, end)

blacklist_rRNA_hvo <- list()
for (i in 1:length(gff_table_hvo$start)){
  blacklist_rRNA_hvo <- unlist(c(blacklist_rRNA_hvo, c(gff_table_hvo$start[i]:gff_table_hvo$end[i])))
}

#....................................ESCHERICHIA
#.................................load id data
load(here("data/tidy_data/ecoli_tex_id_table"))

#.................................use function
ecoli_filtered_ids <- filter_id_table(full_id_table, 
                                    blacklist = blacklist_rRNA_ecoli,
                                    amount_of_reads = 0, 
                                    range_5utr_upstream = 20, 
                                    range_5utr_downstream = 1, 
                                    range_3utr_downstream = 20, 
                                    range_3utr_upstream = 1,
                                    identity_cutoff = 80)

#.................................collapse ids
ecoli_collapsed_ids <- collapse_ids(ecoli_filtered_ids)

#.................................write to file
write_tsv(ecoli_filtered_ids, here("data/operon_data/ecoli_tex_reads_for_operons.tsv"))
write_tsv(ecoli_collapsed_ids,  here("data/operon_data/ecoli_tex_operons.tsv"))

#.................................write to bed file for visualization in IGV
write.table(ecoli_collapsed_ids, file=here("data/operon_data/ecoli_collapsed_reads_above80.bed"), quote=F, sep="\t", row.names=F, col.names=F)



#....................................HALOFERAX
#.................................load id data
load(here("data/tidy_data/hvo_tex_id_table"))

#.................................use function
hvo_filtered_ids <- filter_id_table(full_id_table, 
                                    blacklist = blacklist_rRNA_hvo,
                                    amount_of_reads = 0, 
                                    range_5utr_upstream = 20, 
                                    range_5utr_downstream = 1, 
                                    range_3utr_downstream = 20, 
                                    range_3utr_upstream = 1,
                                    identity_cutoff = 80)

#.................................collapse ids
hvo_collapsed_ids <- collapse_ids(hvo_filtered_ids)

#.................................write to file
write_tsv(hvo_filtered_ids, here("data/operon_data/hvo_tex_reads_for_operons.tsv"))
write_tsv(hvo_collapsed_ids,  here("data/operon_data/hvo_tex_operons.tsv"))

#.................................write to bed file for visualization in IGV
write.table(hvo_collapsed_ids, file=here("data/operon_data/hvo_collapsed_reads_above80.bed"), quote=F, sep="\t", row.names=F, col.names=F)

#....................................PYROCOCCUS

#.................................load id data
load(here("data/tidy_data/pfu_tex_id_table"))

#.................................use function
pfu_filtered_ids <- filter_id_table(full_id_table, 
                                    blacklist = blacklist_rRNA,
                                    amount_of_reads = 0, 
                                    range_5utr_upstream = 20, 
                                    range_5utr_downstream = 1, 
                                    range_3utr_downstream = 20, 
                                    range_3utr_upstream = 1,
                                    identity_cutoff = 80)

#.................................collapse ids
pfu_collapsed_ids <- collapse_ids(pfu_filtered_ids)

#.................................write to file
write_tsv(pfu_filtered_ids, here("data/operon_data/pfu_tex_reads_for_operons.tsv"))
write_tsv(pfu_collapsed_ids,  here("data/operon_data/pfu_tex_operons.tsv"))

#.................................write to bed file for visualization in IGV
write.table(pfu_collapsed_ids, file=here("data/operon_data/pfu_collapsed_reads_above80.bed"), quote=F, sep="\t", row.names=F, col.names=F)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SUB-TU DETECTION
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

#...................................functions
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



