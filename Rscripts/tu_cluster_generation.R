###########################################################################
###########################################################################
###
### TRANSCRIPTIONAL UNIT CLUSTER DETECTION AND ANNOTATION
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

#...................................read in single-read tables, filter by strand 
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
