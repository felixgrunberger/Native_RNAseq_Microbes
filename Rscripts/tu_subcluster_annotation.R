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


#.................................split
splitAt <- function(x, pos) {
  unname(split(x, findInterval(x, pos)))
}

#.................................split full operons based on trnascript coverage on 3´end
divide_operons_to_singles <- function(coverage_forward, coverage_reverse, organism, filtering_table, collapsed_table, wanted_strand){

  if (wanted_strand == "+"){
    coverage_table <- coverage_forward
  }else if(wanted_strand == "-"){
    coverage_table <- coverage_reverse
  }
  stepsize <- 20 # upstream and downstream reads of TSS || TTS
  
  what <- filtering_table %>%
    dplyr::filter(strand == wanted_strand) %>%
    dplyr::select(seqnames, start,end, width, strand, identity, gene, median_utr5, median_utr3) %>%
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
    names(searched_interval)[i] <- collapsed_table$seqnames[i]
  }
  
  list_of_tts   <- sort(what$tts[what$drop_factor_tts > 1.5])
  list_of_tss   <- what$tss[what$drop_factor_tss > 1.5]
  list_of_genes <- what$gene[what$drop_factor_tts > 1.5]
  
  final_list <- NA
  for (i in seq_along(collapsed_table$seqnames)){
    final_list <- c(final_list,list(splitAt(unlist(searched_interval[i]),list_of_tts+1)))
  }
  
  start_positions <- list()
  end_positions   <- list()
  chrom_info      <- NA
  
  if (organism == "hvo"){
    for (i in seq_along(final_list)){
      for(j in seq_along(final_list[[i]])){
        chrom_info      <- c(chrom_info,    as.character(levels(as.factor(substr(names(unlist(final_list[[i]][[j]])), 1,11)))))
        start_positions <- c(start_positions,min(final_list[[i]][[j]]))
        end_positions   <- c(end_positions, max(final_list[[i]][[j]]))
      }
    }
    }else if(organism == "pfu" | organism == "ecoli"){
      for (i in seq_along(final_list)){
        for(j in seq_along(final_list[[i]])){
          chrom_info      <- c(chrom_info,    as.character(levels(as.factor(substr(names(unlist(final_list[[i]][[j]])), 1,8)))))
          start_positions <- c(start_positions,min(final_list[[i]][[j]]))
          end_positions   <- c(end_positions, max(final_list[[i]][[j]]))
        }
      }
    }
  
  split_table <- matrix(nrow = length(start_positions), ncol = 2) %>%
    as_tibble() %>%
    mutate(start = as.numeric(start_positions), 
           end = as.numeric(end_positions),
           width = abs(end - start),
           chr = chrom_info,
           strand = wanted_strand) %>%
    dplyr::filter(!is.na(start), !is.na(end))

  collapsed_table_all_sub <- collapse_ids_subsets(split_table) 
  return(collapsed_table_all_sub)
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
# SUB-TU DETECTION (PYROCOCCUS)
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

#...................................load operon data
pfu_filtered_ids  <- fread(here("data/operon_data/pfu_tex_reads_for_operons.tsv"))
pfu_collapsed_ids <- fread(here("data/operon_data/pfu_tex_operons.tsv"))

#.........................split operons based on 3´end coverage
pfu_sub_operons_plus  <- divide_operons_to_singles(pfu_coverage_forward, pfu_coverage_reverse, "pfu",pfu_filtered_ids, pfu_collapsed_ids, wanted_strand = "+")
pfu_sub_operons_minus <- divide_operons_to_singles(pfu_coverage_forward, pfu_coverage_reverse, "pfu",pfu_filtered_ids, pfu_collapsed_ids, wanted_strand = "-")

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
  dplyr::select(start, end, strand, gene, gene_number, seqid)

#...............................Make the Granges object
regions <- makeGRangesFromDataFrame(pfu_sub_operons %>% mutate(start = start_operon, end = end_operon, strand = strand_operon))

#...............................Make new metadata column called "feature"
mcols(regions)$feature <- ""

#...............................Make Granges object for the feature of interest and annotate with gene name
genes <- makeGRangesFromDataFrame(gff_cds %>% dplyr::rename(seqnames = seqid))
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

#...............................extract relevant information
pfu_sub_operons_annotated_export <- pfu_sub_operons_annotated %>%
  dplyr::rename(chr = seqnames) 

#...............................write to table
writexl::write_xlsx(x = pfu_sub_operons_annotated_export, path = here("tables/tu_tables/tu_pfu.xlsx"))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SUB-TU DETECTION (HALOFERAX)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load read coverage
#.................................forward
hvo_coverage_forward <- fread(here("data/coverage_data/hvo_tex_plus_depth.txt.gz")) %>%
  dplyr::rename(position = 2, depth = 3) %>%
  mutate(depth = depth + 1) # give every position at least one read!

#.................................reverse
hvo_coverage_reverse <- fread(here("data/coverage_data/hvo_tex_minus_depth.txt.gz")) %>%
  dplyr::rename(position = 2, depth = 3) %>%
  mutate(depth = depth + 1) # give every position at least one read!

#...................................load operon data
hvo_filtered_ids  <- fread(here("data/operon_data/hvo_tex_reads_for_operons.tsv"))
hvo_collapsed_ids <- fread(here("data/operon_data/hvo_tex_operons.tsv"))

#.........................split operons based on 3´end coverage
hvo_sub_operons_plus  <- divide_operons_to_singles(hvo_coverage_forward, hvo_coverage_reverse, "hvo",hvo_filtered_ids, hvo_collapsed_ids, wanted_strand = "+")
hvo_sub_operons_minus <- divide_operons_to_singles(hvo_coverage_forward, hvo_coverage_reverse, "hvo",hvo_filtered_ids, hvo_collapsed_ids, wanted_strand = "-")

hvo_sub_operons <- rbindlist(list(hvo_sub_operons_plus, hvo_sub_operons_minus)) %>%
  distinct(start_operon, .keep_all = TRUE) %>%
  arrange(start_operon)

#.................................annotate with gene names using granges objects
#...............................gff
gff_cds <- read.gff(here("data/genome_data/hvo.gff")) %>%
  as_tibble() %>%
  dplyr::filter(type %in% "CDS") %>%
  mutate(gene = str_split_fixed(str_split_fixed(attributes, "ID=",2)[,2],";Parent=",2)[,1],
         gene_number = 1:n()) %>%
  dplyr::select(start, end, strand, gene, gene_number, seqid)

#...............................Make the Granges object
regions <- makeGRangesFromDataFrame(hvo_sub_operons %>% mutate(start = start_operon, end = end_operon, strand = strand_operon))

#...............................Make new metadata column called "feature"
mcols(regions)$feature <- ""

#...............................Make Granges object for the feature of interest and annotate with gene name
genes <- makeGRangesFromDataFrame(gff_cds %>% dplyr::rename(seqnames = seqid))
mcols(genes)$feature <- gff_cds$gene

#...............................Find overlaps and assign feature to regions
hits <- findOverlaps(query = regions, subject = genes, ignore.strand = FALSE)
mcols(hits)$feature <- gff_cds$gene[subjectHits(hits)]
mcols(regions[queryHits(hits)])$feature <- hits@elementMetadata$feature

#...............................combine all
hvo_sub_operons_annotated <- as.data.table(hits) %>%
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

#...............................extract relevant information
hvo_sub_operons_annotated_export <- hvo_sub_operons_annotated %>%
  dplyr::rename(chr = seqnames) 

#...............................write to table
writexl::write_xlsx(x = hvo_sub_operons_annotated_export, path = here("tables/tu_tables/tu_hvo.xlsx"))



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SUB-TU DETECTION (ESCHERICHIA)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load read coverage
#.................................forward
ecoli_coverage_forward <- fread(here("data/coverage_data/ecoli_tex_plus_depth.txt.gz")) %>%
  dplyr::rename(position = 2, depth = 3) %>%
  mutate(depth = depth + 1) # give every position at least one read!

#.................................reverse
ecoli_coverage_reverse <- fread(here("data/coverage_data/ecoli_tex_minus_depth.txt.gz")) %>%
  dplyr::rename(position = 2, depth = 3) %>%
  mutate(depth = depth + 1) # give every position at least one read!

#...................................load operon data
ecoli_filtered_ids  <- fread(here("data/operon_data/ecoli_tex_reads_for_operons.tsv"))
ecoli_collapsed_ids <- fread(here("data/operon_data/ecoli_tex_operons.tsv"))

#.........................split operons based on 3´end coverage
ecoli_sub_operons_plus  <- divide_operons_to_singles(ecoli_coverage_forward, ecoli_coverage_reverse, "ecoli", ecoli_filtered_ids, ecoli_collapsed_ids, wanted_strand = "+")
ecoli_sub_operons_minus <- divide_operons_to_singles(ecoli_coverage_forward, ecoli_coverage_reverse, "ecoli", ecoli_filtered_ids, ecoli_collapsed_ids, wanted_strand = "-")

ecoli_sub_operons <- rbindlist(list(ecoli_sub_operons_plus, ecoli_sub_operons_minus)) %>%
  distinct(start_operon, .keep_all = TRUE) %>%
  arrange(start_operon)

#.................................annotate with gene names using granges objects
#...............................gff
gff_cds <- read.gff(here("data/genome_data/ecoli.gff")) %>%
  as_tibble() %>%
  dplyr::filter(type %in% "CDS") %>%
  mutate(gene = str_split_fixed(str_split_fixed(attributes, "ID=",2)[,2],";Parent=",2)[,1],
         gene_number = 1:n()) %>%
  dplyr::select(start, end, strand, gene, gene_number, seqid)

#...............................Make the Granges object
regions <- makeGRangesFromDataFrame(ecoli_sub_operons %>% mutate(start = start_operon, end = end_operon, strand = strand_operon))

#...............................Make new metadata column called "feature"
mcols(regions)$feature <- ""

#...............................Make Granges object for the feature of interest and annotate with gene name
genes <- makeGRangesFromDataFrame(gff_cds %>% dplyr::rename(seqnames = seqid))
mcols(genes)$feature <- gff_cds$gene

#...............................Find overlaps and assign feature to regions
hits <- findOverlaps(query = regions, subject = genes, ignore.strand = FALSE)
mcols(hits)$feature <- gff_cds$gene[subjectHits(hits)]
mcols(regions[queryHits(hits)])$feature <- hits@elementMetadata$feature

#...............................combine all
ecoli_sub_operons_annotated <- as.data.table(hits) %>%
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

#...............................extract relevant information
ecoli_sub_operons_annotated_export <- ecoli_sub_operons_annotated %>%
  dplyr::rename(chr = seqnames) 

#...............................write to table
writexl::write_xlsx(x = ecoli_sub_operons_annotated_export, path = here("tables/tu_tables/tu_ecoli.xlsx"))

