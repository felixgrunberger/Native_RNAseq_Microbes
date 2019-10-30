###########################################################################
###########################################################################
###
### CONVERT BAM TABLE TO SINGLE_READ COUNT OBJECT INCLUDING ALL GUPPY READ INFORMATION
###
###########################################################################
###########################################################################


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
packages <- c("tidyverse", "here", "ggthemes", "data.table", 
              "ggExtra", "Rsamtools", "GenomicAlignments", "seqTools", 
              "Rsubread", "ape", "DT", "ggpubr", "ggridges", "ggsci")
invisible(lapply(packages, require, character.only = TRUE))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................unlist bam to datatable
.unlist <- function (x){
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#...................................wrapper function for sequencing summary input
wrapper_summary_table <- function (input_summary_file){
  fread(input_summary_file)
}

#...................................wrapper function for bam input (count with featurecounts and allow (!) multi mapping reads)
counts_wrapper <- function (input_bam_file, input_fasta_file, input_gff_file){
  
  # > read in fasta file
  fasta <- readDNAStringSet(filepath = input_fasta_file)
  
  # > use featurecounts to calculate counts of mapped reads to features (CDS, tRNA, rRNA) 
  datalist  <-  list()
  datalist2 <- list()
  datafile <- str_split(input_bam_file, "/")
  interesting_list <- c("CDS", "rRNA", "tRNA")
  
  for (i in seq_along(interesting_list)){
    name <- interesting_list[i]
    dir.create(paste(here("data/featurecounts_data_"),name, sep = ""),showWarnings = FALSE)
    datalist[[i]] <- featureCounts(allowMultiOverlap = T, files = input_bam_file, annot.ext = input_gff_file, isGTFAnnotationFile = T, GTF.featureType = name, GTF.attrType = "ID", isLongRead = T,nthreads = 8, reportReads = "CORE", reportReadsPath = paste(here("data/featurecounts_data_"),name, sep = ""))
    datalist2[[i]] <- fread(paste(here("data/featurecounts_data_"),name, "/",datafile[[1]][length(datafile[[1]])], ".featureCounts", sep = "")) %>%
      dplyr::rename(id = V1, gene = V4) %>%
      dplyr::select(id, gene) %>%
      dplyr::filter(!is.na(gene)) %>%
      mutate(mapped_type = name)
  }
  return(list(datalist, datalist2))
}

#...................................combine bam, summary, gff and featurecounts file
wrapper_bam_to_table <- function (input_bam_file, input_gff_file, input_fasta_file, datalist_input, output = c("read_ids", "gene_ids")){
  
  # > read in fasta file
  fasta <- readDNAStringSet(filepath = input_fasta_file)
  
  # > read in gff file and grep for feature numbers
  interesting_list <- c("CDS", "rRNA", "tRNA")
  gff_table <- read.gff(input_gff_file) %>%
    as_tibble() %>%
    mutate(start_gene = start, end_gene = end,strand_gene = strand) %>%
    dplyr::filter(type %in% interesting_list) %>%
    mutate(id_name = str_split_fixed(str_split_fixed(attributes, ";Parent=", 2)[,1], "ID=", 2)[,2],
           locus_name = ifelse(type == "CDS", str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], ";", 2)[,1],
                               ifelse(type == "rRNA", str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], " ", 2)[,1], 
                                      ifelse(type == "tRNA", str_split_fixed(attributes, ";product=", 2)[,2], NA )))) %>%
    dplyr::select(id_name, locus_name, start_gene, end_gene, strand_gene)
  
  # > for featurecounts calculated transcript abundacies for each gene
  if(output == "gene_ids"){
    feature_list <- c(rep(interesting_list[1],length(datalist_input[[1]]$annotation$Chr)),
                      rep(interesting_list[2],length(datalist_input[[2]]$annotation$Chr)),
                      rep(interesting_list[3],length(datalist_input[[3]]$annotation$Chr)))
    
    big_data_all <- cbind(datalist_input[[1]]$annotation, datalist_input[[1]]$counts) %>%
      rbind(cbind(datalist_input[[2]]$annotation, datalist_input[[2]]$counts)) %>%
      rbind(cbind(datalist_input[[3]]$annotation, datalist_input[[3]]$counts)) %>%
      as_tibble() %>%
      mutate(type = feature_list) %>%
      dplyr::rename(counts = 7) %>%
      rowwise() %>%
      dplyr::filter(Strand == "+" | Strand == "-") 
    
    listofdfs <- list()
    # > enable calculation for different chromosomes
    for(i in 1:length(names(fasta))){
      names(fasta) <-  str_split_fixed(names(fasta), " ", 2)[,1]
      used_chr <- str_split_fixed(names(fasta[i]), " ", 2)[,1]
      
      
      df <- big_data_all %>%
        dplyr::filter(Chr == used_chr) %>%
        rowwise() %>%
        mutate(seq = ifelse(Strand == "+" & End < length(fasta[names(fasta) == used_chr][[1]]), as.character(fasta[names(fasta) == used_chr][[1]][Start:End]),
                            ifelse(Strand == "-" & End < length(fasta[names(fasta) == used_chr][[1]]),as.character(reverseComplement(fasta[names(fasta) == used_chr][[1]][Start:End])), NA)))
      listofdfs[[i]] <- df
    }
    
    full_table_big_data <- data.frame(Reduce(rbind, listofdfs))
    big_data_all_seq_names <- left_join(full_table_big_data, gff_table, by = c("GeneID" = "id_name"))
    return(big_data_all_seq_names)
  }
  
  # > single read table output
  if(output == "read_ids"){
    
    assigned_features <- do.call(rbind,datalist_input)
    
    # >calculate correct start end end positions of mapped reads | read in BAM file with NM tag
    allReads <- readGAlignments(input_bam_file, use.names = T, param = ScanBamParam(tag=c("NM"), what="mapq"))
    
    allReads_table <- GenomicAlignments::as.data.frame(allReads) %>%
      mutate(minion_read_name = names(allReads)) %>%
      left_join(summary_table, by = c("minion_read_name" = "read_id")) 
    
    # > calculate number of aligned reads based on CIGAR operations (M,I)
    allReads_table$aligned_reads <- NA
    allReads_table$aligned_reads <- unlist(lapply(explodeCigarOpLengths(allReads_table$cigar, ops = c("M", "I")), function(x) sum(x)))
    
    # > join featurecounts table | calculate mapping identity | factorize strands
    allReads_table_filtered <- allReads_table %>%
      left_join(assigned_features, by = c("minion_read_name" = "id")) %>%
      mutate(identity = (1 - NM/aligned_reads)*100,
             length_read = qwidth)  %>%
      separate_rows(gene, sep = ",") %>%
      left_join(gff_table, by = c("gene" = "id_name")) %>%
      mutate(shortest_distance_to_gene = ifelse(abs(start-start_gene) <= max(sequence_length_template) & abs(end-end_gene)<=max(sequence_length_template), T, F)) %>%
      dplyr::filter(shortest_distance_to_gene == T) %>%
      mutate(strand = factor(strand, levels = c("+","-"))) %>%
      dplyr::filter(strand == strand_gene) %>%
      dplyr::select(seqnames, strand, qwidth, start, end, width, mapq, NM, minion_read_name, 
                    sequence_length_template, mean_qscore_template, aligned_reads, identity, gene, mapped_type, length_read,
                    start_gene, end_gene, locus_name, strand_gene) 
    
    return(allReads_table_filtered)
  }
}



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FILES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
summary_files <- paste(here("data/summary_data/"), list.files(here("data/summary_data/")), sep = "")
sample_names <- unlist(lapply(summary_files, FUN=function(x){str_split_fixed(str_split_fixed(x, "_seq", 2)[1],"summary_data/",2)[2]}))

#...................................set names
for (i in seq_along(sample_names)){
  working_directory   <- paste(here(), "/data", sep = "")
  sample_name         <- sample_names[i]
  input_summary       <- summary_files[i]
  where_to_store      <- paste(working_directory,"/tidy_data/", sample_name, sep = "")
  input_fasta         <- paste(working_directory, "/genome_data/", str_split_fixed(sample_name, "_",2)[1],".fasta", sep = "")
  input_gff           <- paste(working_directory, "/genome_data/", str_split_fixed(sample_name, "_",2)[1],".gff", sep = "")
  input_bam           <- paste(working_directory, "/mapped_data/", sample_name,".bam", sep = "")
  type_features       <- c("All", "rRNA", "CDS", "tRNA")
  
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # COMPUTE SUMMARY TABLE
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  summary_table <- wrapper_summary_table(input_summary_file = input_summary)
  
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # COMPUTE BAM TABLE
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if("barcode" %in% colnames(summary_table)){
    for (j in seq_along(used_barcodes)){
      # > for BAM plotting (output from featurecounts)
      full_counts_table <- counts_wrapper(input_bam_file   = get(paste("input_bam", used_barcodes[j], sep = "_")),
                                          input_fasta_file = get(paste("input_fasta",used_barcodes[j], sep = "_")),
                                          input_gff_file   = get(paste("input_gff",used_barcodes[j], sep = "_")))
      
      # > for BAM plotting (identity calculation, ...)
      full_id_table <- wrapper_bam_to_table(input_bam_file = get(paste("input_bam", used_barcodes[j], sep = "_")), 
                                            input_gff_file = get(paste("input_gff",used_barcodes[j], sep = "_")),
                                            input_fasta_file = get(paste("input_fasta",used_barcodes[j], sep = "_")),
                                            datalist_input = full_counts_table[[2]],
                                            output = "read_ids") 
      
      # > for BAM plotting (names of genes, ...)
      full_gene_table <- wrapper_bam_to_table(input_bam_file = get(paste("input_bam", used_barcodes[j], sep = "_")), 
                                              input_gff_file = get(paste("input_gff",used_barcodes[j], sep = "_")),
                                              input_fasta_file = get(paste("input_fasta",used_barcodes[j], sep = "_")),
                                              datalist_input = full_counts_table[[1]],
                                              output = "gene_ids")
      
      # > save as R file
      save(full_id_table, file = paste(where_to_store, "/Files/full_id_table", sample_name, "_", used_barcodes[j], sep = ""))
      save(full_gene_table, file = paste(where_to_store, "/Files/full_gene_table", sample_name, "_", used_barcodes[j], sep = ""))
    }
    
    
  } else {
    # > for BAM plotting (output from featurecounts)
    full_counts_table <- counts_wrapper(input_bam_file = input_bam,
                                        input_fasta_file = input_fasta,
                                        input_gff_file = input_gff)
    
    # > for BAM plotting (identity calculation, ...)
    full_id_table <- wrapper_bam_to_table(input_bam_file = input_bam, 
                                          input_gff_file = input_gff,
                                          input_fasta_file = input_fasta,
                                          datalist_input = full_counts_table[[2]],
                                          output = "read_ids")
    
    # > for BAM plotting (names of genes, ...)
    full_gene_table <- wrapper_bam_to_table(input_bam_file = input_bam, 
                                            input_gff_file = input_gff,
                                            input_fasta_file = input_fasta,
                                            datalist_input = full_counts_table[[1]],
                                            output = "gene_ids")
    # > save as R file
    save(full_id_table, file = paste(where_to_store, "_id_table", sep = ""))
    save(full_gene_table, file = paste(where_to_store, "_gene_table", sep = ""))
  }
}




