###########################################################################
###########################################################################
###
### RUN STATISTICS FOR SUPPLEMENTARY TABLE 1
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

#...................................load and mutate summary file
wrapper_summary_table <- function (input_summary_file, identifier, barcodes_used){
  
  # read in summary file from guppy
  summary_table <- fread(input_summary_file) 
  
  # how to deal with different output from poreplex
  summary_table$group <- NA

  # add sequencing data set as 'identifier'
  summary_table_selected <- summary_table %>%
    mutate(group = ifelse(is.na(group),identifier, group)) 
  
  return(summary_table_selected)
}

#...................................calculate bases on CDS
cds_depth_calculator <- function(input_gff){
  gff_length_table <- read.gff(input_gff) %>%
    dplyr::filter(type == "CDS") %>%
    mutate(length = abs(start - end))
  return(sum(gff_length_table$length))
}

#...................................calculate enolase quality from mapped file
enolase_quality_finder <- function(input_bam_file, input_summary_file, seq_set){
  summary_file <- fread(input_summary_file)
  p4 <- ScanBamParam(tag=c("NM", "MD"), what="mapq")
  allReads <- readGAlignments(input_bam_file, use.names = T, param = p4)
  allReads_table <- GenomicAlignments::as.data.frame(allReads) %>%
    as_tibble() %>%
    mutate(minion_read_name = names(allReads))
  bam_summary <- left_join(allReads_table, summary_file, by = c("minion_read_name" = "read_id"))
  bam_summary$aligned_reads <- NA
  bam_summary$aligned_reads <- unlist(lapply(explodeCigarOpLengths(bam_summary$cigar, ops = c("M", "I")), function(x) sum(x)))
  bam_summary$identity <- NA
  bam_summary <- bam_summary %>%
    mutate(identity = (1 - NM/aligned_reads)*100,
           mapped_to = "control",
           sequencing_set = seq_set,
           mapped_type = "CDS")
  return(bam_summary)
}

#...................................modify id table output
modify_id_table <- function(id_table_input, name){
  return(id_table_input %>%
           ungroup() %>%
           distinct(minion_read_name, .keep_all = TRUE) %>%
           mutate(mapped_to = "genome") %>%
           mutate(sequencing_set = name)) 
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA ENOLASE
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................paths to summary files from guppy
summary_files <- paste(here("data/summary_data/"), list.files(here("data/summary_data/"),pattern = ".txt"), sep = "")

#...................................paths to mapped enolase files
enolase_files <- paste(here("data/enolase_data/"), list.files(here("data/enolase_data/"),pattern = ".bam"), sep = "")

#...................................get sample names
sample_names <- unlist(lapply(summary_files, FUN=function(x){str_split_fixed(str_split_fixed(x, "_seq", 2)[1],"summary_data/",2)[2]}))

#...................................calculate enolase tables
enolase_table <- data.frame()

for (i in seq_along(sample_names)){
  enolase_table <- rbindlist(list(enolase_table, enolase_quality_finder(enolase_files[i], summary_files[i],sample_names[i])))
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY RAW READ
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load summary data and modify them
summary_files <- paste(here("data/summary_data/"), list.files(here("data/summary_data/")), sep = "")

for (i in seq_along(summary_files)){
  sample_name <- str_split_fixed(str_split_fixed(summary_files[i], "_seq", 2)[1],"summary_data/",2)[2]
  table_name  <- paste(sample_name, "_table", sep = "")
  assign(table_name, wrapper_summary_table(summary_files[i], sample_name))
}

#...................................merge all tables into one
summary_table <- rbindlist(list(pfu_tex_table,
                                pfu_notex_table,
                                ecoli_tex_table, 
                                ecoli_notex_table, 
                                hvo_tex_table, 
                                hvo_notex_table,
                                hvo_notex_lowsalt_table,
                                hvo_notex_dksga_table)) %>% 
  mutate(sequencing_set = as.factor(group),
         start_time = as.numeric(start_time), 
         sequence_length_template = as.numeric(sequence_length_template),
         mean_qscore_template = as.numeric(mean_qscore_template),
         type = ifelse(calibration_strand_score > 0, "enolase", "genome"))

#...................................reorder levels
summary_table$sequencing_set <-  factor(summary_table$sequencing_set, 
                                        levels = c("ecoli_tex","ecoli_notex",
                                                       "pfu_tex", "pfu_notex",
                                                       "hvo_tex", "hvo_notex",
                                                       "hvo_notex_lowsalt", "hvo_notex_dksga"))

#...................................calculate read statistics
#.................................number of reads/bases, median read lengths/qualities
summary_stats <- summary_table %>%
  group_by(sequencing_set) %>%
  summarise(number_of_reads = n(),
            number_of_bases = sum(sequence_length_template, na.rm = TRUE),
            median_rawread_length = median(sequence_length_template, na.rm = TRUE),
            median_rawread_quality = median(mean_qscore_template, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(-sequencing_set)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA MAPPED TO GENOME
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................input: pre-calculated single-read tidy files
id_files <- paste(here("data/tidy_data/"), list.files(here("data/tidy_data/"),pattern = "_id_table"), sep = "")
sample_names <- str_split_fixed(str_split_fixed(id_files,"tidy_data/",2)[,2], "_id_table",2)[,1]
#...................................calculate genome tables
genome_table <- data.frame()

for (i in seq_along(id_files)){
  load(id_files[i])
  genome_table <- rbindlist(list(genome_table, modify_id_table(full_id_table, sample_names[i])))
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MERGE ALL DATA 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................combine enolase and genome tables
full_table <- rbindlist(list(enolase_table,genome_table), fill = TRUE)

#...................................reorder levels
full_table$sequencing_set <-  factor(full_table$sequencing_set, 
                                     levels = c("ecoli_tex","ecoli_notex",
                                                "pfu_tex", "pfu_notex",
                                                "hvo_tex", "hvo_notex",
                                                "hvo_notex_lowsalt", "hvo_notex_dksga"))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CALCULATE STATISTICS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................calculate mapped read statistics
#.................................number of mapped reads / median identity
mapped_stats <- full_table %>%
  dplyr::filter(mapped_to == "genome") %>%
  group_by(sequencing_set) %>%
  summarise(number_of_mapped_reads = n(),
            median_identity = median(identity)) %>%
  ungroup() %>%
  dplyr::select(-sequencing_set)

#.................................CDS mapping features
#...............................calculate bases on CDS
cds_pfu   <- cds_depth_calculator(here("data/genome_data/pfu.gff"))
cds_hvo   <- cds_depth_calculator(here("data/genome_data/hvo.gff"))
cds_ecoli <- cds_depth_calculator(here("data/genome_data/ecoli.gff"))

cds_stats <- full_table %>%
  mutate(cds_depth = ifelse(sequencing_set == "ecoli_notex" | sequencing_set == "ecoli_tex", cds_ecoli,
                            ifelse(sequencing_set == "pfu_notex" | sequencing_set == "pfu_tex", cds_pfu, cds_hvo))) %>%
  group_by(sequencing_set) %>%
  dplyr::filter(mapped_type == "CDS", mapped_to == "genome") %>%
  summarise(reads_mapped_to_CDS = n(),
            bases_mapped_to_CDS = sum(aligned_reads),
            sequencing_depth_CDS = sum(aligned_reads)/max(cds_depth)) %>%
  ungroup() %>%
  dplyr::select(-sequencing_set)

#.................................rRNA mapping features
rRNA_stats <- full_table %>%
  group_by(sequencing_set) %>%
  dplyr::filter(mapped_type == "rRNA", mapped_to == "genome") %>%
  summarise(reads_mapped_to_rRNA = n()) %>%
  ungroup() %>%
  dplyr::select(-sequencing_set)

#.................................enolase mapping features
enolase_stats <- full_table %>%
  group_by(sequencing_set) %>%
  dplyr::filter(mapped_type == "CDS", mapped_to == "control") %>%
  summarise(reads_mapped_to_enolase = n(),
            median_length_enolase = median(as.numeric(qwidth, na.rm = T)),
            median_identity_enolase = median(as.numeric(identity, na.rm = T))) %>%
  ungroup() %>%
  dplyr::select(-sequencing_set)

#.................................mutliplex_strategy
run_ids <- list()
run_ids$organism <- c(rep("Escherichia coli",2),rep("Pyrococcus furiosus", 2), rep("Haloferax volcanii", 4))
run_ids$strain <- c(rep("k12",2),rep("dsm3638",2),rep("h26",3),"h26/∆ksgA")
run_ids$condition <- c(rep("normal",6), "lowsalt","normal")
run_ids$flowcell_nr <- c("1", "2", "3", "2", "4", "5", "2", "5")
run_ids$tex <- c(rep(c("yes","no"),3), "no","no")
run_ids$multiplexed <- c("no","yes","no", "yes","yes", "yes", "yes", "yes")
run_ids <- as_tibble(run_ids)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# WRITE TO XLSX FILE 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................combine stats
stats_table <- format(cbind(run_ids, summary_stats,mapped_stats,cds_stats,rRNA_stats,enolase_stats) %>%
  mutate(percentage_mapped = number_of_mapped_reads/number_of_reads*100), digits = 3)

#...................................Supplementary Table 1
writexl::write_xlsx(x = stats_table, path = here("tables/Supplementary_Table_1.xlsx"))


