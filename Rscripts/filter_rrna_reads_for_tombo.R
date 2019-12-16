###########################################################################
###########################################################################
###
### TOMBO OUTPUT FRACTION ANALYSIS
###
###########################################################################
###########################################################################


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES AND PLOTTING FUNCTION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(here)
source(here("Rscripts/load_libraries.R"))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load genome gff file
hvo_gff <- read.gff(here("data/genome_data/hvo.gff")) %>%
  dplyr::filter(type == "rRNA")

#...................................load dataset (tidy)
load(here("data/tidy_data/hvo_notex_id_table"))

#...................................filter reads that map to 16S (1 loci), calculated UTR5 and UTR3 and group according to 5´end 
full_id_table_16S_processsing_groups <- full_id_table %>%
  ungroup() %>%
  dplyr::filter(locus_name == "16S", gene == "rna27") %>%
  dplyr::mutate(UTR5 = ifelse(strand_gene == "+" ,start_gene - start,end - end_gene),
                UTR3 = ifelse(strand_gene == "+", end - end_gene, start_gene - start),
                gene = as.factor(gene),
                group = ifelse(UTR5 > -12, "unprocessed", "processed"), # correct for 12 nt´s missing
                minion_read_name = paste(minion_read_name, ".fast5", sep = ""))

#...................................select reads that span a certain region of the 16S rRNA end to enable analysis
full_id_table_16S_processsing_groups_final <- full_id_table_16S_processsing_groups %>%
  dplyr::filter(gene == "rna27", end > 1599660) %>%
  dplyr::select(minion_read_name, group) 

#...................................write minion read names belonging to different groups to files
write.table(x = full_id_table_16S_processsing_groups_final$minion_read_name[full_id_table_16S_processsing_groups_final$group == "unprocessed"], 
            file = here("/data/tombo_data/hvo/unprocessed_hvonotex.txt"),quote = F, row.names = F, col.names = F)
write.table(x = full_id_table_16S_processsing_groups_final$minion_read_name[full_id_table_16S_processsing_groups_final$group == "processed"], 
            file = here("/data/tombo_data/hvo/processed_hvonotex.txt"),quote = F, row.names = F, col.names = F)

#...................................path to all single-read FAST5 files
list_of_all_files_hvo <- list.files(<path_to_file>, recursive = TRUE, full.names = TRUE) %>%
  tibble::enframe(name = NULL) %>%
  dplyr::rename(long_fast5 = 1) %>%
  mutate(short_fast5 = str_split_fixed(long_fast5, "/", 8)[,8])

#...................................
file_processed <- read.delim(here("data/tombo_data/hvo/processed_hvonotex.txt"), header = F) %>%
  mutate(V1 = as.character(V1)) %>%
  dplyr::rename(short_fast5 = V1) %>%
  inner_join(list_of_all_files_hvo, by = "short_fast5") %>%
  dplyr::select(long_fast5) %>%
  mutate(long_fast5 = as.character(long_fast5))

file_unprocessed <- read.delim(here("data/tombo_data/hvo/processed_hvonotex.txt"), header = F) %>%
  mutate(V1 = as.character(V1)) %>%
  dplyr::rename(short_fast5 = V1) %>%
  inner_join(list_of_all_files_hvo, by = "short_fast5") %>%
  dplyr::select(long_fast5) %>%
  mutate(long_fast5 = as.character(long_fast5))

output_folder_processed   <- here("data/tombo_data/hvo/processed_reads/reads")
output_folder_unprocessed <- here("data/tombo_data/hvo/unprocessed_reads/reads")

dir.create(output_folder_processed)
dir.create(output_folder_unprocessed)


for (i in file_processed){
  file.copy(from = i, to = here("/data/tombo_data/hvo/processed_reads/reads"), recursive = TRUE)
}

for (i in file_unprocessed){
  file.copy(from = i, to = here("/data/tombo_data/hvo/unprocessed_reads/reads"), recursive = TRUE)
}


