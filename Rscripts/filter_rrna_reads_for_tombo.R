###########################################################################
###########################################################################
###
### FILTER READS TOMBO
###
###########################################################################
###########################################################################


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(here)
source(here("Rscripts/load_libraries.R"))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load bam file
get_bam_input <- function(org, treatment){
  input_bam_file <- paste(here::here("data/mapped_data/"),org, "_", treatment, ".bam", sep = "")
  allReads <- readGAlignments(input_bam_file, use.names = T, param = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), tag=c("NM"), what=c("mapq", "flag")))
  allReads_table <- GenomicAlignments::as.data.frame(allReads) %>%
    mutate(minion_read_name = names(allReads)) 
  
  # > get read sequences
  param <- ScanBamParam(
    flag=scanBamFlag(isUnmappedQuery=FALSE),
    what="seq")
  res0 <- scanBam(input_bam_file,param = param)[[1]] # always list-of-lists
  allReads_sequence <- res0[["seq"]]                 # query widths
  allReads_sequence_table <- as.list(as.character(allReads_sequence))
  allReads_table$sequence <- unlist(allReads_sequence_table)
  allReads_table$n_char   <- nchar(allReads_table$sequence[1:length(allReads_table$sequence)])
  left  <- paste(str_split_fixed(string = allReads_table$cigar, pattern = "M", n = 2)[,1],"M", sep = "")
  right <- paste(str_split_fixed(string = allReads_table$cigar, pattern = "M", n = 2)[,2],"1M", sep = "")
  
  #................................calculate cigar tables / SOFT AND HARD CLIPPING!!!
  allReads_table$soft_l <- as_tibble(cigarOpTable(left))$S
  allReads_table$hard_l <- as_tibble(cigarOpTable(left))$H
  allReads_table$soft_r <- as_tibble(cigarOpTable(right))$S
  allReads_table$hard_r <- as_tibble(cigarOpTable(right))$H
  return(allReads_table)
  
}

#...................................load circ bam file
get_circ_bam_input <- function(org){
  input_bam_file <- paste(here::here("data/nanopore_reanalysis/hvo/"),org, "_sorted.bam", sep = "")
  allReads <- readGAlignments(input_bam_file, use.names = T, param = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), tag=c("NM"), what=c("mapq", "flag")))
  allReads_table <- GenomicAlignments::as.data.frame(allReads) %>%
    mutate(minion_read_name = names(allReads)) 
  
  # > get read sequences
  param <- ScanBamParam(
    flag=scanBamFlag(isUnmappedQuery=FALSE),
    what="seq")
  res0 <- scanBam(input_bam_file,param = param)[[1]] # always list-of-lists
  allReads_sequence <- res0[["seq"]]                 # query widths
  allReads_sequence_table <- as.list(as.character(allReads_sequence))
  allReads_table$sequence <- unlist(allReads_sequence_table)
  allReads_table$n_char   <- nchar(allReads_table$sequence[1:length(allReads_table$sequence)])
  left  <- paste(str_split_fixed(string = allReads_table$cigar, pattern = "M", n = 2)[,1],"M", sep = "")
  right <- paste(str_split_fixed(string = allReads_table$cigar, pattern = "M", n = 2)[,2],"1M", sep = "")
  
  #................................calculate cigar tables / SOFT AND HARD CLIPPING!!!
  allReads_table$soft_l <- as_tibble(cigarOpTable(left))$S
  allReads_table$hard_l <- as_tibble(cigarOpTable(left))$H
  allReads_table$soft_r <- as_tibble(cigarOpTable(right))$S
  allReads_table$hard_r <- as_tibble(cigarOpTable(right))$H
  return(allReads_table)
}



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#...................................gff
gff <- read.gff(here::here("data/genome_data/hvo.gff")) %>%
  dplyr::filter(type == "rRNA")

#...................................fasta
fasta <- readDNAStringSet(here::here("data/genome_data/hvo.fasta")) 

#...................................load bam file / normal & circ
bam_notex          <- get_bam_input(org = "hvo", treatment = "notex")
bam_notex_circ    <- get_circ_bam_input(org = "hvo/16S_opencirc_hvo")

#...................................detect_primary_pre_rRNA
detect_primary_pre_rRNA <- function(input_bam){
  output_readnames <- input_bam %>%
    ungroup() %>%
    dplyr::filter(start > (gff$start[1]-500), 
                  start < (gff$start[1]+12),
                  hard_l == 0,
                  soft_l < 100,
                  hard_r == 0,
                  njunc == 0,
                  end > (gff$end[1]-5)) %>%
    distinct(minion_read_name, .keep_all = T) %>%
    mutate(read_group = "primary_pre_rRNA") %>%
    dplyr::select(minion_read_name, read_group)
}

#...................................detect_mature_rRNA
detect_mature_rRNA <- function(input_bam){
  output_readnames <- input_bam %>%
    ungroup() %>%
    dplyr::filter(start > (gff$start[1]), 
                  start < (gff$end[1]),
                  hard_l == 0,
                  soft_l < 50,
                  hard_l == 0,
                  hard_r == 0,
                  njunc == 0,
                  end > (gff$end[1]-100),
                  end < (gff$end[1]+10)) %>%
    distinct(minion_read_name, .keep_all = T) %>%
    mutate(read_group = "mature_rRNA") %>%
    dplyr::select(minion_read_name, read_group)
}

#...................................detect_closed_circular
detect_closed_circ_rRNA <- function(input_bam){
  output_readnames <- input_bam %>%
    dplyr::filter(flag == 0,
                  hard_r == 0,
                  start < (570-100),
                  end > (570+108+5)) %>%
    distinct(minion_read_name, .keep_all = TRUE) %>%
    mutate(read_group = "closed_circ_rRNA") %>%
    dplyr::select(minion_read_name, read_group)
}

#...................................detect_open_circular
detect_open_circ_rRNA <- function(input_bam){
  output_readnames <- input_bam %>%
    dplyr::filter(flag == 0,
                  hard_r == 0,
                  start < (570-100),
                  end %in% (570+108-5):(570+108+5)) %>%
    distinct(minion_read_name, .keep_all = TRUE) %>%
    mutate(read_group = "open_circ_rRNA") %>%
    dplyr::select(minion_read_name, read_group)
}


#.................................primary pre-rRNA
notex_primary_pre_rRNA <- detect_primary_pre_rRNA(bam_notex)

#.................................closed circ pre-rRNA
notex_closed_circ_rRNA   <- detect_closed_circ_rRNA(bam_notex_circ)

#.................................open circ pre-rRNA
notex_open_circ_rRNA <- detect_open_circ_rRNA(bam_notex_circ)

#.................................mature rRNA
notex_mature_rRNA <- detect_mature_rRNA(bam_notex)

#.................................write to list notex
write_lines(notex_primary_pre_rRNA$minion_read_name, 
            path = here("data/reanalysis_tombo/hvo/primary_pre_rRNA_wt.lst"))
write_lines(notex_closed_circ_rRNA$minion_read_name, 
            path = here("data/reanalysis_tombo/hvo/closed_circ_rRNA_wt.lst"))
write_lines(notex_open_circ_rRNA$minion_read_name, 
            path = here("data/reanalysis_tombo/hvo/open_circ_rRNA_wt.lst"))
write_lines(notex_mature_rRNA$minion_read_name, 
            path = here("data/reanalysis_tombo/hvo/mature_rRNA_wt.lst"))


#.................................subset fast5 
# >>>>>>>>>>> mk dirs
dir_delete(here::here("data/reanalysis_tombo/hvo/primary_pre_rRNA_reads"))
dir_create(here::here("data/reanalysis_tombo/hvo/primary_pre_rRNA_reads"))
dir_delete(here::here("data/reanalysis_tombo/hvo/closed_circ_rRNA_reads"))
dir_create(here::here("data/reanalysis_tombo/hvo/closed_circ_rRNA_reads"))
dir_delete(here::here("data/reanalysis_tombo/hvo/open_circ_rRNA_reads"))
dir_create(here::here("data/reanalysis_tombo/hvo/open_circ_rRNA_reads"))
dir_delete(here::here("data/reanalysis_tombo/hvo/mature_rRNA_reads"))
dir_create(here::here("data/reanalysis_tombo/hvo/mature_rRNA_reads"))
dir_create(here::here("data/reanalysis_tombo/hvo/mature_rRNA_reads2"))
dir_create(here::here("data/reanalysis_tombo/hvo/primary_pre_rRNA_reads2"))

# >>>>>>>>>>> get all files 
list_of_all_files_hvo <- list.files("..list_to_directory_with_single_fast5_files..", recursive = TRUE, full.names = TRUE)

list_fast5_hvo <- list_of_all_files_hvo %>%
  tibble::enframe(name = NULL) %>%
  dplyr::rename(long_fast5 = 1) %>%
  mutate(short_fast5 = str_split_fixed(long_fast5, "/", 9)[,9])

# >>>>>>>>>>> add path to read categories
mod_file_path <- function(input_dataset){
  paste(input_dataset$minion_read_name, ".fast5", sep = "") %>%
    tibble::enframe(name = NULL) %>%
    dplyr::rename(short_fast5 = value) %>%
    inner_join(list_fast5_hvo, by = "short_fast5") %>%
    dplyr::select(long_fast5) %>%
    mutate(long_fast5 = as.character(long_fast5)) 
}

#..........NOTEX
# > primary pre rRNA
notex_primary_pre_rRNA_filepaths <- mod_file_path(notex_primary_pre_rRNA)

# > closed circ rRNA
notex_closed_circ_rRNA_filepaths <- mod_file_path(notex_closed_circ_rRNA)

# > open circ rRNA
notex_open_circ_rRNA_filepaths   <- mod_file_path(notex_open_circ_rRNA)

# > mature rRNA
notex_mature_rRNA_filepaths      <- mod_file_path(notex_mature_rRNA)


# >>>>>>>>>>> copy files
# > primary pre rRNA
for (i in notex_primary_pre_rRNA_filepaths){
  fs::file_copy(path = i, new_path = here::here("data/reanalysis_tombo/hvo/primary_pre_rRNA_reads2"),overwrite = T) 
}

# > closed circ rRNA
for (i in notex_closed_circ_rRNA_filepaths){
  fs::file_copy(path = i, new_path = here::here("data/reanalysis_tombo/hvo/closed_circ_rRNA_reads"),overwrite = T) 
}

# > open circ rRNA
for (i in notex_open_circ_rRNA_filepaths){
  fs::file_copy(path = i, new_path = here::here("data/reanalysis_tombo/hvo/open_circ_rRNA_reads"),overwrite = T) 
}

# > mature rRNA
for (i in notex_mature_rRNA_filepaths){
  fs::file_copy(path = i, new_path = here::here("data/reanalysis_tombo/hvo/mature_rRNA_reads"),overwrite = T) 
}


#.................................subset fastq outside of R
#...............................NOTEX
seqtk subseq /data/fastq_data/hvo_notex.fastq /data/reanalysis_tombo/hvo/primary_pre_rRNA_wt.lst > /data/reanalysis_tombo/hvo/primary_pre_rRNA_wt.fastq
seqtk subseq /data/fastq_data/hvo_notex.fastq /data/reanalysis_tombo/hvo/primary_pre_rRNA2_wt.lst > /data/reanalysis_tombo/hvo/primary_pre_rRNA2_wt.fastq
seqtk subseq /data/fastq_data/hvo_notex.fastq /data/reanalysis_tombo/hvo/closed_circ_rRNA_wt.lst > /data/reanalysis_tombo/hvo/closed_circ_rRNA_wt.fastq
seqtk subseq /data/fastq_data/hvo_notex.fastq /data/reanalysis_tombo/hvo/open_circ_rRNA_wt.lst > /data/reanalysis_tombo/hvo/open_circ_rRNA_wt.fastq

#.................................mapping outside of R            
directory=/data/reanalysis_tombo/hvo
fasta=/data/genome_data/hvo.fasta
output=$directory/mapped_data

for file in ${directory}/*.fastq
do 
xbase=${file##*/}
  filename=${xbase%%.*}
  minimap2 -p 0.99 -ax splice -k14  --MD -uf ${fasta} $file > ${output}"/"$filename".sam"
  samtools view -bS ${output}"/"$filename".sam" -o ${output}"/"$filename".bam"
  samtools sort ${output}"/"$filename".bam" -o ${output}"/"$filename"_sorted.bam"
  samtools index ${output}"/"$filename"_sorted.bam"
  echo $filename "mapping finished"  
done

      