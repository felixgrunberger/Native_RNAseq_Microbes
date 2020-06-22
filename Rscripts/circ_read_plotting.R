###########################################################################
###########################################################################
###
### CIRC-rRNA PLOTTING
###
###########################################################################
###########################################################################


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(here)
source(here("Rscripts/load_libraries.R"))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load bam file
get_bam_input_filepath <- function(path_to_file){
  input_bam_file <- path_to_file
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

#...................................filter bam files by certain criteria that only leave circ-reads mapping
#1-500:16S_end-500:16S_end
#500-570: 16S_end-16S_3´bhb
#to overlap the read start has to be < 570 and the end > 570
filter_set <- function(input){
  input %>%
    dplyr::filter(flag == 0,
                  hard_r == 0,
                  start < (570-100),
                  #start < 570,
                  end > (570+5)) %>%
    distinct(minion_read_name, .keep_all = TRUE)
}




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................genomic data
gff <- read.gff(here("data/genome_data/hvo.gff")) %>%
  dplyr::filter(type == "rRNA")

fasta <- readDNAStringSet(here("data/genome_data/hvo.fasta"))
names(fasta)[1] <- "chr"

#...................................reads mapping to total bam
hvo_notex <- get_bam_input_filepath(here("data/mapped_data/hvo_notex.bam"))

#...................................read in the 16S region
hvo_notex_f <- hvo_notex %>%
  dplyr::filter(start %in% (gff$start[1]-500):(gff$end[1]+100),
                end %in% (gff$start[1]-500):(gff$end[1]+100)) %>%
  distinct(minion_read_name)

#...................................load bam file of reads mapping to permuted rRNA sequence
bam16circ_notex         <- get_bam_input_filepath(here("data/mapped_data/16S_circ_hvo.bam"))

#...................................filter reads to make sure that they are circular
bam16circ_notex_f             <- filter_set(bam16circ_notex)

bam16circ_tex_f2 <- bam16circ_notex_f %>%
  mutate(width = end - start - soft_l) %>%
  mutate(ldf= 1:nrow(bam16circ_notex_f))
   
#...................................plot
pdf(here("figures/single_read_hvo_notex_permuted.pdf"),
    width = 6, height = 6, paper = "special",onefile=FALSE)
ggplot(data = bam16circ_tex_f2, aes(x= ldf, ymax = end, ymin = start-soft_l)) +
  geom_linerange() +
  coord_flip() +
  geom_hline(yintercept = c(570,570+106, 500-1467), linetype = "dashed") +
  scale_y_continuous(limits = c(500-1467-10, 500+500+70+108)) +
  theme_Publication_white() +
  ylab("") +
  xlab("")
dev.off()




