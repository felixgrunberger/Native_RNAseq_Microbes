###########################################################################
###########################################################################
###
### UTR 16S/23S PROCESSING SITES IN ESCHERICHIA COLI                                    
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
get_bam_input <- function(org){
  input_bam_file <- paste(here::here("data/mapped_data/"),org, ".bam", sep = "")
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

#...................................two color
two_color_npg <- rev(c("#3C5488", "#00A087"))

#...................................modify bam table
mod_bam_all <- function(input_table, set_name){
  return(input_table %>%
           mutate(UTR5_16 = ifelse(hard_l == 0 & hard_r == 0 & soft_l < 500, (start - 12) - gff$start[11],(start) - gff$start[11]),
                  UTR3_16 = (end - gff$end[11])) %>%
           mutate(UTR5_23 = ifelse(hard_l == 0 & hard_r == 0 & soft_l < 500, (start - 12) - gff$start[12],(start) - gff$start[12]),
                  UTR3_23 = (end - gff$end[12])) %>%
           mutate(UTR5_5 = ifelse(hard_l == 0 & hard_r == 0 & soft_l < 500, (start - 12) - gff$start[13],(start) - gff$start[13]),
                  UTR3_5 = (end - gff$end[13])) %>%
           mutate(sample = set_name)
  )
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................genome data
gff <- read.gff(here("data/genome_data/ecoli.gff")) %>%
  dplyr::filter(type == "rRNA")

fasta <- readDNAStringSet(here("data/genome_data/ecoli.fasta"))

#...................................input data
tex_data           <- get_bam_input("ecoli_tex")
notex_data         <- get_bam_input("ecoli_notex")

#...................................add 12 or not depended on circular or not
tex_data_f2   <- mod_bam_all(tex_data, "tex")
notex_data_f2 <- mod_bam_all(notex_data, "notex")

#...................................merge_data
all_data <- rbind(tex_data_f2, notex_data_f2)

#...................................reorder
all_data$sample <-  factor(all_data$sample, levels = rev(c("tex", "notex")))

#...................................plot 16S start position nucleotide
pdf(here("figures/200503_ecoli_16S_TSS.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = all_data, aes(x = UTR5_16 , fill = sample, y = sample)) +
  geom_density_ridges2(aes(height = ..ndensity..),alpha = 1, size = 0.1, scale = 0.95,
                       stat = "binline", draw_baseline = FALSE, bins = 150, color = NA) +
  scale_x_continuous(limits = c(-320,100), 
                     breaks = c(-293,-175,-115,-66,0), 
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("5`utr [nt]") +
  ggtitle("") +
  scale_fill_manual(values = two_color_npg) +
  scale_color_manual(values = two_color_npg) +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off() 

#...................................plot 16S end position nucleotide
pdf(here("figures/200503_ecoli_16S_TTS.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = all_data, aes(x = UTR3_16 , fill = sample, y = sample)) +
  geom_density_ridges2(aes(height = ..ndensity..),alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 200, color = NA) +
  scale_x_continuous(limits = c(-100,100), 
                     breaks = c(0,33,100), 
                     labels = c(0,33,100),
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("3`utr [nt]") +
  ggtitle("") +
  scale_fill_manual(values = two_color_npg) +
  scale_color_manual(values = two_color_npg) +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()

#...................................plot 23 start position nucleotide
pdf(here("figures/200503_ecoli_23S_TSS.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = all_data, aes(x = UTR5_23 , fill = sample, y = sample)) +
  geom_density_ridges2(aes(height = ..ndensity..),alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 150, color = NA) +
  scale_x_continuous(limits = c(-50,100), 
                     breaks = c(-7,0,50), 
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("5`utr [nt]") +
  ggtitle("") +
  scale_fill_manual(values = two_color_npg) +
  scale_color_manual(values = two_color_npg) +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()

#...................................plot 23S end position nucleotide
pdf("figures/200503_ecoli_23S_TTS.pdf", 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = all_data, aes(x = UTR3_23 , fill = sample, y = sample)) +
  geom_density_ridges2(aes(height = ..ndensity..),alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 150, color = NA) +
  scale_x_continuous(limits = c(-50,50), 
                     breaks = c(0,7,9), 
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("5`utr [nt]") +
  ggtitle("") +
  scale_fill_manual(values = two_color_npg) +
  scale_color_manual(values = two_color_npg) +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()
