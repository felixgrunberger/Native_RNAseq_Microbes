###########################################################################
###########################################################################
###
### ANALYSIS OF MODIFIED BASES BASED ON WRONGLY ASSIGNED POSTITIONS
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

#...................................pileup base frequency
pile_freq_base <- function(input_bam_file_location, coord_left, coord_right, fasta_file, selected_strand){
  
  fasta <- readDNAStringSet(fasta_file)
  genomename <- str_split_fixed(names(fasta)[1], " ", 2)[,1]
  names(fasta) <- "genome"
  
  param <- ScanBamParam(which=GRanges(genomename, IRanges(start=coord_left, end=coord_right)))
  allReads_pileup <- pileup(input_bam_file_location, scanBamParam = param, pileupParam=PileupParam(max_depth=8000000, min_mapq=0, min_base_quality=0))
  
  allReads_pileup_table <- allReads_pileup %>%
    dplyr::filter(strand == selected_strand) %>%
    rowwise() %>%
    mutate(correct_base = as.character(fasta$genome[pos])) %>%
    group_by(pos) %>%
    mutate(nucleotide = as.character(nucleotide)) %>%
    mutate(nuc_freq = count/sum(count) * 100, 
           category = ifelse(nucleotide == correct_base, "correct", 
                             ifelse(nucleotide == "-", "deletion", "wrong")))
  allReads_pileup_table$category <- factor(allReads_pileup_table$category, levels(as.factor(allReads_pileup_table$category))[c(3,2,1)])
  
  return(allReads_pileup_table)
}

base_annotation <- function(coord_left, coord_right, fasta_file){
  fasta <- readDNAStringSet(fasta_file)
  names(fasta) <- "genome"
  sequence <- as.character(fasta$genome[coord_left:coord_right])
  return(sequence)
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CALCULATE FREQUENCIES FOR SMALL PLOTS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................HVO WT m6A
hvo_wt_freq_base <- pile_freq_base(input_bam_file_location = here("data/mapped_data/hvo_tex.bam"),
                                   coord_left = 1599630, 
                                   coord_right = 1599650, 
                                   fasta_file = here("data/genome_data/hvo.fasta"),
                                   selected_strand = "+")

#...................................HVO dKSGA m6A
hvo_dksga_freq_base <- pile_freq_base(input_bam_file_location = here("data/mapped_data/hvo_notex_dksga.bam"),
                                      coord_left = 1599630, 
                                      coord_right = 1599650, 
                                      fasta_file = here("data/genome_data/hvo.fasta"),
                                      selected_strand = "+")



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CALCULATE FREQUENCIES FOR COMPLETE 16S
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#.................................heatmap plotting
three_color_npg <- c(pal_npg()(10)[1],
                     pal_npg()(10)[4],
                     pal_npg()(10)[7]) 

#...................................HVO WT 
hvo_whole16 <- pile_freq_base(input_bam_file_location = here("data/mapped_data/hvo_tex.bam"),
                              coord_left = 1598192, 
                              coord_right = 1599672, 
                              fasta_file = here("data/genome_data/hvo.fasta"),
                              selected_strand = "+") 
n <- hvo_whole16 %>% dplyr::filter(category %in% c("correct")) %>% dplyr::filter(nuc_freq < 70)
nrow(n)/(1599672-1598192)*100


#...................................PFU WT 
pfu_whole16 <- pile_freq_base(input_bam_file_location = here("data/mapped_data/pfu_tex.bam"),
                              coord_left = 120675, 
                              coord_right = 122193, 
                              fasta_file = here("data/genome_data/pfu.fasta"),
                              selected_strand = "+") 
n <- pfu_whole16 %>% dplyr::filter(category %in% c("correct")) %>% dplyr::filter(nuc_freq < 70)
nrow(n)/(122193-120675)*100

#...................................ECOLI WT 
ecoli_whole16 <- pile_freq_base(input_bam_file_location = here("data/mapped_data/ecoli_tex.bam"),
                                coord_left = 223771, 
                                coord_right = 225312, 
                                fasta_file = here("data/genome_data/ecoli.fasta"),
                                selected_strand = "+") 
n <- ecoli_whole16 %>% dplyr::filter(category %in% c("correct")) %>% dplyr::filter(nuc_freq < 70)
nrow(n)/(225312-223771)*100


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# GET SEQUENCE
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................HVO m6A
m6a_seq <- base_annotation(coord_left = 1599630, 
                           coord_right = 1599650, 
                           fasta_file = here("data/genome_data/hvo.fasta"))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOT SMALL 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...............................HVO m6A WT (Supplementary Fig. 15a)
pdf(here("figures/hvo_wt_modified_m6A_wrong_nucleotide.pdf"),
    width = 16, height = 3, paper = "special",onefile=FALSE)
ggplot(data = hvo_wt_freq_base, aes(x = pos, y = nuc_freq, fill = category)) +
  geom_col(width = 0.97) +
  scale_x_continuous(expand = c(0,0), breaks = 1599630:1599650, labels = strsplit(m6a_seq, "*")[[1]]) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("") +
  ylab("Frequency (%)") +
  scale_fill_manual(values = three_color_npg) +
  theme_Publication_white()
dev.off()

#...............................HVO m6A dKSGA (Supplementary Fig. 15b)
pdf(here("figures/hvo_dkska_modified_m6A_wrong_nucleotide.pdf"),
    width = 16, height = 3, paper = "special",onefile=FALSE)
ggplot(data = hvo_dksga_freq_base, aes(x = pos, y = nuc_freq, fill = category)) +
  geom_col(width = 0.97) +
  xlab("") +
  ylab("Frequency (%)") +
  scale_x_continuous(expand = c(0,0), breaks = 1599630:1599650, labels = strsplit(m6a_seq, "*")[[1]]) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = three_color_npg) +
  theme_Publication_white()
dev.off()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SUMMARY BARPLOT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................SUMMARY BARPLOT 
compare_table <- bind_rows(hvo_whole16 %>% group_by(category) %>% summarise(n = sum(count)) %>% ungroup() %>% mutate(freq = 100*n/(sum(n))) %>% mutate(group = "hvo"),
                           pfu_whole16 %>% group_by(category) %>% summarise(n = sum(count)) %>% ungroup() %>% mutate(freq = 100*n/(sum(n))) %>% mutate(group = "pfu"),
                           ecoli_whole16 %>% group_by(category) %>% summarise(n = sum(count)) %>% ungroup() %>% mutate(freq = 100*n/(sum(n))) %>% mutate(group = "ecoli"))

#...................................regroup
compare_table$group <- factor(compare_table$group, levels(as.factor(compare_table$group))[c(3,2,1)])

#...................................calculate
compare_table %>%
  dplyr::filter(category == "correct") %>%
  group_by(group, category) %>%
  summarise(wrong = 100 - freq)

#...................................plot (Supplementary Fig. 15c)  
pdf(here("figures/whole16S_wrong_base_all.pdf"), 
    width = 8, height = 3, paper = "special",onefile=FALSE)
ggplot(data = compare_table, aes(x = group, y = freq, fill = category)) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("") +
  ylab("Frequency (%)") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  scale_fill_manual(values = three_color_npg) 
dev.off()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOT FULL
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................HVO WT (Supplementary Fig. 15e)  
pdf(here("figures/hvo_whole16S_wrong_base.pdf"), 
    width = 16, height = 3, paper = "special",onefile=FALSE)
ggplot(data = hvo_whole16, aes(x = pos, y = nuc_freq, fill = category, color = category)) +
  geom_histogram(stat="identity", size = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Genomic coordinate") +
  ylab("Frequency (%)") +
  ggtitle("H. volcanii") +
  scale_fill_manual(values = three_color_npg) +
  scale_color_manual(values = three_color_npg) +
  theme_Publication_white()
dev.off()

#...................................PFU WT (Supplementary Fig. 15f)  
pdf(here("figures/pfu_whole16S_wrong_base.pdf"), 
    width = 16, height = 3, paper = "special",onefile=FALSE)
ggplot(data = pfu_whole16, aes(x = pos, y = nuc_freq, fill = category, color = category)) +
  geom_histogram(stat="identity", size = 1) +geom_bar(stat="identity", size = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Genomic coordinate") +
  ylab("Frequency (%)") +
  ggtitle("P. furiosus") +
  scale_fill_manual(values = three_color_npg) +
  scale_color_manual(values = three_color_npg) +
  theme_Publication_white()
dev.off()

#...................................ECOLI WT (Supplementary Fig. 15d)  
pdf(here("figures/ecoli_whole16S_wrong_base.pdf"), 
    width = 16, height = 3, paper = "special",onefile=FALSE)
ggplot(data = ecoli_whole16, aes(x = pos, y = nuc_freq, fill = category, color = category)) +
  geom_histogram(stat="identity", size = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Genomic coordinate") +
  ylab("Frequency (%)") +
  ggtitle("E. coli") +
  scale_fill_manual(values = three_color_npg) +
  scale_color_manual(values = three_color_npg) +
  theme_Publication_white()
dev.off()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# COMPARISON TO TEXTBOOK MODIFIED SITES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#..................................real modification sites enriched for wrong base
#................................E. COLI
interesting_postions_ecoli <- c(223771+516, 223771+527, 223771 + 966, 223771 + 967, 223771 + 1207, 223771 + 1402, 223771 + 1407, 223771 + 1498, 223771 + 1516, 223771 + 1518, 223771 + 1519)

list_int_ecoli <- list()
for (i in 1:length(interesting_postions_ecoli)){
  list_int_ecoli <- append(list_int_ecoli, (interesting_postions_ecoli[i]-4):(interesting_postions_ecoli[i]+4))
}

ecoli_whole16_groups <- ecoli_whole16 %>%
  dplyr::mutate(group = ifelse(pos %in% list_int_ecoli, "mod", "else")) %>%
  dplyr::filter(category %in% c("wrong", "deletion")) %>%
  mutate(organism = "Ecoli")


#................................P. FURIOSUS
interesting_postions_pfu <- c(120675+938,120675+1379,120675+1469, 120675+1486,120675+1487)

list_int_pfu <- list()

for (i in 1:length(interesting_postions_pfu)){
  list_int_pfu <- append(list_int_pfu, (interesting_postions_pfu[i]-4):(interesting_postions_pfu[i]+4))
}

pfu_whole16_groups <- pfu_whole16 %>%
  dplyr::mutate(group = ifelse(pos %in% list_int_pfu, "mod", "else")) %>%
  dplyr::filter(category %in% c("wrong", "deletion")) %>%
  mutate(organism = "Pfu")

#................................H. VOLCANII
interesting_postions_hvo <- c(1598192+910, 1598192+1352, 1598192+1432, 1598192+1450,1598192+1451)
list_int_hvo <- list()
for (i in 1:length(interesting_postions_hvo)){
  list_int_hvo <- append(list_int_hvo, (interesting_postions_hvo[i]-4):(interesting_postions_hvo[i]+4))
}
hvo_whole16_groups <- hvo_whole16 %>%
  dplyr::mutate(group = ifelse(pos %in% list_int_hvo, "mod", "else")) %>%
  dplyr::filter(category %in% c("wrong", "deletion")) %>%
  mutate(organism = "Hvo")

#..................................rbind all data
whole_16_all <- rbindlist(list(ecoli_whole16_groups, pfu_whole16_groups, hvo_whole16_groups))

#..................................plot (Supplementary Fig. 15g)
color2_npg <- c("#00A087", "#4DBBD5")

pdf(here("figures/wrong_base_statistics.pdf"),
    width = 8, height = 8, paper = "special",onefile=FALSE)
ggplot(data = whole_16_all, aes(y = nuc_freq, color = group, fill = group, x = category)) +
  geom_boxplot(notch = TRUE, notchwidth = 0.8, outlier.alpha = 0.2, color = "black") +
  facet_grid(~organism) +
  stat_compare_means(label =  "p.signif", label.x = 1.5, method = "t.test") +
  theme_Publication_white() +
  ylab("Frequency (%)") +
  scale_fill_manual(values = color2_npg) +
  scale_color_manual(values = color2_npg) +
  theme(panel.grid.major.x = element_blank())
dev.off()
