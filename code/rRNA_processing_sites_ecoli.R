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
packages <- c("tidyverse", "here", "ggthemes", "data.table", 
              "ggExtra", "Rsamtools", "GenomicAlignments", 
              "seqTools", "Rsubread", "ape", "DT", "ggpubr", 
              "ggridges", "ggsci")
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


#...................................function for 5´analysis
analyse_5utr <- function(id_table){
  full_id_table_genes <- id_table %>%
    ungroup() %>%
    dplyr::filter(mapped_type == "rRNA") %>%
    mutate(UTR5 = ifelse(strand_gene == "+" ,start_gene - (start - 12), (end + 12) - end_gene),
           UTR3 = ifelse(strand_gene == "+", end - end_gene, start_gene - start))
  return(full_id_table_genes)
}

#...................................get rRNA (C) locus E. coli
interesting_positions <- 3939348:3945138

analyse_rDNA_all <- function(id_table){
  full_id_table_genes <- id_table %>%
    ungroup() %>%
    dplyr::filter(start > min(interesting_positions), end < max(interesting_positions)) %>%
    mutate(UTR5 = (start-12) - 3939831,
           UTR3 = (end - 3939831))
  return(full_id_table_genes)
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................input data
load(here("data/tidy_data/ecoli_tex_id_table"))
ecoli_tex_data <- full_id_table
load(here("data/tidy_data/ecoli_notex_id_table"))
ecoli_notex_data <- full_id_table

#...................................calculate 5´utr
ecoli_tex_data_mod <- analyse_5utr(id_table = ecoli_tex_data) %>%
  mutate(sample = "Ecoli\n(TEX)")
ecoli_notex_data_mod <- analyse_5utr(id_table = ecoli_notex_data) %>%
  mutate(sample = "Ecoli\n(NOTEX)")

ecoli_tex_data_all <- analyse_rDNA_all(id_table = ecoli_tex_data) %>%
  mutate(sample = "ecoli\n(TEX)") 

ecoli_notex_data_all <- analyse_rDNA_all(id_table = ecoli_notex_data) %>%
  mutate(sample = "ecoli\n(NOTEX)")

#...................................collapsed and filter data
combined_16 <- rbindlist(list(ecoli_tex_data_mod, ecoli_notex_data_mod), fill = T) %>%
  dplyr::filter(locus_name %in% c("16S")) %>%
  mutate(sample = as.factor(sample))

combined_23 <- rbindlist(list(ecoli_tex_data_mod, ecoli_notex_data_mod), fill = T) %>%
  dplyr::filter(locus_name %in% c("23S")) %>%
  mutate(sample = as.factor(sample))

combined_5 <- rbindlist(list(ecoli_tex_data_mod, ecoli_notex_data_mod), fill = T) %>%
  dplyr::filter(locus_name %in% c("5S")) %>%
  mutate(sample = as.factor(sample))

combined_positions_all_ecoli <- rbindlist(list(ecoli_tex_data_all, ecoli_notex_data_all)) %>%
  mutate(sample = as.factor(sample))


#...................................for single loci analysis
combined_16_single <- rbindlist(list(ecoli_tex_data_mod, ecoli_notex_data_mod), fill = T) %>%
  dplyr::filter(gene %in% c("rna-b3756")) %>%
  mutate(sample = as.factor(sample))
combined_23_single <- rbindlist(list(ecoli_tex_data_mod, ecoli_notex_data_mod), fill = T) %>%
  dplyr::filter(gene %in% c("rna-b3758")) %>%
  mutate(sample = as.factor(sample))
combined_5_single <- rbindlist(list(ecoli_tex_data_mod, ecoli_notex_data_mod), fill = T) %>%
  dplyr::filter(gene %in% c("rna-b3759")) %>%
  mutate(sample = as.factor(sample))

#...................................reorder levels
combined_16$sample <-  factor(combined_16$sample, 
                              levels = rev(c("Ecoli\n(TEX)", "Ecoli\n(NOTEX)")))

combined_23$sample <-  factor(combined_23$sample, 
                              levels = rev(c("Ecoli\n(TEX)", "Ecoli\n(NOTEX)")))


combined_16_single$sample <-  factor(combined_16_single$sample, 
                                     levels = rev(c("Ecoli\n(TEX)", "Ecoli\n(NOTEX)")))

combined_23_single$sample <-  factor(combined_23_single$sample, 
                                     levels = rev(c("Ecoli\n(TEX)", "Ecoli\n(NOTEX)")))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOT WITH GGPLOT2
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................set colors
npg2 <- c(pal_npg()(10)[7],
          pal_npg()(10)[4])

#...................................single locus 5´UTR 16S
gg_utr_16S_ecoli <- ggplot(data = combined_positions_all_ecoli, aes(x = UTR5 , fill = sample, y = sample)) +
  geom_density_ridges2(alpha = 0.3, size = 1, scale = 1, color = NA) +
  geom_density_ridges2(alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 200, color = NA) +
  scale_x_continuous(limits = c(-500,100), 
                     breaks = c(-293,-175,-115,-66,0), 
                     labels = c(-293,-175,-115,-66,0),
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("read start position of 16S \n(nt relative to mature rRNA start)") +
  ggtitle("") +
  scale_fill_manual(values = npg2) +
  scale_color_manual(values = npg2) +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 

pdf(here("/figures/16S_utr5_ecoli"), 
    width = 10, height = 5, paper = "special",onefile=FALSE)
gg_utr_16S_ecoli
dev.off()

#...................................single locus 3´UTR 16S
as.numeric(levels(as_factor(combined_16_single$end_gene)))-as.numeric(levels(as_factor(combined_16_single$start_gene)))

pdf(here("/figures/190607_16S_utr3_all_upstream.pdf"), 
    width = 10, height = 5, paper = "special",onefile=FALSE)
ggplot(data = combined_positions_all_ecoli, aes(x = UTR3 , fill = sample, y = sample)) +
  geom_density_ridges2(alpha = 0.3, size = 1, scale = 1, color = NA) +
  geom_density_ridges2(alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 200, color = NA) +
  scale_x_continuous(limits = c(1541-100,1541+100), 
                     breaks = c(1541-100,1541,1541+33,1541+100), 
                     labels = c(-100,0,33,100),
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("5`utr [nt]") +
  ggtitle("") +
  scale_fill_npg() +
  scale_color_npg() +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()


#...................................single locus 5´UTR 23S
as.numeric(levels(as_factor(combined_23_single$start_gene)))-as.numeric(levels(as_factor(combined_16_single$start_gene)))

pdf(here("/figures/190607_23S_utr5_all_upstream.pdf"), 
    width = 10, height = 5, paper = "special",onefile=FALSE)
ggplot(data = combined_positions_all_ecoli, aes(x = UTR5 , fill = sample, y = sample)) +
  geom_density_ridges2(alpha = 0.3, size = 1, scale = 1, color = NA) +
  geom_density_ridges2(alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 100, color = NA) +
  scale_x_continuous(limits = c(1896-50,1896 + 150), 
                     breaks = c(1896 - 50,1896 - 7,1896 -0,1896 + 50,1896 + 150), 
                     labels = c(50,7,0,-50,-150),
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("5`utr [nt]") +
  ggtitle("") +
  scale_fill_npg() +
  scale_color_npg() +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()

#...................................single locus 3´UTR 23S
as.numeric(levels(as_factor(combined_23_single$end_gene)))-as.numeric(levels(as_factor(combined_16_single$start_gene)))


pdf(here("/figures/190607_23S_utr3_all_upstream.pdf"), 
    width = 10, height = 5, paper = "special",onefile=FALSE)
ggplot(data = combined_positions_all_ecoli, aes(x = UTR3 , fill = sample, y = sample)) +
  geom_density_ridges2(alpha = 0.3, size = 1, scale = 1, color = NA) +
  geom_density_ridges2(alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 100, color = NA) +
  scale_x_continuous(limits = c(4799-50,4799+50), 
                     breaks = c(4799-50,4799,4799+7,4799+9,4799+50), 
                     labels = c(-50,0,7,9,50),
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("5`utr [nt]") +
  ggtitle("") +
  scale_fill_npg() +
  scale_color_npg() +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()


#...................................single locus 5´UTR 5S
as.numeric(levels(as_factor(combined_5_single$start_gene)))-as.numeric(levels(as_factor(combined_16_single$start_gene)))


pdf(here("/figures/190607_5S_utr5_all_upstream.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = combined_positions_all_ecoli, aes(x = UTR5 , fill = sample, y = sample)) +
  geom_density_ridges2(alpha = 0.3, size = 1, scale = 1, color = NA) +
  geom_density_ridges2(alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 100, color = NA) +
  scale_x_continuous(limits = c(4892-50,4892+50), 
                     breaks = c(4892-50,4892-3,4892,4892+50), 
                     labels = c(-50,-3,0,50),
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("5`utr [nt]") +
  ggtitle("") +
  scale_fill_npg() +
  scale_color_npg() +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()


#...................................single locus 3´UTR 5S
as.numeric(levels(as_factor(combined_5_single$end_gene)))-as.numeric(levels(as_factor(combined_16_single$start_gene)))

pdf(here("/figures/190607_5S_utr3_all_upstream.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = combined_positions_all_ecoli, aes(x = UTR3 , fill = sample, y = sample)) +
  geom_density_ridges2(alpha = 0.3, size = 1, scale = 1, color = NA) +
  geom_density_ridges2(alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 100, color = NA) +
  scale_x_continuous(limits = c(5011-50,5011+50), 
                     breaks = c(5011-50,5011+0,5011+3,5011+50), 
                     labels = c(-50,0,3,50),
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("5`utr [nt]") +
  ggtitle("") +
  scale_fill_npg() +
  scale_color_npg() +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# COVERAGE CALCULATED PLOT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
interesting_positions <- 3939348:3945138

ecoli_16_coverage_tex   <- fread("/Volumes/Lacie/direct_rna_seq/data/mapped_data/190123_all_forward_sorted_depth.txt") %>%
  dplyr::filter(V2 %in% interesting_positions) %>%
  dplyr::rename(position = V2, coverage = V3) %>%
  mutate(dataset = "TEX")

ecoli_16_coverage_notex <- fread("/Volumes/Lacie/190509_notex4/data/mapped_data/BC4_combined_forward_sorted_depth.txt") %>%
  dplyr::filter(V2 %in% interesting_positions) %>%
  dplyr::rename(position = V2, coverage = V3) %>%
  mutate(dataset = "NOTEX")


rescale_factor <- max(ecoli_16_coverage_tex$coverage)/max(ecoli_16_coverage_notex$coverage)
ecoli_16_coverage_notex$coverage <- ecoli_16_coverage_notex$coverage * rescale_factor


#...................................position specific coverage area plot
pdf(here("/figures/190605_coverage_area_ecoli_all.pdf"), 
    width = 26, height = 6, paper = "special",onefile=FALSE)
ggplot(data = ecoli_16_coverage_tex, aes(x = position, y = coverage)) +
  geom_area(alpha = 0.6) +
  geom_area(data = ecoli_16_coverage_notex,fill = "red1", alpha = 0.6) +
  theme_Publication_white() +
  #geom_vline(data = known_positions, aes(xintercept = position, color = as.factor(position))) +
  #geom_text(data = known_positions,aes(label = feature, y = 10000, x = position, color = as.factor(position)), angle = 90, vjust = 1) +
  theme(panel.grid.major = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_npg() +
  xlab("") +
  ylab("coverage") +
  guides(color = F) 
theme(axis.text.y = element_blank())
dev.off()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# GC CALCULATIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load gc table
dna_gc <- read.table(file = "/Volumes/Lacie/SanDisk_save//DNAprops/scales/gc.tsv",
                     fill = TRUE, quote = "", header = TRUE, sep = "\t")

#...................................load sequence
ecoli_fasta <- readDNAStringSet("/Volumes/Lacie/direct_rna_seq/data/genome_data/e_coli_k12.fasta")
interesting_positions_plus1 <- 3939348:3945139
sequence_rdna <- ecoli_fasta$`U00096.2 Escherichia coli str. K-12 substr. MG1655, complete genome`[interesting_positions_plus1]

gc <- matrix(ncol = 2, nrow = length(interesting_positions))

for (i in 1:length(interesting_positions)){
  gc[i,1] <- interesting_positions[i]
  gc[i,2] <- dna_gc$Parameter[dna_gc$Dinucleotide == as.character(substring(sequence_rdna, i, i + 1))]
}


gc_datatable <- gc %>%
  as_tibble() %>%
  dplyr::rename(position = 1, gc_scale = 2)


heat_color_npg <- c(pal_npg()(10)[4],
                    
                    #pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    #pal_npg()(10)[5],
                    #pal_npg()(10)[1],
                    pal_npg()(10)[8])

pdf(here("/figures/190607_gc_rDNA_ecoli.pdf"), 
    width = 26, height = 2, paper = "special",onefile=FALSE)
ggplot(data = gc_datatable, aes(x = position, y = 1, fill = as.factor(gc_scale))) +
  geom_tile(aes(color = as.factor(gc_scale), size = 1)) +
  scale_fill_manual(values = heat_color_npg) +
  scale_color_manual(values = heat_color_npg) +
  theme_Publication_white() +
  theme(panel.grid.major = element_blank()) +
  xlab("") +
  ylab("gc") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank())
dev.off()
