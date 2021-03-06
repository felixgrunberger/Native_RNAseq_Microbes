###########################################################################
###########################################################################
###
### READ IDENTITY [%] TO GENOMIC FEATURES AND ENOLASE SPIKE-IN 
###
###########################################################################
###########################################################################

# (identity = (1 - NM/aligned_reads)*100)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES AND PLOTTING FUNCTION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(here)
source(here("Rscripts/load_libraries.R"))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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
# LOAD & TIDY DATA GENOME
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................input: pre-calculated single-read tidy files
id_files <- paste(here("data/tidy_data/"), list.files(here("data/tidy_data/"),pattern = "_id_table"), sep = "")

#...................................calculate genome tables
genome_table <- data.frame()

for (i in seq_along(sample_names)){
  load(id_files[i])
  genome_table <- rbindlist(list(genome_table, modify_id_table(full_id_table, sample_names[i])))
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MERGE ALL DATA 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................combine enolase and genome tables
full_table <- rbindlist(list(enolase_table, genome_table), fill = TRUE)

#...................................reorder levels
full_table$sequencing_set <-  factor(full_table$sequencing_set, 
                                     levels = rev(c("ecoli_tex","ecoli_notex",
                                                    "pfu_tex", "pfu_notex",
                                                    "hvo_tex", "hvo_notex")))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CALCULATE STATISTICS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................calculate mapped read statistics
#.................................number of mapped reads / median identity
full_table %>%
  dplyr::filter(mapped_to == "genome") %>%
  group_by(sequencing_set) %>%
  summarise(number_of_mapped_reads = n(),
            median_identity = median(identity))

#.................................CDS mapping features
#...............................calculate bases on CDS
cds_pfu   <- cds_depth_calculator(here("data/genome_data/pfu.gff"))
cds_hvo   <- cds_depth_calculator(here("data/genome_data/hvo.gff"))
cds_ecoli <- cds_depth_calculator(here("data/genome_data/ecoli.gff"))

full_table %>%
  mutate(cds_depth = ifelse(sequencing_set == "ecoli_notex" | sequencing_set == "ecoli_tex", cds_ecoli,
                            ifelse(sequencing_set == "pfu_notex" | sequencing_set == "pfu_tex", cds_pfu, cds_hvo))) %>%
  group_by(sequencing_set) %>%
  dplyr::filter(mapped_type == "CDS", mapped_to == "genome") %>%
  summarise(reads_mapped_to_CDS = n(),
            bases_mapped_to_CDS = sum(aligned_reads),
            sequencing_depth_CDS = sum(aligned_reads)/max(cds_depth))

#.................................rRNA mapping features
full_table %>%
  group_by(sequencing_set) %>%
  dplyr::filter(mapped_type == "rRNA", mapped_to == "genome") %>%
  summarise(reads_mapped_to_rRNA = n())

#.................................enolase mapping features
full_table %>%
  group_by(sequencing_set) %>%
  dplyr::filter(mapped_type == "CDS", mapped_to == "control") %>%
  summarise(reads_mapped_to_enolase = n(),
            median_length_enolase = median(as.numeric(sequence_length_template)),
            median_quality_enolase = median(as.numeric(mean_qscore_template)))



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................change colors
heat_color_npg <- rev(c(pal_npg()(10)[4],
                        pal_npg()(10)[1]))


#...................................for different group layout
full_table_rearranged <- full_table %>%
  ungroup() %>%
  mutate(new_group = ifelse(mapped_to == "control", "control",
                            ifelse(mapped_to == "genome" & mapped_type == "CDS", "CDS",
                                   ifelse(mapped_to == "genome" & mapped_type == "rest", "rest",
                                          ifelse(mapped_to == "genome" & mapped_type == "rRNA", "rRNA", NA)))))

#...................................plot mapped read identity comparing genome-derived reads with enolase reads (Supplementary Fig. 3b)
gg_identity <- ggplot(data = full_table_rearranged, aes(x = identity, color = mapped_to, fill = mapped_to, linetype = mapped_type, y = sequencing_set, alpha = mapped_to)) +
  geom_density_ridges2(size = 1, scale = 0.9) +
  scale_alpha_discrete(range = c(0.5, 0.2)) +
  theme_Publication_white() +
  scale_x_continuous(limits = c(65,100), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0, 0.2, 0)) +
  ylab("") +
  xlab("Mapped read identity (%)") +
  scale_linetype_manual(values = c("CDS" = "solid", "tRNA" = "dashed", "rRNA" = "dotted")) +
  scale_color_manual(values = heat_color_npg) +
  scale_fill_manual(values = heat_color_npg) +
  theme(axis.text.y = element_text(face = "italic")) +
  guides(linetype = guide_legend(title="")) +
  guides(fill = guide_legend(title="")) +
  guides(color = guide_legend(title=""))

pdf(here("figures/identity_control_vs_genome.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_identity
dev.off()


#...................................plot mapped read lengths comparing genome-derived reads with enolase reads (Supplementary Fig. 3c)
gg_length <- ggplot(data = full_table, aes(x = aligned_reads, color = mapped_to, fill = mapped_to, linetype = mapped_type, y = sequencing_set, alpha = mapped_to)) +
  geom_density_ridges2(aes(height =..ndensity..), scale = 0.9, size = 1) +
  theme_Publication_white() +
  scale_alpha_discrete(range = c(0.5, 0.2)) +
  scale_x_continuous(trans = "log10", limits = c(50,5000), breaks = c(100,1000,1314,3000,5000),expand = c(0, 0)) +
  scale_y_discrete(expand = c(0,0, 0.2, 0)) +
  ylab("") +
  xlab("Log10 Mapped read length (nt)") +
  scale_linetype_manual(values = c("CDS" = "solid", "tRNA" = "dashed", "rRNA" = "dotted")) +
  scale_color_manual(values = heat_color_npg) +
  scale_fill_manual(values = heat_color_npg) +
  theme(axis.text.y = element_text(face = "italic")) +
  guides(linetype = guide_legend(title="")) +
  guides(fill = guide_legend(title="")) +
  guides(color = guide_legend(title=""))

pdf(here("figures/aligned_length_control_vs_genome.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_length
dev.off()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MINIMUM MAPPED READ LENGHTS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................set heatmap like colors
heat_color_npg <- c(pal_npg()(10)[4],
                    pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    pal_npg()(10)[5],
                    pal_npg()(10)[1],
                    pal_npg()(10)[8])

#...................................filter genome table for short reads
genome_table_short <- genome_table %>%
  dplyr::filter(sequence_length_template < 300)

#...................................calculate frequencies for short reads
calculate_frequency <- function(set_org){
  
  # > filter for organism
  chosen_set <- genome_table %>%
    dplyr::filter(sequencing_set == set_org) %>%
    dplyr::filter(sequence_length_template < 300)
  
  # > list read lengths
  list_lengths <- NULL
  list_lengths <- list(as.numeric(chosen_set$sequence_length_template))
  
  # > make data table 
  heat_data <- as.data.frame(table(unlist(list_lengths))) %>%
    mutate(coordinate = as.numeric(as.character(Var1))) %>%
    mutate(group = set_org) %>% 
    mutate(Freq = Freq/sum(Freq))
}


#...................................calculate heat tables for frequency of short redas
short_read_table <- data.frame()

for (i in seq_along(sample_names)){
  short_read_table <- rbindlist(list(short_read_table, calculate_frequency(sample_names[i])))
}

#...................................reorder levels
short_read_table$group <-  factor(short_read_table$group, 
                                     levels = rev(c("ecoli_tex","ecoli_notex",
                                                    "pfu_tex", "pfu_notex",
                                                    "hvo_tex", "hvo_notex")))

#...................................plot frequencies of short reads (Supplementary Fig. 3d)
gg_heatmap_shortreads <- ggplot() +
  geom_tile(data = short_read_table, aes(y = group, x = coordinate, color = Freq, fill= Freq), size = 0.5) +
  scale_x_continuous(limits = c(0,300),expand = c(0,0)) +
  scale_fill_gradientn(colours = heat_color_npg) +
  scale_color_gradientn(colours = heat_color_npg) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("Minimum read length (nt)") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "white") +
  guides(fill = guide_colorbar(title = "counts",barwidth = 15, barheight = 0.5, ticks = T, label = T)) +
  guides(color = F)

pdf(here("figures/short_reads_heatmap.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_heatmap_shortreads
dev.off()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# EXAMPLE PLOTS FOR PYROCOCCUS FURIOSUS SET
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................extract Pyrococcus furiosus TEX data (only CDS)
pfu_genome_table <- genome_table %>%
  dplyr::filter(sequencing_set == "pfu_tex", mapped_to == "genome", mapped_type == "CDS") 



#...................................Mapped identity of reads vs. mapped read length (including spearman correlation, Supplementary Fig. 3e)
identity_length <- ggplot(data = pfu_genome_table, aes(x = identity, y = as.numeric(aligned_reads))) +
  xlab("Mapped identity (%)") +
  ylab("Mapped read length (nt)") +
  ggtitle("") +
  stat_binhex(bins = 30, aes(fill = ..count.., alpha = ..count..), color = "white") + 
  scale_fill_gradientn(colours = heat_color_npg) +
  scale_alpha(range = c(0.7,1)) +
  theme_Publication_white() +
  guides(alpha = F) +
  guides(fill = guide_colorbar(title = "counts",barwidth = 15, barheight = 0.5, ticks = T, label = T)) +
  stat_cor(method = "spearman", color = "black")

pdf(here("figures/identity_vs_length.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
identity_length
dev.off()

#...................................Mapped identity of reads vs. read quality (including spearman correlation, Supplementary Fig. 3f)
identity_quality <- ggplot(data = pfu_genome_table, aes(x = identity, y = as.numeric(mean_qscore_template))) +
  xlab("Mapped identity (%)") +
  ylab("Read quality (Phred-like score)") +
  ggtitle("") +
  stat_binhex(bins = 30, aes(fill = ..count.., alpha = ..count..), color = "white") + 
    scale_fill_gradientn(colours = heat_color_npg) +
  scale_alpha(range = c(0.7,1)) +
  theme_Publication_white() +
  guides(alpha = F) +
  guides(fill = guide_colorbar(title = "counts",barwidth = 15, barheight = 0.5, ticks = T, label = T)) +
  stat_cor(method = "spearman", color = "black")

pdf(here("figures/identity_vs_quality.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
identity_quality
dev.off()




