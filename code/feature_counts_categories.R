###########################################################################
###########################################################################
###
### COUNTS IN CATEGORIES (READS MAPPING TO WHICH FEATURES?)
###
###########################################################################
###########################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
packages <- c("tidyverse", "here", "ggthemes", "data.table", 
              "ggExtra", "Rsamtools", "GenomicAlignments", 
              "seqTools", "Rsubread", "ape", "DT", "ggpubr", "ggridges", "ggsci")
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

#...................................calculate number of reads mapping to features protein coding genes (CDS) and rRNA
category_calculator <- function(gene_table, gff, identifier){

  # > calculate for CDS and rRNA 
  interesting_list <- c("CDS", "rRNA")
  
  # > read in gff file and modify gene names
  gff_table <- read.gff(gff) %>%
    as_tibble() %>%
    dplyr::filter(type %in% interesting_list) %>%
    mutate(id_name = str_split_fixed(str_split_fixed(attributes, ";Parent=", 2)[,1], "ID=", 2)[,2],
           locus_name = ifelse(type == "CDS", str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], ";", 2)[,1],
                               ifelse(type == "rRNA", str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], " ", 2)[,1], 
                                      ifelse(type == "tRNA", str_split_fixed(attributes, ";product=", 2)[,2], NA )))) %>%
    dplyr::select(id_name, locus_name)
  
  # > read in output from tidy_data code, filter out tRNAs and calculate number of reads mapping to a group
  dataset <- gene_table %>%
    dplyr::filter(type != "tRNA") %>%
    mutate(split_group = ifelse(type == "rRNA", locus_name, type)) %>%
    group_by(split_group) %>%
    summarise(total_count = sum(counts)) %>%
    mutate(percentage = round(total_count/(sum(total_count))*100, digits = 1),
           sequencing_set = identifier)
  
  return(dataset)
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................paths to _gene_tables from tidy_data input
gene_files <- paste(here("data/tidy_data/"), list.files(here("data/tidy_data/"),pattern = "_gene_table"), sep = "")

#...................................get sample names
sample_names <- unlist(lapply(gene_files, FUN=function(x){str_split_fixed(str_split_fixed(x, "_gene", 2)[1],"tidy_data/",2)[2]}))

#...................................create data.frame
gene_table_counts <- data.frame()

for (i in seq_along(sample_names)){
  
  # > set sample fasta, gff and input gene table file with counts information
  working_directory   <- paste(here(), "/data", sep = "")
  sample_name         <- sample_names[i]
  input_fasta         <- paste(working_directory, "/genome_data/", str_split_fixed(sample_name, "_",2)[1],".fasta", sep = "")
  input_gff           <- paste(working_directory, "/genome_data/", str_split_fixed(sample_name, "_",2)[1],".gff", sep = "")
  input_gene_table    <- gene_files[i]
  
  # > load saved R file
  load(input_gene_table)
  
  # > add gff information to count tables
  gene_table_counts <- rbind(gene_table_counts, category_calculator(full_gene_table, input_gff, sample_name))

}

#...................................reorder levels
gene_table_counts$sequencing_set <-  factor(gene_table_counts$sequencing_set, 
                                            levels = rev(c("ecoli_tex","ecoli_notex",
                                                           "pfu_tex", "pfu_notex",
                                                           "hvo_tex", "hvo_notex")))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOT WITH GGPLOT2
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#...................................change colors
heat_color_npg <- c(pal_npg()(10)[4],
                    pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    pal_npg()(10)[1])

#...................................% mapped to features CDS and 5S rRNA, 16S rRNA, 23S rRNA (Fig. 1b)
gg_percentage <- ggplot(data = gene_table_counts, aes(x = sequencing_set, y = percentage, fill = split_group)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = heat_color_npg) +
  theme_Publication_white() +
  xlab("") +
  ylab("Mapped reads to feature (%)") +
  ggtitle("") +
  labs(fill = "") +
  coord_flip() +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(panel.grid.major.y = element_blank()) +
  scale_y_continuous(expand = c(0, 0))

pdf(here("figures/mapping_to_features_percentage.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
gg_percentage
dev.off()

#...................................total_count mapped to features CDS and 5S rRNA, 16S rRNA, 23S rRNA (Supplementary Fig. 3a)
gg_total_counts <- ggplot(data = gene_table_counts, aes(x = sequencing_set, y = total_count, fill = split_group, color = split_group)) +
  geom_bar(stat="identity", 
           position=position_dodge(width=0.5), width=0.1) +
  geom_point(position=position_dodge(width = 0.5), 
             mapping = aes(group = split_group), size = 3, fill = "white", shape = 21, stroke = 2) +
  scale_fill_manual(values = heat_color_npg) +
  scale_color_manual(values = heat_color_npg) +
  theme_Publication_white() +
  xlab("") +
  ylab("Mapped reads to feature [M counts]") +
  ggtitle("") +
  labs(fill = "") +
  coord_flip() +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(panel.grid.major.y = element_blank()) +
  scale_y_log10(expand = c(0,0, 0.1, 0), breaks = c(0, 100, 1000, 10000, 100000, 1000000))

pdf(here("figures/mapping_to_features_total_counts.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
gg_total_counts
dev.off()


