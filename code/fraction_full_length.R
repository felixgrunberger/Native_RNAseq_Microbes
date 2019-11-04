###########################################################################
###########################################################################
###
### FRACTION OF FULL-LENGTH TRANSCRIPTS
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


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load data
#.................................ecoli

#..............................load gff annotation
interesting_list <- c("CDS", "rRNA", "tRNA")
ecoli_gff <- here("data/genome_data/ecoli.gff")

gff_table_ecoli <- read.gff(ecoli_gff) %>%
  as_tibble() %>%
  mutate(start_gene = start, end_gene = end,strand_gene = strand) %>%
  dplyr::filter(type %in% interesting_list) %>%
  mutate(id_name = str_split_fixed(str_split_fixed(attributes, ";Parent=", 2)[,1], "ID=", 2)[,2],
         locus_name = ifelse(type == "CDS", str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], ";", 2)[,1],
                             ifelse(type == "rRNA", str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], " ", 2)[,1], 
                                    ifelse(type == "tRNA", str_split_fixed(attributes, ";product=", 2)[,2], NA ))),
         parent = ifelse(type == "CDS", str_split_fixed(str_split_fixed(str_split_fixed(attributes, ";Dbxref=", 2)[,1], ";", 2)[,2],"Parent=gene-",2)[,2], NA),
         gene_name = ifelse(type == "rRNA", str_split_fixed(str_split_fixed(attributes, "ribosomal RNA of ", 2)[,2], " operon", 2)[,1], NA),
         length_gene = end_gene - start_gene) %>%
  dplyr::select(start_gene, end_gene, strand_gene, parent, locus_name, id_name, type, length_gene)

#..............................load tu table
load(here("/data/tidy_data/ecoli_tex_id_table"))

ecoli_id <- full_id_table %>%
  dplyr::filter(mapped_type == "CDS") %>%
  mutate(length_gene = abs(start_gene - end_gene),
         fraction_read = (aligned_reads/length_gene*100)) %>%
  mutate(template_length = end_gene - start_gene,
         group = ifelse(template_length < 500, "[0-500]",
                        ifelse(template_length < 1000 & template_length >= 500, "[500-1000]",
                               ifelse(template_length < 1500 & template_length >= 1000, "[1000-1500]", 
                                      ifelse(template_length < 2000 & template_length >= 1500, "[1500-2000]", 
                                             ifelse(template_length >= 2000, "[>=2000]", NA)))))) %>%
  dplyr::select(fraction_read, group)

#.................................hvo

#..............................load gff annotation
hvo_gff <- here("data/genome_data/hvo.gff")

gff_table_hvo <- read.gff(hvo_gff) %>%
  as_tibble() %>%
  dplyr::filter(type %in% interesting_list) %>%
  mutate(id_name = str_split_fixed(str_split_fixed(attributes, ";Parent=", 2)[,1], "ID=", 2)[,2],
         locus_name = ifelse(type == "CDS", str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], ";", 2)[,1],
                             ifelse(type == "rRNA", str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], " ", 2)[,1], 
                                    ifelse(type == "tRNA", str_split_fixed(attributes, ";product=", 2)[,2], NA ))),
         parent = ifelse(type == "CDS", str_split_fixed(str_split_fixed(str_split_fixed(attributes, ";Dbxref=", 2)[,1], ";", 2)[,2],"Parent=gene-",2)[,2], NA),
         gene_name = ifelse(type == "rRNA", str_split_fixed(str_split_fixed(attributes, "ribosomal RNA of ", 2)[,2], " operon", 2)[,1], NA),
         geneID = ifelse(type == "rRNA", str_split_fixed(str_split_fixed(attributes, "GeneID:", 2)[,2],";gbkey", 2)[,1], NA)) %>%
  mutate(strand_gene = strand,
         start_gene = start,
         end_gene = end,
         length_gene = end_gene - start_gene) %>%
  dplyr::select(start_gene, end_gene, strand_gene, id_name,locus_name, gene_name, geneID, length_gene)

#..............................load tu table
load(here("/data/tidy_data/hvo_tex_id_table"))

hvo_id <- full_id_table %>%
  dplyr::filter(mapped_type == "CDS") %>%
  mutate(length_gene = abs(start_gene - end_gene),
         fraction_read = (aligned_reads/length_gene*100)) %>%
  mutate(template_length = end_gene - start_gene,
         group = ifelse(template_length < 500, "[0-500]",
                        ifelse(template_length < 1000 & template_length >= 500, "[500-1000]",
                               ifelse(template_length < 1500 & template_length >= 1000, "[1000-1500]", 
                                      ifelse(template_length < 2000 & template_length >= 1500, "[1500-2000]", 
                                             ifelse(template_length >= 2000, "[>=2000]", NA)))))) %>%
  dplyr::select(fraction_read, group)

#.................................pfu

#..............................load gff annotation
pfu_gff <- here("data/genome_data/pfu.gff")

gff_table_pfu <- read.gff(pfu_gff) %>%
  dplyr::filter(type == "CDS") %>%
  mutate(start_gene = start, end_gene = end, strand_gene = strand,
         id_gene = str_split_fixed(str_split_fixed(attributes, ".p",2)[,1], "ID=",2)[,2],
         length_gene = end_gene - start_gene) %>%
  dplyr::select(start_gene, end_gene, strand_gene, id_gene,length_gene)

#..............................load tu table
load(here("/data/tidy_data/pfu_tex_id_table"))

pfu_id <- full_id_table %>%
  dplyr::filter(mapped_type == "CDS") %>%
  mutate(length_gene = abs(start_gene - end_gene),
         fraction_read = (aligned_reads/length_gene*100)) %>%
  mutate(template_length = end_gene - start_gene,
         group = ifelse(template_length < 500, "[0-500]",
                        ifelse(template_length < 1000 & template_length >= 500, "[500-1000]",
                               ifelse(template_length < 1500 & template_length >= 1000, "[1000-1500]", 
                                      ifelse(template_length < 2000 & template_length >= 1500, "[1500-2000]", 
                                             ifelse(template_length >= 2000, "[>=2000]", NA)))))) %>%
  dplyr::select(fraction_read, group)

#...................................combine all data
all_id <- bind_rows(ecoli_id %>%
                      mutate(sequencing_set = "Ecoli (TEX)"), 
                    hvo_id %>% 
                      mutate(sequencing_set = "Hvo (TEX)"), 
                    pfu_id %>%
                      mutate(sequencing_set = "Pfu (TEX)"))

#.................................reorder Caccording to protein coding gene lengths
all_id$group <-  factor(all_id$group, 
                        levels = c("[0-500]","[500-1000]","[1000-1500]","[1500-2000]","[>=2000]"))

#.................................plot fraction of full-length transcripts (Supplementary Fig. 8a)
fraction_percentage_plot <- ggplot(data = all_id, aes(x = group, y = fraction_read,  fill = sequencing_set)) +
  geom_violin(alpha = 0.15,position="identity", size = 0.9, color = NA) + 
  theme_Publication_white() +
  xlab("CDS length (nt)") +
  ylab("Fraction of full-length \n(%, log10 scale)") +
  scale_y_log10(expand = c(0,0), breaks = c(10,25,50,100,1000)) +
  scale_fill_npg() +
  guides(fill = guide_legend(title = ""), color = guide_legend(title = "")) +
  geom_boxplot(alpha = 1, size = 0.7, width = 0.4, color = "black", position = position_dodge(width = 0.5), outlier.alpha = 0.1) 

pdf(here("figures/fraction_full_length.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
fraction_percentage_plot
dev.off()

