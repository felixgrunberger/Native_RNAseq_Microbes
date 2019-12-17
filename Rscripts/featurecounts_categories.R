###########################################################################
###########################################################################
###
### COUNTS IN CATEGORIES (READS MAPPING TO WHICH FEATURES?)
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
# only for wt tex vs notex
gene_files <- gene_files[-c(3,5)]

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
  ylab("Log10 Mapped reads") +
  ggtitle("") +
  labs(fill = "") +
  coord_flip() +
  guides(color = FALSE, fill = guide_legend("")) +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(panel.grid.major.y = element_blank()) +
  scale_y_log10(expand = c(0,0, 0.1, 0), breaks = c(0, 100, 1000, 10000, 100000, 1000000))

pdf(here("figures/mapping_to_features_total_counts.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
gg_total_counts
dev.off()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# WRITE TO TABLES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................ONT transcript abundance data for TEX samples 
counts_files_tex <- paste(here("data/tidy_data/"), list.files(here("data/tidy_data/"),pattern = "_tex_gene_table"), sep = "")

#...................................get sample names
sample_names <- unlist(lapply(counts_files_tex, FUN=function(x){str_split_fixed(str_split_fixed(x, "_gene", 2)[1],"tidy_data/",2)[2]}))

#...................................get annotation files
ecoli_annotation <- fread(here("data/genome_data/ecoli_annotation.tsv"))
pfu_annotation <- fread(here("data/genome_data/pfu_annotation.tsv"))
hvo_annotation <- fread(here("data/genome_data/hvo_annotation.tsv"))

#...................................write to tsv files
#.................................ecoli
load(counts_files_tex[1])
full_gene_table_ecoli <- full_gene_table %>%
    left_join(ecoli_annotation, by = c("GeneID" = "id")) %>%
    dplyr::select(GeneID, Chr, counts, locus_name, locus_tag) %>%
  arrange(desc(counts))

writexl::write_xlsx(x = full_gene_table_ecoli, 
                    path = here("tables/counts_tables/counts_ecoli_tex.xlsx"))

#.................................pfu
load(counts_files_tex[3])
full_gene_table_pfu <- full_gene_table %>%
  rowwise() %>%
  mutate(GeneID = str_split_fixed(GeneID, ".p01", 2)[1]) %>%
  left_join(pfu_annotation, by = c("GeneID" = "gene")) %>%
  dplyr::rename(locus_tag = old_name) %>%
  dplyr::select(GeneID, Chr, counts, locus_name, locus_tag) %>%
  arrange(desc(counts))

writexl::write_xlsx(x = full_gene_table_pfu, 
                    path = here("tables/counts_tables/counts_pfu_tex.xlsx"))

#.................................hvo
load(counts_files_tex[2])
full_gene_table_hvo <- full_gene_table %>%
  left_join(hvo_annotation, by = c("GeneID" = "id_name")) %>%
  dplyr::rename(locus_tag = old_name) %>%
  dplyr::select(GeneID, Chr, counts, locus_name, locus_tag) %>%
  arrange(desc(counts))

writexl::write_xlsx(x = full_gene_table_hvo, 
                    path = here("tables/counts_tables/counts_hvo_tex.xlsx"))





