###########################################################################
###########################################################################
###
### FRACTION OF FULL-LENGTH TRANSCRIPTS
###
###########################################################################
###########################################################################


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES AND PLOTTING FUNCTION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(here)
source(here("Rscripts/load_libraries.R"))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................tidy table to gene_length groups
tidy_to_fraction <- function(id_table_input, seq_name){
  return(id_table_input) %>%
  dplyr::filter(mapped_type == "CDS") %>%
    mutate(length_gene = abs(start_gene - end_gene),
           fraction_read = (aligned_reads/length_gene*100)) %>%
    mutate(template_length = end_gene - start_gene,
           group = ifelse(template_length < 500, "[0-500]",
                          ifelse(template_length < 1000 & template_length >= 500, "[500-1000]",
                                 ifelse(template_length < 1500 & template_length >= 1000, "[1000-1500]", 
                                        ifelse(template_length < 2000 & template_length >= 1500, "[1500-2000]", 
                                               ifelse(template_length >= 2000, "[>=2000]", NA)))))) %>%
    dplyr::select(fraction_read, group) %>%
    mutate(sequencing_set = seq_name)
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load data

#.................................ecoli
load(here("/data/tidy_data/ecoli_tex_id_table"))
ecoli_id <- tidy_to_fraction(full_id_table, "Ecoli (TEX)")
  
#.................................pfu
load(here("/data/tidy_data/pfu_tex_id_table"))
pfu_id <- tidy_to_fraction(full_id_table, "Pfu (TEX)")

#.................................hvo
load(here("/data/tidy_data/hvo_tex_id_table"))
hvo_id <- tidy_to_fraction(full_id_table, "Hvo (TEX)")

#...................................combine all data
all_id <- bind_rows(ecoli_id, pfu_id, hvo_id)

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

