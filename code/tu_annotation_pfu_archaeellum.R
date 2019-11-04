###########################################################################
###########################################################################
###
### PLOT FLA-OPERON IN P. FURIOSUS
###
###########################################################################
###########################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
packages <- c("tidyverse", "here", "ggthemes", "data.table", 
              "ggExtra", "Rsamtools", "GenomicAlignments", "seqTools", 
              "Rsubread", "ape", "DT", "ggpubr", "ggridges", "ggsci", 
              "GenomicRanges","IRanges")
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

#.................................load samtools depth data
pfu_coverage_plus  <- fread(here("data/coverage_data/pfu_tex_plus_depth.txt.gz"))
pfu_coverage_minus <- fread(here("data/coverage_data/pfu_tex_minus_depth.txt.gz"))

#.................................filter for interesting region
fla_positions <- c(324688:337796)

#...............................NANOPORE READ DEPTH
pfu_coverage_plus_fla <- pfu_coverage_plus %>%
  dplyr::rename(position = 2, depth = 3) %>%
  dplyr::filter(position %in% fla_positions)

pfu_coverage_minus_fla <- pfu_coverage_minus %>%
  dplyr::rename(position = 2, depth = 3) %>%
  dplyr::filter(position %in% fla_positions)

#...............................NANOPORE DETECTED OPERONS
pfu_sub_operons_annotated_fla <- read_tsv(here("data/operon_data/pfu_tex_suboperons.tsv")) %>%
  dplyr::filter(start_operon > min(fla_positions), end_operon < max(fla_positions))

#...............................CURRENT ANNOTATION
gff_cds_fla <- read.gff(here("data/genome_data/pfu.gff")) %>%
  dplyr::filter(start > min(fla_positions), end < max(fla_positions))

#...............................DOOR2 OPERONS
door <- fread(here("data/operon_data/pfu_database_operons.tsv"))

# --> region spanning from 01600 to 1660 
door_table <- gff_cds_fla 
door_table$start <- min(door_table$start)
door_table$end <- max(door_table$end)
door_table <- door_table %>%
  head(n= 1)

#.................................plot coverage of Nanopore reads and annotation (compare Fig. 3)
area_plot <- ggplot() +
  geom_area(data = pfu_coverage_plus_fla, aes(x = position, y = -depth), size = 3,color = "black", fill = "black") +
  geom_area(data = pfu_coverage_minus_fla, aes(x = position, y = depth), size = 3,color = "black", fill = "black") +
  theme_Publication_white() +
  theme(panel.grid.major = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(min(fla_positions), max(fla_positions)), expand = c(0,0)) +
  scale_color_npg() +
  xlab("") +
  ylab("Coverage") +
  guides(color = F) 

annotation_plot <- ggplot() +
  geom_rect(data = pfu_sub_operons_annotated_fla, aes(xmin = start_operon, xmax = end_operon, ymin = 0, ymax = 1), color = "black", fill = "black", alpha = 0.5) +
  geom_rect(data = gff_cds_fla, aes(xmin = start, xmax = end, ymin = -1, ymax = 0), color = "black", fill = "grey", alpha = 0.5) +
  geom_rect(data = door_table, aes(xmin = start, xmax = end, ymin = -2, ymax = -1), color = "red", fill = "red", alpha = 0.5) +
  theme_Publication_white() +
  theme(panel.grid.major = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(min(fla_positions), max(fla_positions)), expand = c(0,0)) +
  scale_color_npg() +
  xlab("") +
  ylab("") +
  theme(axis.text.y = element_blank()) +
  guides(color = F) 

pdf(here("/figures/fla_operon_pfu"), 
    width = 23, height = 6, paper = "special",onefile=FALSE)
ggarrange(area_plot, annotation_plot, nrow = 2, heights = c(80,20))
dev.off()

