###########################################################################
###########################################################################
###
### TOMBO OUTPUT ANALYSIS
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

#...................................for adding the sequence
base_annotation <- function(coord_left, coord_right, fasta_file){
  fasta <- readDNAStringSet(fasta_file)
  names(fasta) <- "genome"
  sequence <- as.character(fasta$genome[coord_left:coord_right])
  return(sequence)
}

#...................................calculation % of modified reads based on new model calculated
calculate_fraction_after_tombo <- function(organism = c("pfu", "hvo","ecoli"),
                                           model_type = c("de_novo", "model_testing"), 
                                           intermediate_type = c("primary_pre", "closed_cird", "open_circ", "mature"),
                                           data_folder){
  
  # > set positions for 16S
  if (organism == "pfu"){
    interesting_zone <- 120675:122193
  }else if(organism == "hvo"){
    interesting_zone <- hvo_gff$start[1]:hvo_gff$end[1]
  }else if(organism == "ecoli"){
    interesting_zone <- 223771:225312
  }
  
  # > calculated fractions 
  fraction <- fread(paste(data_folder, model_type, "_", intermediate_type, ".fraction_modified_reads.plus.wig", sep = ""), header = F,fill = F)  %>%
    dplyr::rename(position = 1, mod = 2) %>%
    dplyr::filter(position %in% interesting_zone) %>%
    mutate(group_name = intermediate_type)
  
  # >  return values
  return(fraction)
}

#.................................heatmap plotting
four_color_npg <- c(pal_npg()(10)[3],
                    pal_npg()(10)[4],
                    pal_npg()(10)[2],
                    pal_npg()(10)[1]) 

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................genome data
hvo_gff <- read.gff(here("data/genome_data/hvo.gff")) %>%
  dplyr::filter(type == "rRNA")

hvo_fasta <- readDNAStringSet(here("data/genome_data/hvo.fasta"))
names(hvo_fasta) <- "genome"

#...................................input files: 16S modification positions, fasta file
#.................................HVO
interesting_postions_hvo <- c(1598192+910, 1598192+1352, 1598192+1432, 1598192+1450,1598192+1451)
fasta_file_hvo <- here("data/genome_data/hvo.fasta")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# STATISTICS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................calculate % of probability that a position is modified based on de_novo and model_compare tombo calculations

#.................................HVO
#................................DE NOVO NOTEX
#...............................primary_pre
de_novo_primary_pre_fraction   <- calculate_fraction_after_tombo(organism = "hvo", 
                                                                 model_type = "de_novo", 
                                                                 intermediate_type = "primary_pre",
                                                                 data_folder = paste(here(), "/data/reanalysis_tombo/hvo/wig/", sep = ""))

#...............................closed_circ
de_novo_closed_circ_fraction   <- calculate_fraction_after_tombo(organism = "hvo", 
                                                                 model_type = "de_novo", 
                                                                 intermediate_type = "closed_circ",
                                                                 data_folder = paste(here(), "/data/reanalysis_tombo/hvo/wig/", sep = ""))

#...............................open_circ
de_novo_open_circ_fraction   <- calculate_fraction_after_tombo(organism = "hvo", 
                                                               model_type = "de_novo", 
                                                               intermediate_type = "open_circ",
                                                               data_folder = paste(here(), "/data/reanalysis_tombo/hvo/wig/", sep = ""))

#...............................mature
de_novo_mature_fraction   <- calculate_fraction_after_tombo(organism = "hvo", 
                                                            model_type = "de_novo", 
                                                            intermediate_type = "mature",
                                                            data_folder = paste(here(), "/data/reanalysis_tombo/hvo/wig/", sep = ""))

#...............................combine all de_novo
de_novo_all <- rbind(de_novo_primary_pre_fraction,
                     de_novo_closed_circ_fraction,
                     de_novo_open_circ_fraction,
                     de_novo_mature_fraction)


#................................sample_compare NOTEX VS PRIMARY_PRE
#...............................primary_pre
sample_compare_primary_pre_fraction   <- calculate_fraction_after_tombo(organism = "hvo", 
                                                                        model_type = "sample_compare", 
                                                                        intermediate_type = "primary_pre",
                                                                        data_folder = paste(here(), "/data/reanalysis_tombo/hvo/wig/", sep = ""))

#...............................closed_circ
sample_compare_closed_circ_fraction   <- calculate_fraction_after_tombo(organism = "hvo", 
                                                                        model_type = "sample_compare", 
                                                                        intermediate_type = "closed_circ",
                                                                        data_folder = paste(here(), "/data/reanalysis_tombo/hvo/wig/", sep = ""))

#...............................open_circ
sample_compare_open_circ_fraction   <- calculate_fraction_after_tombo(organism = "hvo", 
                                                                      model_type = "sample_compare", 
                                                                      intermediate_type = "open_circ",
                                                                      data_folder = paste(here(), "/data/reanalysis_tombo/hvo/wig/", sep = ""))

#...............................mature
sample_compare_mature_fraction   <- calculate_fraction_after_tombo(organism = "hvo", 
                                                                   model_type = "sample_compare", 
                                                                   intermediate_type = "mature",
                                                                   data_folder = paste(here(), "/data/reanalysis_tombo/hvo/wig/", sep = ""))

#...............................combine all sample_compare
sample_compare_all <- rbind(sample_compare_primary_pre_fraction,
                            sample_compare_closed_circ_fraction,
                            sample_compare_open_circ_fraction,
                            sample_compare_mature_fraction)

#................................sample_compare NOTEX VS DKSGA MATURE
#...............................primary_pre
sample_compare_dksga_primary_pre_fraction   <- calculate_fraction_after_tombo(organism = "hvo", 
                                                                              model_type = "sample_compare_dksga", 
                                                                              intermediate_type = "primary_pre",
                                                                              data_folder = paste(here(), "/data/reanalysis_tombo/hvo/wig/", sep = ""))

#...............................closed_circ
sample_compare_dksga_closed_circ_fraction   <- calculate_fraction_after_tombo(organism = "hvo", 
                                                                              model_type = "sample_compare_dksga", 
                                                                              intermediate_type = "closed_circ",
                                                                              data_folder = paste(here(), "/data/reanalysis_tombo/hvo/wig/", sep = ""))

#...............................open_circ
sample_compare_dksga_open_circ_fraction   <- calculate_fraction_after_tombo(organism = "hvo", 
                                                                            model_type = "sample_compare_dksga", 
                                                                            intermediate_type = "open_circ",
                                                                            data_folder = paste(here(), "/data/reanalysis_tombo/hvo/wig/", sep = ""))

#...............................mature
sample_compare_dksga_mature_fraction   <- calculate_fraction_after_tombo(organism = "hvo", 
                                                                         model_type = "sample_compare_dksga", 
                                                                         intermediate_type = "mature",
                                                                         data_folder = paste(here(), "/data/reanalysis_tombo/hvo/wig/", sep = ""))

#...............................combine all sample_compare
sample_compare_dksga_all <- rbind(sample_compare_dksga_primary_pre_fraction,
                                  sample_compare_dksga_closed_circ_fraction,
                                  sample_compare_dksga_open_circ_fraction,
                                  sample_compare_dksga_mature_fraction)


#...............................combine sample compare and de novo model
all_models <- rbind(notex_and_dksga_de_novo %>% mutate(model = "de_novo"),
                    sample_compare_all %>% mutate(model = "sample_compare", seqset = "wt"),
                    sample_compare_dksga_all%>% mutate(model = "sample_compare_dksga", seqset = "wt"))

# > annotation seq
left_b  <- -17
right_b <- +11
annotation_seq <- base_annotation(coord_left = interesting_postions_hvo[3]+left_b, 
                                  coord_right = interesting_postions_hvo[3]+right_b,
                                  fasta_file = here("data/genome_data/hvo.fasta"))

all_models$group_name <- factor(all_models$group_name, levels = rev(c("primary_pre", "closed_circ", "open_circ", "mature")))

all_models_simple <- all_models %>%
  dplyr::filter(model == "sample_compare",seqset == "wt", group_name != "primary_pre")

#............................simply de novo m62A hvo wt
pdf(here("figures/m6a_hvo_sample_compare_ccg.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggplot(data = all_models_simple, aes(x = position, y = mod, color = group_name)) +
  facet_grid(rows = vars(seqset), cols = vars(model)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_line(size = 2, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0),
                     limits = c((interesting_postions_hvo[3]+left_b),(interesting_postions_hvo[3]+right_b)), 
                     breaks = c((interesting_postions_hvo[3]+left_b):(interesting_postions_hvo[3]+right_b)), 
                     labels = strsplit(annotation_seq, "*")[[1]]) +
  scale_color_manual(values = pal_npg()(10)[c(4,8,2,7)]) +
  scale_fill_manual(values = pal_npg()(10)[c(4,8,2,7)]) +
  theme_Publication_white() 
dev.off()



#............................simply de novo m62A hvo wt
all_models_simple <- all_models %>%
  dplyr::filter(seqset == "wt", group_name == "mature")

pdf(here("figures/m6a_hvo_samplecompare_ksga.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggplot(data = all_models_simple, aes(x = position, y = mod, color = model)) +
  facet_grid(rows = vars(seqset)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_line(size = 2, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0),
                     limits = c((interesting_postions_hvo[4]+left_b),(interesting_postions_hvo[4]+right_b)), 
                     breaks = c((interesting_postions_hvo[4]+left_b):(interesting_postions_hvo[4]+right_b)), 
                     labels = strsplit(annotation_seq, "*")[[1]]) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_color_manual(values = pal_npg()(10)[c(4,8,2,7)]) +
  scale_fill_manual(values = pal_npg()(10)[c(4,8,2,7)]) +
  theme_Publication_white() 
dev.off()


