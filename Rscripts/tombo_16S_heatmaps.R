###########################################################################
###########################################################################
###
### TOMBO OUTPUT HEATMAP ANALYSIS
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

#...................................calculate ggplot2 datatable of known modifications site
calc_mod_heatmap <- function(inputset_positions, tombo_model){
  nl <- NA
  for (i in seq_along(inputset_positions)){
    nl <- c(nl, seq(inputset_positions[i]-1, inputset_positions[i]+1, by = 1))
    nl <- nl[!is.na(nl)]
  }
  nl2 <- tombo_model %>%
    mutate(mod = ifelse(position %in% nl, 1, 0))
  return(nl2)
}

#...................................plot 16S heatmaps
plot_16S_heatmap <- function(input_model, input_denovo, input_known_positions){
  ggplot(data = input_model, aes(x = position, y = "tombo model", fill = mod, color = mod)) +
    geom_tile(size = 1) +
    geom_tile(data = input_denovo, aes(x = position, y = "tombo de novo",  fill = mod, color = mod), size = 1) +
    geom_tile(data = input_known_positions, aes(x = position, y = "literature",  fill = mod, color = mod), size = 1) +
    scale_fill_gradientn(colors = heat_color_npg) +
    scale_color_gradientn(colors = heat_color_npg) +
    scale_y_discrete(expand = c(0,0)) +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank()) +
    xlab("") +
    guides(fill = F, color = F)
}

#...................................calculation % of modified reads based on new model calculated
calculate_fraction_after_tombo <- function(organism = c("pfu", "hvo","ecoli"),model_type = c("de_novo", "model_testing"), data_folder){
  
  # > set positions for 16S
  if (organism == "pfu"){
    interesting_zone <- 120675:122193
  }else if(organism == "hvo"){
    interesting_zone <- 1599672:1598192
  }else if(organism == "ecoli"){
    interesting_zone <- 223771:225312
  }
  
  # > calculated fractions 
  fraction <- fread(paste(data_folder, model_type,".fraction_modified_reads.plus.wig", sep = ""), header = F,fill = F)  %>%
    dplyr::rename(position = 1, mod = 2) %>%
    dplyr::filter(position %in% interesting_zone)   
  
  # >  return values
  return(fraction)
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................input files: 16S modification positions, fasta file
#.................................HVO
interesting_postions_hvo <- c(1598192+910, 1598192+1352, 1598192+1432, 1598192+1450,1598192+1451)
fasta_file_hvo <- here("data/genome_data/hvo.fasta")

#.................................HVO_dKSGA
interesting_postions_hvo <- c(1598192+910, 1598192+1352, 1598192+1432, 1598192+1450,1598192+1451)
fasta_file_hvo <- here("data/genome_data/hvo.fasta")

#.................................PFU
interesting_postions_pfu <- c(120675+938,120675+1379,120675+1469, 120675+1486,120675+1487)
fasta_file_pfu <- here("data/genome_data/pfu.fasta")

#.................................ECOLI
interesting_postions_ecoli <- c(223771+516, 223771+527, 223771 + 966, 223771 + 967, 223771 + 1207, 223771 + 1402, 223771 + 1407, 223771 + 1498, 223771 + 1516, 223771 + 1518, 223771 + 1519)
fasta_file_ecoli <- here("data/genome_data/ecoli.fasta")


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# STATISTICS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................calculate % of probability that a position is modified based on de_novo and model_compare tombo calculations
#.................................HVO
hvo_model   <- calculate_fraction_after_tombo(organism = "hvo", model_type = "model_testing", data_folder = paste(here(), "/data/tombo_data/hvo/wt/", sep = ""))
hvo_denovo  <- calculate_fraction_after_tombo(organism = "hvo", model_type = "de_novo", data_folder = paste(here(), "/data/tombo_data/hvo/wt/", sep = ""))

sum(hvo_model$mod > 0.5)/nrow(hvo_model)*100
which(hvo_model$mod > 0.50)
sum(hvo_denovo$mod > 0.5)/nrow(hvo_denovo)*100

#.................................HVO_dKSGA
hvo_dksga_model   <- calculate_fraction_after_tombo(organism = "hvo", model_type = "model_testing", data_folder = paste(here(), "/data/tombo_data/hvo/dksga/", sep = ""))
hvo_dksga_denovo  <- calculate_fraction_after_tombo(organism = "hvo", model_type = "de_novo", data_folder = paste(here(), "/data/tombo_data/hvo/dksga/", sep = ""))
sum(hvo_dksga_model$mod > 0.5)/nrow(hvo_dksga_model)
sum(hvo_dksga_denovo$mod > 0.5)/nrow(hvo_dksga_denovo)

#.................................PFU
pfu_model  <- calculate_fraction_after_tombo(organism = "pfu", model_type = "model_testing", data_folder = paste(here(), "/data/tombo_data/pfu/", sep = ""))
pfu_denovo <- calculate_fraction_after_tombo(organism = "pfu", model_type = "de_novo", data_folder = paste(here(), "/data/tombo_data/pfu/", sep = ""))
sum(pfu_model$mod > 0.5)/nrow(pfu_model)*100
sum(pfu_denovo$mod > 0.5)/nrow(pfu_denovo)*100

#.................................ECOLI
ecoli_model  <- calculate_fraction_after_tombo(organism = "ecoli", model_type = "model_testing", data_folder = paste(here(), "/data/tombo_data/ecoli/", sep = ""))
ecoli_denovo <- calculate_fraction_after_tombo(organism = "ecoli", model_type = "de_novo", data_folder = paste(here(), "/data/tombo_data/ecoli/", sep = ""))
sum(ecoli_model$mod > 0.5)/nrow(ecoli_model)*100
which(ecoli_model$mod > 0.25)

sum(ecoli_denovo$mod > 0.5)/nrow(ecoli_denovo)*100

#...................................heatmap calculations
#.................................HVO
hvo_known_positions <- calc_mod_heatmap(interesting_postions_hvo, hvo_model)

#.................................PFU
pfu_known_positions <- calc_mod_heatmap(interesting_postions_pfu, pfu_model)

#.................................ECOLI
ecoli_known_positions <- calc_mod_heatmap(interesting_postions_ecoli, ecoli_model)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOT HEATMAPS WITH GGPLOT2
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#.................................heatmap plotting
heat_color_npg <- c(pal_npg()(10)[4],
                    pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    pal_npg()(10)[5],
                    pal_npg()(10)[1],
                    pal_npg()(10)[8])

#...............................HVO (Fig. 6e)
pdf(here("figures/whole16S_tombo_hvo.pdf"),
    width = 16, height = 3, paper = "special",onefile=FALSE)
plot_16S_heatmap(hvo_model, hvo_denovo, hvo_known_positions)
dev.off()

#...............................PFU (Fig. 6g)
pdf(here("figures/whole16S_tombo_pfu.pdf"),
    width = 16, height = 3, paper = "special",onefile=FALSE)
plot_16S_heatmap(pfu_model, pfu_denovo, pfu_known_positions)
dev.off()

#...............................ECOLI (Fig. 6f)
pdf(here("figures/whole16S_tombo_ecoli.pdf"),
    width = 16, height = 3, paper = "special",onefile=FALSE)
plot_16S_heatmap(ecoli_model, ecoli_denovo, ecoli_known_positions)
dev.off()

