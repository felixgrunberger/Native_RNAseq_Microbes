###########################################################################
###########################################################################
###
### TOMBO OUTPUT FRACTION ANALYSIS
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
library("patchwork")
#...................................plot positions after tombo calculations
plot_after_tombo <- function(fasta, position, zone, model_type, data_folder){
  fasta <- readDNAStringSet(fasta)
  names(fasta)[1] <- "genome"
  interesting_zone <- (position - zone):(position + zone)
  selected_sequence <- fasta$genome[interesting_zone]
  type <- model_type
  
  sample_wig_mean <- fread(paste(data_folder, type,".signal.sample.plus.wig", sep = ""), header = F,fill = F) %>%
    dplyr::rename(position = 1, value = 2) %>%
    mutate(group = "sample")
  
  control_wig_mean <- fread(paste(data_folder, type, ".signal.control.plus.wig", sep = ""), header = F,fill = F)  %>%
    dplyr::rename(position = 1, value = 2) %>%
    mutate(group = "control")   
  
  # modified fraction
  fraction <- fread(paste(data_folder, type,".fraction_modified_reads.plus.wig", sep = ""), header = F,fill = F)  %>%
    dplyr::rename(position = 1, mod = 2) %>%
    dplyr::filter(position %in% interesting_zone)   
  
  # de novo modified fraction
  fraction_de_novo <- fread(paste(data_folder, "de_novo",".fraction_modified_reads.plus.wig", sep = ""), header = F,fill = F)  %>%
    dplyr::rename(position = 1, mod = 2) %>%
    dplyr::filter(position %in% interesting_zone)  
  
  # difference 
  difference <- fread(paste(data_folder, type,".difference.plus.wig", sep = ""), header = F,fill = F)  %>%
    dplyr::rename(position = 1, value = 2) %>%
    mutate(group = "difference") %>%
    dplyr::filter(position %in% interesting_zone) 
  
  # sequence
  sequence <- fread(paste(data_folder, type, ".signal.sample.plus.wig", sep = ""), header = F,fill = F) %>%
    dplyr::rename(position = 1, value = 2) %>%
    dplyr::filter(position %in% interesting_zone) %>%
    rowwise() %>%
    mutate(value = as.character(fasta$genome[position])) 
  
  # combine
  control_wig_zone <- rbind(sample_wig_mean, control_wig_mean) %>%
    dplyr::filter(position %in% interesting_zone) %>%
    left_join(fraction, by = "position")
  
  # plot sequence
  plot1 <- ggplot(data = control_wig_zone, aes(x = position, y = value, fill = group)) +
    geom_vline(xintercept = position, size = 8, alpha = 0.4, color = "yellow") +
    theme_Publication_white() +
    scale_alpha_continuous(range = c(1,1)) +
    geom_text(data = sequence, aes(x = position, y = 1,label = value, color = value,fill = NA, alpha = 1), size = 6, fontface = "bold") +
    scale_fill_npg() +
    scale_color_npg() +
    scale_x_continuous(minor_breaks = seq(min(interesting_zone), max(interesting_zone), by = 1),
                       breaks = seq(min(interesting_zone), max(interesting_zone), by = 5),
                       labels = seq(min(interesting_zone), max(interesting_zone), by = 5)) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    xlab("") +
    ylab("Signal") +
    guides(label = F, alpha = F, color = F, fill = guide_legend(title = "")) +
    theme(plot.margin = unit(c(0,0,-.8,-.8), "cm")) +
    rremove("y.text") +
    rremove("y.ticks") +
    ylab("") +
    rremove("x.text") +
    rremove("x.ticks")
  
  # plot fraction
  plot2 <- ggplot(data = fraction, aes(x = position, y = mod)) +
    geom_vline(xintercept = position, size = 6, alpha = 0.4, color = "yellow") +
    geom_line(size = 1.5) +
    geom_line(data = fraction_de_novo, aes(x = position, y = mod), size = 1.5, linetype = "dashed", alpha = 0.7) +
    scale_fill_npg() +
    scale_color_npg() +
    theme_Publication_white() +
    scale_x_continuous(minor_breaks = seq(min(interesting_zone), max(interesting_zone), by = 1),
                       breaks = seq(min(interesting_zone), max(interesting_zone), by = 5),
                       labels = seq(min(interesting_zone), max(interesting_zone), by = 5)) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    xlab("") +
    theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
    ylab("modified \nfraction") +
    rremove("x.text") +
    rremove("x.ticks")
  
  # plot points
  plot3 <- ggplot(data = control_wig_zone, aes(x = position, y = value, fill = group)) +
    geom_vline(xintercept = position, size = 6, alpha = 0.4, color = "yellow") +
    theme_Publication_white() +
    geom_point(shape = 21, size = 10, aes(alpha = mod, color = group), color = "white") +
    scale_alpha_continuous(range = c(0.5,1)) +
    #geom_text(data = sequence, aes(x = position, y = -2, label = value, color = value,fill = NA, alpha = 1), size = 6, fontface = "bold") +
    scale_fill_npg() +
    scale_color_npg() +
    scale_x_continuous(minor_breaks = seq(min(interesting_zone), max(interesting_zone), by = 1),
                       breaks = seq(min(interesting_zone), max(interesting_zone), by = 5),
                       labels = seq(min(interesting_zone), max(interesting_zone), by = 5)) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    xlab("") +
    ylab("Signal") +
    guides(label = F, alpha = F, color = F, fill = guide_legend(title = "")) +
    theme(plot.margin = unit(c(0,0,0,0), "cm")) 
  
  position_1 <- plot1 + plot2 + plot3 + plot_layout(ncol = 1, heights = c(1,1,3)) + xlab("Position")
  return(position_1)
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
# PLOT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#.................................how to 
plot_after_tombo(fasta_file_pfu, # fasta file of the organims interested in 
                 1598192+1450,   # position of modified nucleotide
                 20,             # position left and right adjacent to modified base
                 "model_testing",# model calulated (de novo or sample compare (=model_testing))
                 data_folder = paste(here(), "/data/tombo_data/hvo/wt/", sep = "")) # where are data stored
                 
