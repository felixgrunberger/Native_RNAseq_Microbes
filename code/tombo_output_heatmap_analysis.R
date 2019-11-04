###########################################################################
###########################################################################
###
### TOMBO OUTPUT ANALYSIS
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

#...................................calculation % of modified reads based on new model calculated
calculate_fraction_after_tombo <- function(organism = c("pfu", "hvo","ecoli"),model_type = c("de_novo", "model_testing"), data_folder){
  # > set positions
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
interesting_postions_pfu <- c(120675+938,120675+1379,120675+1469, 120675+1487,120675+1488)
fasta_file_hvo <- here("data/genome_data/pfu.fasta")

#.................................ECOLI
interesting_postions_ecoli <- c(223771+516, 223771+527, 223771 + 966, 223771 + 967, 223771 + 1207, 223771 + 1402, 223771 + 1407, 223771 + 1498, 223771 + 1516, 223771 + 1518, 223771 + 1519)
fasta_file_hvo <- here("data/genome_data/ecoli.fasta")


#...................................calculate % of probability that a position is modified based on de_novo and model_compare tombo calculations
#.................................HVO
hvo_model   <- calculate_fraction_after_tombo(organism = "hvo", model_type = "model_testing", data_folder = paste(here(), "/data/tombo_data/hvo/wt/", sep = ""))
hvo_denovo  <- calculate_fraction_after_tombo(organism = "hvo", model_type = "de_novo", data_folder = paste(here(), "/data/tombo_data/hvo/wt/", sep = ""))
hvo_wt_vs_ksga  <- calculate_fraction_after_tombo(organism = "hvo", model_type = "model_testing_wt_vs_ksga", data_folder = paste(here(), "/data/tombo_data/hvo/wt_vs_dksga/", sep = ""))

sum(hvo_model$mod > 0.5)/nrow(hvo_model)*100
which(hvo_model$mod > 0.50)
sum(hvo_denovo$mod > 0.5)/nrow(hvo_denovo)*100
sum(hvo_wt_vs_ksga$mod > 0.5)/nrow(hvo_wt_vs_ksga)*100
which(hvo_wt_vs_ksga$mod > 0.5)

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
ggplot(data = hvo_model, aes(x = position, y = "tombo model", fill = mod, color = mod)) +
  geom_tile(size = 1) +
  geom_tile(data = hvo_denovo, aes(x = position, y = "tombo de novo",  fill = mod, color = mod), size = 1) +
  geom_tile(data = hvo_known_positions, aes(x = position, y = "literature",  fill = mod, color = mod), size = 1) +
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
dev.off()

#...............................PFU (Fig. 6g)
pdf(here("figures/whole16S_tombo_pfu.pdf"),
    width = 16, height = 3, paper = "special",onefile=FALSE)
ggplot(data = pfu_model, aes(x = position, y = "tombo model", fill = mod, color = mod)) +
  geom_tile(size = 1) +
  geom_tile(data = pfu_denovo, aes(x = position, y = "tombo de novo",  fill = mod, color = mod), size = 1) +
  geom_tile(data = pfu_known_positions, aes(x = position, y = "literature",  fill = mod, color = mod), size = 1) +
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
dev.off()

#...............................ECOLI (Fig. 6f)
pdf(here("figures/whole16S_tombo_ecoli.pdf"),
    width = 16, height = 3, paper = "special",onefile=FALSE)
ggplot(data = ecoli_model, aes(x = position, y = "tombo model", fill = mod,color = mod)) +
  geom_tile(size = 1) +
  geom_tile(data = ecoli_denovo, aes(x = position, y = "tombo de novo",  color = mod, fill = mod), size = 1) +
  geom_tile(data = ecoli_known_positions, aes(x = position, y = "literature",  color = mod,fill = mod), size = 1) +
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
dev.off()
