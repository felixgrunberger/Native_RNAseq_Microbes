###########################################################################
###########################################################################
###
### STATISTIC PLOTS FROM SUMMARY FILE (GUPPY/POREPLEX)
###
###########################################################################
###########################################################################


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
packages <- c("tidyverse", "here", "ggthemes", "data.table", "ggExtra", 
              "Rsamtools", "GenomicAlignments", "seqTools", "Rsubread", 
              "ape", "DT", "ggpubr", "ggridges", "ggsci", "R.utils")
invisible(lapply(packages, require, character.only = TRUE))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................publication_theme_white for plotting
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

#...................................load and mutate summary file
wrapper_summary_table <- function (input_summary_file, identifier, barcodes_used){
  
  # read in summary file from guppy
  summary_table <- fread(input_summary_file) 
  
  # how to deal with different output from poreplex
  summary_table$group <- NA
  if("barcode" %in% colnames(summary_table)){
    summary_table <- summary_table %>%
      dplyr::filter(barcode %in% barcodes_used) %>%
      mutate(sequence_length_template = sequence_length,
             mean_qscore_template = mean_qscore,
             group = paste(identifier, barcode, sep = "_"))
  }
  
  # add sequencing data set as 'identifier'
  summary_table_selected <- summary_table %>%
    mutate(group = ifelse(is.na(group),identifier, group)) 
  
  return(summary_table_selected)
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................summary files were first compressed using gzip (-9 option for optimal compression) and are stored under summary data

#...................................load summary data and modify them
summary_files <- paste(here("data/summary_data/"), list.files(here("data/summary_data/")), sep = "")

for (i in seq_along(summary_files)){
  sample_name <- str_split_fixed(str_split_fixed(summary_files[i], "_seq", 2)[1],"summary_data/",2)[2]
  table_name  <- paste(sample_name, "_table", sep = "")
  assign(table_name, wrapper_summary_table(summary_files[i], sample_name))
}

#...................................merge all tables into one
summary_table <- rbindlist(list(pfu_tex_table,
                                pfu_notex_table,
                                ecoli_tex_table, 
                                ecoli_notex_table, 
                                hvo_tex_table, 
                                hvo_notex_table)) %>% 
  mutate(group = as.factor(group),
         start_time = as.numeric(start_time), 
         sequence_length_template = as.numeric(sequence_length_template),
         mean_qscore_template = as.numeric(mean_qscore_template),
         type = ifelse(calibration_strand_score > 0, "enolase", "genome"))

#...................................calculate read statistics
#.................................number of reads/bases, median read lengths/qualities
summary_table %>%
  group_by(group) %>%
  summarise(number_of_reads = n(),
            number_of_bases = sum(sequence_length_template, na.rm = TRUE),
            median_rawread_length = median(sequence_length_template, na.rm = TRUE),
            median_rawread_quality = median(mean_qscore_template, na.rm = TRUE))

#...................................regroup
summary_table$group <-  factor(summary_table$group, 
                                     levels = rev(c("ecoli_tex","ecoli_notex",
                                                    "pfu_tex", "pfu_notex",
                                                    "hvo_tex", "hvo_notex")))
                                     
                                     
#...................................set colors for plotting
my_npg_colors <- sample(c("#E64B35",
                          "#DC0000",
                          "#91D1C2",
                          "#00A087",
                          "#8491B4",
                          "#4DBBD5",
                          "#3C5488"))
my_2npg_colors  <- c("#E64B35", "#3C5488")

#...................................prepare for yield plot
summary_table_yield <- summary_table %>%
  mutate(timeindex = round((start_time + 3600)/1800)) %>%
  group_by(group, timeindex) %>% 
  summarise(yield = sum(sequence_length_template)/1000000000) %>%
  mutate(cumulative_yield = cumsum(yield))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOT WITH GGPLOT2
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................total yield plot over time (Supplementary Fig. 2a)
yield_plot <- ggplot(summary_table_yield, aes(x=timeindex/2, y = cumulative_yield, color = group, fill = group)) +
  geom_line(alpha = 0.8, size = 3 ) +
  geom_area(alpha = 0.1, position = "identity", color = NA) +
  xlab("Run time (h)") +
  ylab("Cumulative yield (Gbases)") +
  theme_Publication_white() +
  scale_color_manual(values = my_npg_colors) +
  scale_fill_manual(values = my_npg_colors) +
  scale_linetype_manual(values = c(4,1)) +
  scale_y_continuous(expand = c(0,0, 0.1, 0), breaks = c(0.0,0.2,0.4,0.6,0.8)) +
  scale_x_continuous(limits = c(0,48), expand = c(0,0)) +
  guides(color = guide_legend(title = ""), fill = guide_legend(title = ""), linetype = guide_legend(title = ""))

pdf(here("figures/total_yield.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
yield_plot
dev.off()

#...................................raw read length genome vs control (Supplementary Fig. 2b)
gg_length <- ggplot(data = summary_table, aes(x = sequence_length_template, y = group, fill = type, color = type)) +
  geom_density_ridges2(aes(height =..ndensity..), scale = 0.9, alpha = 0.8, size = 1) +
  theme_Publication_white() +
  scale_x_continuous(trans = "log10", limits = c(50,5000), breaks = c(100,1000,1314,3000,3000),expand = c(0, 0)) +
  ylab("") +
  scale_y_discrete(expand = c(0,0, 0.2, 0)) +
  xlab("Log10 Read length (nt)") +
  scale_fill_manual(values = my_2npg_colors) +
  scale_color_manual(values = my_2npg_colors) +
  guides(group = F, fill= F, color = F) +
  theme(axis.text.y = element_text(face = "italic"))  

pdf(here("figures/lengths_control_vs_genome.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_length
dev.off()


#...................................quality of raw reads genome vs control (Supplementary Fig. 2c)
gg_quality <- ggplot(data = summary_table, aes(x = mean_qscore_template, y = group, fill = type, color = type)) +
  geom_density_ridges2(aes(height =..ndensity..), scale = 0.9, alpha = 0.8, size = 1) +
  theme_Publication_white() +
  scale_x_continuous(limits = c(0,16), expand = c(0, 0)) +
  ylab("") +
  scale_y_discrete(expand = c(0,0, 0.2, 0)) +
  xlab("Read quality") +
  scale_fill_manual(values = my_2npg_colors) +
  scale_color_manual(values = my_2npg_colors) +
  guides(group = F, fill= F, color = F) +
  theme(axis.text.y = element_text(face = "italic"))  

pdf(here("figures/quality_control_vs_genome.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_quality
dev.off()



