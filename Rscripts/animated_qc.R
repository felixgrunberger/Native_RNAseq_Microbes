###########################################################################
###########################################################################
###
### ANIMATED QUALITY CONTROLS
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

#...................................read in and modify summary file
wrapper_summary_table <- function (input_summary_file, identifier, barcodes_used){
  summary_table <- fread(input_summary_file) 
  summary_table$group <- NA
  if("barcode" %in% colnames(summary_table)){
    summary_table <- summary_table %>%
      dplyr::filter(barcode %in% barcodes_used) %>%
      mutate(sequence_length_template = sequence_length,
             mean_qscore_template = mean_qscore,
             group = paste(identifier, barcode, sep = "_"))
  }
  summary_table_selected <- summary_table %>%
    mutate(group = ifelse(is.na(group),identifier, group)) 
  #dplyr::select(channel, start_time, sequence_length_template, mean_qscore_template, group, minion_read_name)
  return(summary_table_selected)
}

#...................................set heatmap colors
heat_color_npg <- c(pal_npg()(10)[4],
                    pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    pal_npg()(10)[5],
                    pal_npg()(10)[1],
                    pal_npg()(10)[8])

my_npg_colors <- c("#E56650", 
                   "#9BADDE", 
                   "#F2A5A3",
                   "#356758",
                   "#3D819A",
                   "#5FC2AD")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CUMULATIVE YIELD OVER THE TIME
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................read in summary data from guppy
summary_files <- paste(here("data/summary_data/"), list.files(here("data/summary_data/")), sep = "")

for (i in seq_along(summary_files)){
  sample_name <- str_split_fixed(str_split_fixed(summary_files[i], "_seq", 2)[1],"summary_data/",2)[2]
  table_name  <- paste(sample_name, "_table", sep = "")
  assign(table_name, wrapper_summary_table(summary_files[i], sample_name))
}

#...................................merge data
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

#...................................culculate yield values 
summary_table_yield <- summary_table %>%
  mutate(timeindex = round((start_time + 3600)/1800)) %>%
  group_by(group, timeindex) %>% 
  summarise(yield = sum(sequence_length_template)/1000000000) %>%
  mutate(cumulative_yield = cumsum(yield),
         timeindex= as.integer(timeindex)) %>%
  dplyr::filter(!is.na(timeindex)) %>%
                mutate(time = timeindex/2)

#...................................animated yield plot over time using gganimate 
hours_yield <- ggplot(summary_table_yield, aes(x=time, y = as.numeric(cumulative_yield), group = group, color = group, fill = group)) +
  geom_line(alpha = 0.8, size = 3 ) +
  xlab("Run time (h)") +
  geom_text(aes(x = 49, label = group), hjust = 0) + 
  geom_segment(aes(xend = 49, yend = cumulative_yield), linetype = 2, alpha = 0.5, colour = 'white') +
  ylab("Cumulative yield (Gbases)") +
  theme_Publication_white() +
  scale_color_manual(values = my_npg_colors) +
  scale_fill_manual(values = my_npg_colors) +
  scale_y_continuous(expand = c(0,0, 0.1, 0), breaks = c(0.0,0.2,0.4,0.6,0.8)) +
  scale_x_continuous(limits = c(0,60), expand = c(0,0)) +
  transition_reveal(time) 

hours_yield_gif <- animate(hours_yield)
anim_save(animation = hours_yield_gif, here("figures/animated_hours_yield.gif"))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MinKNOW rebuild using 1 example dataset
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................set channels 1:512
list_a <- c(1,15,13,11,9,7,5,3)
platelayout_a <- list()
platelayout_b <- list()
for(i in 1:8){
  platelayout_a <- rbind(platelayout_a,data.frame (rown = rep (8:1, 1),  coln = rep ((29-(4*(i-1))):(32-(4*(i-1))), each = 8), channel = (32*list_a[i]):(1 + 32*(list_a[i]-1))))
  platelayout_b <- rbind(platelayout_b,data.frame (rown = rep (16:9, 1), coln = rep ((29-(4*(i-1))):(32-(4*(i-1))), each = 8), channel = (32*(list_a[i]+1)):(1 + 32*((list_a[i]+1)-1))))
}
platelayout <- rbind(platelayout_a, platelayout_b)

#...................................get ecoli tex data as an example
data <- ecoli_tex_table %>% 
  mutate(timeindex = round((start_time + 3600)/1800)) %>%
  group_by(channel, timeindex) %>% 
  summarise(n = n()) %>%
  arrange(channel)

#...................................add data to channel info
numbers <- 1:512
datahelper <- list()
for (i in seq_along(numbers)){
  datahelper$channel <- c(datahelper$channel,rep(i, 98))
}
datahelper$timeindex <- rep(c(1:98), 512)
datahelper$value <- 1
datahelper <- as_tibble(datahelper)

for (i in seq_along(data$channel)){
  datahelper[which(datahelper$channel == data$channel[i] & datahelper$timeindex == data$timeindex[i]),3] <- data$n[i]
}

datahelper[datahelper$value == 1,3] <- 0

plot_data <- datahelper %>%
  group_by(channel) %>%
  mutate(new_n = cumsum(value)) %>%
  left_join(platelayout, by = "channel") %>%
  ungroup()

#...................................plot and animate ussing gganimate
read_heatmap <- ggplot(data = plot_data, aes(y = factor(rown),x = factor(coln))) +
  geom_tile(color  = "black", aes(fill = new_n),size = 0.25)  +
  scale_fill_gradientn(colours = heat_color_npg) +
  ggtitle("") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  labs(x=NULL, y = NULL) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 0.5, ticks = F, label = T)) +
  coord_fixed(ratio=1) +
  transition_states(states = timeindex, wrap = TRUE) +
  theme_Publication_white() +
  theme(axis.text = element_blank(), axis.title = element_blank()) +
  theme(legend.title = element_blank())

read_heatmap_animated <- animate(read_heatmap)
anim_save(animation = anim, here("figures/animated_throughput.gif"))
