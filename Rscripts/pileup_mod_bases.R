###########################################################################
###########################################################################
###
### ANALYSIS OF MODIFIED BASES BASED ON WRONGLY ASSIGNED POSTITIONS
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

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#.................................heatmap plotting
four_color_npg <- c(pal_npg()(10)[3],
                    pal_npg()(10)[4],
                    pal_npg()(10)[2],
                    pal_npg()(10)[1]) 

#...................................for adding the sequence
base_annotation <- function(coord_left, coord_right, fasta_file){
  fasta <- readDNAStringSet(fasta_file)
  names(fasta) <- "genome"
  sequence <- as.character(fasta$genome[coord_left:coord_right])
  return(sequence)
}

calculate_window <- function(which_mod = c("n4", "other", "else", "homo_g", "all_ccg", "high_n4"), which_dataset){
  if(which_mod == "n4"){
    wanted <- n4_pfu  
  }
  else if(which_mod == "other"){
    wanted <- other_pfu 
  }
  else if(which_mod == "else"){
    wanted <- gff$start[1]:gff$end[1]
  }
  else if(which_mod == "homo_g"){
    wanted <- gggg_table$pos - 1
  }
  else if(which_mod == "all_ccg"){
    wanted <- ccg_table_16_table$pos
  }
  else if(which_mod == "high_n4"){
    wanted <- high_n4
  }
  
  
  # > list
  list_int <- list()
  for (i in 1:length(wanted)){
    list_int <- append(list_int, (wanted[i]-5):(wanted[i]+5))
  }
  sel_tab <- data.table(pos = unlist(list_int), 
                        type = "sel", 
                        pos_n = rep(-5:5,length(wanted)),
                        identifier = rep(wanted, each = 11)) %>%
    as_tibble()
  
  # > get dataset
  return(which_dataset %>%
           dplyr::select(pos, nuc_freq, category) %>%
           group_by(pos, category) %>%
           mutate(nuc_freq_total = sum(nuc_freq, na.rm = T)) %>%
           
           dplyr::left_join(sel_tab, by = "pos") %>%
           dplyr::filter(type == "sel") %>%
           ungroup() %>%
           dplyr::select(-pos, -type) %>%
           group_by(pos_n, category) %>%
           dplyr::mutate(mean_pos = mean(nuc_freq_total, na.rm = T),
                         sd_down =  mean_pos - sd(nuc_freq_total, na.rm = T),
                         sd_up =  mean_pos + sd(nuc_freq_total, na.rm = T),
                         median_pos = quantile(nuc_freq_total)[3],
                         q1 = quantile(nuc_freq_total)[2],
                         q3 = quantile(nuc_freq_total)[4]) %>%
           dplyr::mutate(group = which_mod))
  
}


#...................................gff
gff <- read.gff(here::here("data/genome_data/hvo.gff")) %>%
  dplyr::filter(type == "rRNA")

#...................................modify input data
calc_pilup_freq <- function(input, grouptype){
  total <- input %>% 
    dplyr::filter(pos == 1599640) %>%
    dplyr::select(reads_all)
  output <- input %>% 
    #dplyr::filter(pos > gff$start[1], pos < gff$end[1]) %>%
    dplyr::select(pos,A_fwd, T_fwd, G_fwd, C_fwd, insertions,deletions,ref,reads_all) %>%
    dplyr::rename(correct_base = ref) %>%
    as_tibble() %>%
    gather(nucleotide, count, A_fwd:deletions, factor_key=TRUE) %>%
    mutate(nucleotide = substr(nucleotide,1,1)) %>%
    arrange(pos) %>%
    mutate(strand = "+") %>%
    mutate(category = ifelse(nucleotide == correct_base, "correct", 
                             ifelse(nucleotide == "i", "insertion",
                                    ifelse(nucleotide == "d", "deletion","wrong")))) %>%
    group_by(pos) %>%
    mutate(nuc_freq = count/sum(count) * 100) %>%
    dplyr::select(pos,strand, reads_all,nucleotide ,count, correct_base, nuc_freq, category) %>%
    mutate(type = grouptype) %>%
    mutate(total_n = max(total$reads_all))
  
  return(output)
}

pysamstats_per_base_freq <- function(input, dataset_name){
  input %>% 
    dplyr::filter(pos %in% gff$start[1]:gff$end[1]) %>%
    dplyr::select(pos,A_fwd, T_fwd, G_fwd, C_fwd, insertions,deletions,ref) %>%
    dplyr::rename(correct_base = ref) %>%
    as_tibble() %>%
    gather(nucleotide, count, A_fwd:deletions, factor_key=TRUE) %>%
    mutate(nucleotide = substr(nucleotide,1,1)) %>%
    arrange(pos) %>%
    mutate(seqnames = "CP023154", which_label = "CP023154:120675-122193") %>%
    mutate(strand = "+") %>%
    mutate(category = ifelse(nucleotide == correct_base, "correct", 
                             ifelse(nucleotide == "i", "insertion",
                                    ifelse(nucleotide == "d", "deletion","wrong")))) %>%
    group_by(pos) %>%
    mutate(nuc_freq = count/sum(count) * 100) %>%
    mutate(set = dataset_name) %>%
    dplyr::select(set,seqnames, pos,strand, nucleotide ,count, which_label, correct_base, nuc_freq, category) 
  # dplyr::mutate(mod_type = ifelse(pos %in% n4_pfu, "n4", 
  #                                ifelse(pos %in% other_pfu, "other","else")))
}



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load pysamstats data
#.................................NOTEX
wt_allReads_pre           <- fread(here("data/pileups/primary_pre_rRNA_wt_pileup.tsv"))
wt_allReads_mature        <- fread(here("data/pileups/mature_rRNA_wt_pileup.tsv"))
wt_allReads_closed_circ   <- fread(here("data/pileups/closed_circ_rRNA_wt_pileup.tsv"))
wt_allReads_open_circ     <- fread(here("data/pileups/open_circ_rRNA_wt_pileup.tsv"))


#.................................calculate pileup frequencies
#...............................NOTEX
wt_pre_rRNA_pileup             <- calc_pilup_freq(wt_allReads_pre, "pre_rRNA")
wt_closed_circ_rRNA_pileup     <- calc_pilup_freq(wt_allReads_closed_circ, "closed_circ_rRNA")
wt_open_circ_rRNA_pileup       <- calc_pilup_freq(wt_allReads_open_circ, "open_circ_rRNA")
wt_mature_rRNA_pileup          <- calc_pilup_freq(wt_allReads_mature, "mature_rRNA")

#...................................per base frequencies
wt_allReads_pre_freq          <- pysamstats_per_base_freq(wt_allReads_pre, "pre")
wt_allReads_closed_circ_freq  <- pysamstats_per_base_freq(wt_allReads_closed_circ, "closed_circ")
wt_allReads_open_circ_freq    <- pysamstats_per_base_freq(wt_allReads_open_circ, "open_circ")
wt_allReads_mature_freq       <- pysamstats_per_base_freq(wt_allReads_mature, "mature")

#...................................frequencies in window -5/+5 - PRE
wt_allReads_pre_freq_window_n4 <- calculate_window("n4", wt_allReads_pre_freq)
wt_allReads_pre_freq_window_other <- calculate_window("other", wt_allReads_pre_freq)
wt_allReads_pre_freq_window_else <- calculate_window("else", wt_allReads_pre_freq)
wt_allReads_pre_freq_window_all <- rbind(wt_allReads_pre_freq_window_n4,
                                         wt_allReads_pre_freq_window_other,
                                         wt_allReads_pre_freq_window_else)


#...................................frequencies in window -5/+5 - CLOSED CIRC
wt_allReads_closed_circ_freq_window_n4 <- calculate_window("n4", wt_allReads_closed_circ_freq)
wt_allReads_closed_circ_freq_window_other <- calculate_window("other", wt_allReads_closed_circ_freq)
wt_allReads_closed_circ_freq_window_else <- calculate_window("else", wt_allReads_closed_circ_freq)
wt_allReads_closed_circ_freq_window_all <- rbind(wt_allReads_closed_circ_freq_window_n4,
                                                 wt_allReads_closed_circ_freq_window_other,
                                                 wt_allReads_closed_circ_freq_window_else)

#...................................frequencies in window -5/+5 - OPEN CIRC
wt_allReads_open_circ_freq_window_n4 <- calculate_window("n4", wt_allReads_open_circ_freq)
wt_allReads_open_circ_freq_window_other <- calculate_window("other", wt_allReads_open_circ_freq)
wt_allReads_open_circ_freq_window_else <- calculate_window("else", wt_allReads_open_circ_freq)
wt_allReads_open_circ_freq_window_all <- rbind(wt_allReads_open_circ_freq_window_n4,
                                               wt_allReads_open_circ_freq_window_other,
                                               wt_allReads_open_circ_freq_window_else)


#...................................frequencies in window -5/+5 - MATURE
wt_allReads_mature_freq_window_n4 <- calculate_window("n4", wt_allReads_mature_freq)
wt_allReads_mature_freq_window_other <- calculate_window("other", wt_allReads_mature_freq)
wt_allReads_mature_freq_window_else <- calculate_window("else", wt_allReads_mature_freq)
wt_allReads_mature_freq_window_all <- rbind(wt_allReads_mature_freq_window_n4,
                                            wt_allReads_mature_freq_window_other,
                                            wt_allReads_mature_freq_window_else)

wt_all_group_hvo <- rbind(wt_allReads_pre_freq_window_all %>% 
                            mutate(name = "pre"),
                          wt_allReads_closed_circ_freq_window_all %>% 
                            mutate(name = "closed_circ"),
                          wt_allReads_open_circ_freq_window_all %>% 
                            mutate(name = "open_circ"),
                          wt_allReads_mature_freq_window_all %>% 
                            mutate(name = "mature"))

#.................................how many reads in group
#...............................NOTEX
wt_pre_rRNA_pileup$total_n[1]
wt_closed_circ_rRNA_pileup$total_n[1]
wt_open_circ_rRNA_pileup$total_n[1]
wt_mature_rRNA_pileup$total_n[1]


#.................................combine groups
wt_pileup <- rbind(wt_pre_rRNA_pileup, 
                   wt_closed_circ_rRNA_pileup, 
                   wt_open_circ_rRNA_pileup,
                   wt_mature_rRNA_pileup) %>%
  mutate(mutant = as.factor("wt"))


#.................................modified 16S HVO rRNA positions
interesting_postions_hvo <- c(1598192+910, 1598192+1352, 1598192+1432, 1598192+1450,1598192+1451, 1598192+1442)

#...................................HVO m6A
selector <- 6
(interesting_postions_hvo[selector]-10)-gff$start[1]
(interesting_postions_hvo[selector]+20)-gff$start[1]

m6a_seq <- base_annotation(coord_left = (interesting_postions_hvo[selector]-10), 
                           coord_right = (interesting_postions_hvo[selector]+20),
                           fasta_file = here("data/genome_data/hvo.fasta"))

wt_pileup$category <- factor(wt_pileup$category, levels = c("correct",  "insertion","deletion","wrong"))
wt_pileup$type <- factor(wt_pileup$type, levels = c("pre_rRNA", "pre_rRNA2", "closed_circ_rRNA", "open_circ_rRNA", "mature_rRNA", "mature_rRNA2"))


#...................................plot
pdf(here("figures/pileup_wt_16S_hvo.pdf"),
    width = 3.5, height = 3.5, paper = "special",onefile=FALSE)
ggplot(data = wt_pileup, aes(x = pos, y = nuc_freq, color = category, fill = category)) +
  geom_col(width = 0.9) +
  facet_grid(rows = vars(type), cols = vars(mutant)) +
  scale_x_continuous(expand = c(0,0),
                     breaks = (interesting_postions_hvo[selector]-10):(interesting_postions_hvo[selector]+20), 
                     limits = c((interesting_postions_hvo[selector]-10),(interesting_postions_hvo[selector]+20)),
                     labels = strsplit(m6a_seq, "*")[[1]]) +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = c(interesting_postions_hvo[selector], interesting_postions_hvo[4])) +
  xlab("") +
  ylab("Frequency (%)") +
  scale_fill_manual(values = four_color_npg) +
  scale_color_manual(values = four_color_npg) +
  theme_Publication_white()
dev.off()

###group plots 
#.................................m6A position
wt_allReads_all_freq_window_all_m6A <- wt_pileup %>%
  dplyr::filter(identifier == interesting_postions_hvo[4]) %>%
  distinct(category, nuc_freq_total, pos_n, identifier, name, mean_pos, sd_up, sd_down) %>%
  dplyr::filter(category == "wrong") %>%
  ungroup()

pdf(here::here("figures/m6A_hvo_basecall.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggplot(data = wt_allReads_all_freq_window_all_m6A, 
       aes(x = pos_n, y = nuc_freq_total, fill = name, color = name)) +
  facet_grid(~category) +
  geom_line(size = 1.5) +
  theme_Publication_white() +
  scale_fill_manual(values = pal_npg()(10)[c(3,8,4, 10)]) +
  scale_color_manual(values = pal_npg()(10)[c(3,8,4, 10)]) +
  scale_x_continuous(breaks = -5:5,expand = c(0,0)) +
  scale_y_continuous(limits = c(0,50), expand = c(0,0)) +
  xlab("") +
  ylab("")
dev.off()

#.................................n4AC position
wt_allReads_all_freq_window_all_n4 <- wt_all_group_hvo %>%
  dplyr::filter(identifier == (interesting_postions_hvo[4]-8)) %>%
  distinct(category, nuc_freq_total, pos_n, identifier, name, mean_pos, sd_up, sd_down) %>%
  dplyr::filter(category == "wrong") %>%
  ungroup()

pdf(here::here("figures/N4AC_hvo_basecall.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggplot(data = wt_allReads_all_freq_window_all_n4, 
       aes(x = pos_n, y = nuc_freq_total, fill = name, color = name)) +
  facet_grid(~category) +
  geom_line(size = 1.5) +
  theme_Publication_white() +
  scale_fill_manual(values = pal_npg()(10)[c(3,8,4, 10)]) +
  scale_color_manual(values = pal_npg()(10)[c(3,8,4, 10)]) +
  scale_x_continuous(breaks = -5:5,expand = c(0,0)) +
  scale_y_continuous(limits = c(0,50), expand = c(0,0)) +
  xlab("") +
  ylab("")
dev.off()


#.................................region
wt_allReads_all_freq_window_all_n4 <- wt_all_group_hvo %>%
  dplyr::filter(pos_n == 0) %>%
  dplyr::filter(identifier %in% ((interesting_postions_hvo[4]-17):(interesting_postions_hvo[4]+11))) %>%
  distinct(category, nuc_freq_total, pos_n, identifier, name, mean_pos, sd_up, sd_down) %>%
  dplyr::filter(category == "wrong") %>%
  ungroup()

pdf(here::here("figures/N4AC_hvo_basecall_line.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
length(min(wt_allReads_all_freq_window_all_n4$identifier):max(wt_allReads_all_freq_window_all_n4$identifier))
ggplot(data = wt_allReads_all_freq_window_all_n4, 
       aes(x = identifier, y = nuc_freq_total, fill = name, color = name)) +
  facet_grid(~category) +
  geom_line(size = 1.5) +
  theme_Publication_white() +
  scale_fill_manual(values = pal_npg()(10)[c(3,8,4, 10)]) +
  scale_color_manual(values = pal_npg()(10)[c(3,8,4, 10)]) +
  scale_x_continuous(breaks = -5:5,expand = c(0,0)) +
  scale_y_continuous(limits = c(0,50), expand = c(0,0)) +
  xlab("") +
  ylab("")
dev.off()


