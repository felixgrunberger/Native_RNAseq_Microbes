###########################################################################
###########################################################################
###
### UTR 16S/23S PROCESSING SITES IN ESCHERICHIA COLI - CO-OCCURENCE ANALYSIS                                    
###
###########################################################################
###########################################################################


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(here)
source(here("Rscripts/load_libraries.R"))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................coloring for density
heat_color_npg <- c(pal_npg()(10)[4],
                    pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    pal_npg()(10)[5],
                    pal_npg()(10)[1],
                    pal_npg()(10)[8])

#...................................load bam file
get_bam_input <- function(org){
  input_bam_file <- paste(here("/data/mapped_data/"),org, ".bam", sep = "")
  allReads <- readGAlignments(input_bam_file, use.names = T, param = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), tag=c("NM"), what=c("mapq", "flag")))
  allReads_table <- GenomicAlignments::as.data.frame(allReads) %>%
    mutate(minion_read_name = names(allReads)) 
  
  # > get read sequences
  param <- ScanBamParam(
    flag=scanBamFlag(isUnmappedQuery=FALSE),
    what="seq")
  res0 <- scanBam(input_bam_file,param = param)[[1]] # always list-of-lists
  allReads_sequence <- res0[["seq"]]                 # query widths
  allReads_sequence_table <- as.list(as.character(allReads_sequence))
  allReads_table$sequence <- unlist(allReads_sequence_table)
  allReads_table$n_char   <- nchar(allReads_table$sequence[1:length(allReads_table$sequence)])
  left  <- paste(str_split_fixed(string = allReads_table$cigar, pattern = "M", n = 2)[,1],"M", sep = "")
  right <- paste(str_split_fixed(string = allReads_table$cigar, pattern = "M", n = 2)[,2],"1M", sep = "")
  
  #................................calculate cigar tables / SOFT AND HARD CLIPPING!!!
  allReads_table$soft_l <- as_tibble(cigarOpTable(left))$S
  allReads_table$hard_l <- as_tibble(cigarOpTable(left))$H
  allReads_table$soft_r <- as_tibble(cigarOpTable(right))$S
  allReads_table$hard_r <- as_tibble(cigarOpTable(right))$H
  return(allReads_table)
  
}

#...................................modify bam table
mod_bam <- function(input_table, set_name){
  return(input_table %>%
           mutate(UTR5_16 = ifelse(hard_l == 0 & hard_r == 0 & soft_l < 10, (start - 12) - gff$start[11],(start) - gff$start[11]),
                  UTR3_16 = (end - gff$end[11])) %>%
           mutate(UTR5_23 = ifelse(hard_l == 0 & hard_r == 0 & soft_l < 10, (start - 12) - gff$start[12],(start) - gff$start[12]),
                  UTR3_23 = (end - gff$end[12])) %>%
           mutate(sample = set_name)
  )
}
#...................................two color
two_color_npg <- rev(c("#3C5488", "#00A087"))

#...................................for density plots
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................genome data
gff <- read.gff(here("data/genome_data/ecoli.gff")) %>%
  dplyr::filter(type == "rRNA")

fasta <- readDNAStringSet(here("data/genome_data/ecoli.fasta"))

#...................................input data
tex_data           <- get_bam_input("ecoli_tex")
notex_data         <- get_bam_input("ecoli_notex")

#...................................filter reads within rrnC
tex_data_f <- tex_data %>%
  dplyr::filter(start > (gff$start[11] - 400),
                start < (gff$end[11] + 100),
                end > (gff$start[11] - 400),
                end < (gff$end[11] + 100)) 

notex_data_f <- notex_data %>%
  dplyr::filter(start > (gff$start[11] - 400),
                start < (gff$end[11] + 100),
                end > (gff$start[11] - 400),
                end < (gff$end[11] + 100)) 

#...................................select read categories (window = +/-1)
window_s <- 1
tex_data_f_bins <- tex_data %>%
  dplyr::filter(njunc == 0, soft_l < 1, hard_l == 0, hard_r == 0) %>%
  mutate(group = ifelse(start %in% (gff$start[11]-293+12-window_s):(gff$start[11]-293+12+window_s), "START", 
                        ifelse(start %in% (gff$start[11]-175+12-window_s):(gff$start[11]-175+12+window_s), "START II", 
                               ifelse(start %in% (gff$start[11]-115+12-window_s):(gff$start[11]-115+12+window_s), "RNASE_III", 
                                      ifelse(start %in% (gff$start[11]-66+12-window_s):(gff$start[11]-66+12+window_s),"RNASE_E",
                                             ifelse(start >= (gff$start[11]-0+12-window_s) & end < (gff$end[11]+50),"RNASE_G","rest")))))) %>%
  dplyr::filter(group != "rest") %>%
  group_by(group) %>%
  mutate(dens = approxfun(density(end))(end))

#...................................how many reads in category
tex_data_f_bins %>% 
  group_by(group) %>%
  summarise(n = n())

#...................................regroup
tex_data_f_bins$group <- factor(tex_data_f_bins$group, levels = rev(c("START", "START II", "RNASE_III", "RNASE_E", "RNASE_G", "rest")))

#...................................16S plotting 
pdf(here("figures/ecoli_16S_start_end_groups.pdf"), 
    width = 3, height = 3, paper = "special",onefile=FALSE)
ggplot(data = tex_data_f_bins, aes(x = group, y = end, color = dens, fill = dens)) +
  geom_hline(yintercept = gff$start[11]-293, linetype = "dashed") +
  geom_hline(yintercept = gff$start[11]-175, linetype = "dashed") +
  geom_hline(yintercept = gff$start[11]-115, linetype = "dashed") +
  geom_hline(yintercept = gff$start[11]-66, linetype = "dashed") +
  geom_hline(yintercept = gff$start[11]-0, linetype = "dashed") +
  geom_hline(yintercept = gff$end[11]+33) +
  geom_hline(yintercept = gff$end[11]) +
  geom_tile(size = 1) +
  scale_y_continuous(limits = c(gff$start[11]-300, (gff$end[11]+33))) +
  coord_flip() +
  theme_Publication_white() +
  scale_fill_gradientn(colours = heat_color_npg) +
  scale_color_gradientn(colours = heat_color_npg) +
  ylab("5´position") +
  xlab("3´position")
dev.off()


#...................................exact position for single-read plots
pdf(here("figures/ecoli_16S_start_end_groups_hist.pdf"), 
    width = 5, height = 5, paper = "special",onefile=FALSE)
ggplot(data = tex_data_f_bins, aes( x = end, y = group, color = group, fill = dens)) +
  geom_density_ridges2(aes(height = ..ndensity..),
                       alpha = 1, size = 0.1, scale = 0.95, 
                       stat = "binline", draw_baseline = FALSE, 
                       binwidth = 1, color = NA, fill = pal_npg()(10)[4]) +
  geom_vline(xintercept = gff$end[11]+33) +
  geom_vline(xintercept = gff$end[11]+0) +
  scale_x_continuous(limits = c(gff$end[11]-40, (gff$end[11]+40)), expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("3`utr [nt]") +
  ggtitle("") +
  scale_y_discrete(expand = c(0,0)) +
  #scale_fill_npg() +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()

#...................................23S: grouping and plotting
window_s <- 1
tex_23_ecoli <- tex_data %>%
  dplyr::filter(njunc == 0, soft_l < 1, hard_l == 0, hard_r == 0) %>%
  mutate(group = ifelse(start %in% (gff$start[12]-7+12-window_s):(gff$start[12]-7+12+window_s), "RNASE_III", 
                        ifelse(start >= (gff$start[12]-0+12-window_s) & end < (gff$end[12]+50), "mature", "rest"))) %>%
  dplyr::filter(group != "rest") %>%
  group_by(group) %>%
  mutate(dens = approxfun(density(end))(end))

#...................................how many in category?
tex_23_ecoli %>% 
  group_by(group) %>%
  summarise(n = n())

#...................................plot
pdf(here("figures/ecoli_23S_start_end_groups.pdf"), 
    width = 3, height = 3, paper = "special",onefile=FALSE)
ggplot(data = tex_23_ecoli, aes(x = group, y = end, color = dens, fill = dens)) +
  geom_tile(size = 1) +
  scale_y_continuous(limits = c(gff$start[12]-100, (gff$end[12]+33))) +
  coord_flip() +
  geom_hline(yintercept = gff$start[12], linetype = "dashed") +
  geom_hline(yintercept = gff$end[12]-7, linetype = "dashed") +
  geom_hline(yintercept = gff$end[12]-0, linetype = "dashed") +
  theme_Publication_white() +
  scale_fill_gradientn(colours = heat_color_npg) +
  scale_color_gradientn(colours = heat_color_npg) +
  ylab("5´position") +
  xlab("3´position")
dev.off()

#...................................exact position for single-read plots
pdf(here("figures/ecoli_23S_start_end_groups_hist.pdf"), 
    width = 5, height = 5, paper = "special",onefile=FALSE)
ggplot(data = tex_23_ecoli, aes( x = end, y = group, color = group, fill = group)) +
  geom_density_ridges2(aes(height = ..ndensity..),
                       alpha = 1, size = 0.1, scale = 0.95, 
                       stat = "binline", draw_baseline = FALSE, 
                       binwidth = 1, color = NA, fill = pal_npg()(10)[4]) +
  geom_vline(xintercept = gff$end[12]+9) +
  geom_vline(xintercept = gff$end[12]+0) +
  scale_x_continuous(limits = c(gff$end[12]-40, (gff$end[12]+40)), expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("3`utr [nt]") +
  ggtitle("") +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_npg() +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()

#...................................category number
tex_23_ecoli %>%
  dplyr::filter(end > gff$end[12]-40, end < gff$end[12]+40) %>%
  group_by(group) %>%
  summarise(n = n())

#...................................top 5 enriched positions 3´end --> 16S
tex_data_f_bins %>%
  mutate(end = end - gff$end[11]) %>%
  group_by(group, end) %>%
  mutate(n = n()) %>%
  group_by(group) %>%
  mutate(total= sum(n()),
         perc = n/total*100) %>%
  arrange(desc(n)) %>%
  select(n, total,perc, group, end) %>%
  group_by(group) %>%
  top_n(10, perc) %>%
  distinct(end, group, .keep_all = T) %>%
  arrange(group)

#...................................top 5 enriched positions 3´end --> 23S
tex_23_ecoli %>%
  mutate(end = end - gff$end[12]) %>%
  group_by(group, end) %>%
  mutate(n = n()) %>%
  group_by(group) %>%
  mutate(total= sum(n()),
         perc = n/total*100) %>%
  arrange(desc(n)) %>%
  select(n, total,perc, group, end) %>%
  group_by(group) %>%
  top_n(10, perc) %>%
  distinct(end, group, .keep_all = T) %>%
  arrange(group)


#...................................top 5 enriched positions 3´end --> 5S
tex_5_ecoli %>%
  mutate(end = end - gff$end[13]) %>%
  group_by(group, end) %>%
  mutate(n = n()) %>%
  group_by(group) %>%
  mutate(total= sum(n()),
         perc = n/total*100) %>%
  arrange(desc(n)) %>%
  select(n, total,perc, group, end) %>%
  group_by(group) %>%
  top_n(200, perc) %>%
  distinct(end, group, .keep_all = T) %>%
  arrange(group)

#...................................5S: category & plotting 
window_s <- 1
tex_5_ecoli <- tex_data %>%
  dplyr::filter(njunc == 0, soft_l == 0, hard_l == 0, hard_r == 0) %>%
  mutate(group =ifelse(start %in% (gff$start[13]-3+12-window_s):(gff$start[13]-3+12+window_s), "RNASE_E", 
                       ifelse(start >= (gff$start[13]-0+12-window_s) & end <= (gff$end[13]+50), "mature", "rest"))) %>%
  dplyr::filter(group != "rest") %>%
  group_by(group) %>%
  mutate(dens = approxfun(density(end))(end))

#...................................how many in group
tex_5_ecoli %>%
  group_by(group) %>%
  summarise(n = n())

#...................................plotting
pdf(here("figures/ecoli_5S_start_end_groups.pdf"), 
    width = 3, height = 3, paper = "special",onefile=FALSE)
ggplot(data = tex_5_ecoli, aes(x = group, y = end, fill = dens, color = dens)) +
  geom_tile(size = 1) +
  scale_y_continuous(limits = c(gff$start[13]-3, (gff$end[13]+40))) +
  coord_flip() +
  geom_hline(yintercept = gff$end[13]-0, linetype = "dashed") +
  geom_hline(yintercept = gff$start[13]-0, linetype = "dashed") +
  geom_hline(yintercept = gff$end[13]+3, linetype = "dashed") +
  geom_hline(yintercept = gff$end[13]-40, linetype = "dashed") +
  theme_Publication_white() +
  scale_fill_gradientn(colours = heat_color_npg) +
  scale_color_gradientn(colours = heat_color_npg) +
  ylab("5´position") +
  xlab("3´position")
dev.off()


#...................................exact position for single-read plots
pdf(here("figures/200502_ecoli_5S_start_end_groups_hist.pdf"), 
    width = 5, height = 5, paper = "special",onefile=FALSE)
ggplot(data = tex_5_ecoli, aes( x = end, y = group, color = group, fill = group)) +
  geom_density_ridges2(aes(height = ..ndensity..),
                       alpha = 1, size = 0.1, scale = 0.95, 
                       stat = "binline", draw_baseline = FALSE, 
                       binwidth = 1, color = NA, fill = pal_npg()(10)[4]) +
  geom_vline(xintercept = gff$end[13]+3) +
  geom_vline(xintercept = gff$end[13]+0) +
  scale_x_continuous(limits = c(gff$end[13]-40, (gff$end[13]+40)), expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("3`utr [nt]") +
  ggtitle("") +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_npg() +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()
