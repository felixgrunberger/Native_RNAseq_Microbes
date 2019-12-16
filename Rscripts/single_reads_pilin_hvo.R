###########################################################################
###########################################################################
###
### SINGLE-READ ANALYSIS - PILIN HVO
###
###########################################################################
###########################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES AND PLOTTING FUNCTION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(here)
source(here("Rscripts/load_libraries.R"))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DATA ANALYIS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................extract reads mapping to a location (cds1962 = pilin)
load(here("data/tidy_data/hvo_notex_id_table"))

cds1962_table <- full_id_table %>%
  dplyr::filter(gene == "cds1962")

#...................................list can be extracted and used to extract reads based on read ids - visualized in a genome browser
list_names <- list(cds1962_table$minion_read_name)


# > coverage plot at this positions

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# COVERAGE CALCULATED PLOT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................genome positions
interesting_positions <- 1911045:1912101

#...................................get coverage files
coverage_tex   <- fread(here("data/coverage_data/hvo_tex_minus_depth.txt.gz")) %>%
  dplyr::filter(V2 %in% interesting_positions) %>%
  dplyr::rename(position = V2, coverage = V3) %>%
  mutate(dataset = "TEX")

coverage_notex <- fread(here("data/coverage_data/hvo_notex_minus_depth.txt.gz")) %>%
  dplyr::filter(V2 %in% interesting_positions) %>%
  dplyr::rename(position = V2, coverage = V3) %>%
  mutate(dataset = "NOTEX")

#...................................scale them 
rescale_factor <- max(coverage_tex$coverage)/max(coverage_notex$coverage)
coverage_notex$coverage <- coverage_notex$coverage * rescale_factor

#...................................change colors
heat_color_npg <- rev(c(pal_npg()(10)[4],
                        pal_npg()(10)[1]))

#...................................position specific coverage area plot
pdf(here("/figures/single_reads_pilin_hvo_area.pdf"), 
    width = 26, height = 6, paper = "special",onefile=FALSE)
gg_area_plot <- ggplot(data = coverage_notex, aes(x = position, y = coverage)) +
  geom_line(color = pal_npg()(10)[4],  size = .5) +
  geom_area(width = 1, fill = pal_npg()(10)[4], color = NA, alpha = 0.7) +
  geom_area(data = coverage_tex,width = 1, fill = pal_npg()(10)[1], alpha = 0.7, color = pal_npg()(10)[1]) +
  geom_line(data = coverage_tex,color = pal_npg()(10)[1],  size = .5) +
  theme_Publication_white() +
  theme(panel.grid.major = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), trans = "reverse") +
  geom_rect(aes(xmin = 1911379, xmax = 1911816, ymin = -10, ymax = -5), color = "black", fill = "black", alpha = 0.5) +
  xlab("") +
  ylab("coverage") +
  geom_vline(xintercept = c(1911908,1911364))
gg_area_plot
dev.off()

#...................................get ranges
ggp <- ggplot_build(gg_area_plot)
my.ggp.yrange <- ggp$layout$panel_scales_y[[1]]$range$range  # data range!
my.ggp.xrange <- ggp$layout$panel_scales_x[[1]]$range$range  # data range!

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# EXTRACT READ AND CHECK WHAT IS TRIMMED OFF DURING BASECALLING
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#................................modify for plotting
plot_group <- cds1962_table %>%
  arrange(start) %>%
  #arrange(desc(width)) %>%
  ungroup() %>%
  mutate(lfd_n = 1:n()) 

#................................add annotation information
#................................get TSS/TTS
operon_annotation <- fread(here("operon_data/hvo_notex_reads_for_operons.tsv")) %>%
  dplyr::filter(gene %in% cds1962_table$gene)

#................................get gff info
gff_annotation <- read.gff(here("data/genome_data/hvo.gff")) %>%
  dplyr::filter(type == "CDS") %>%
  rowwise() %>%
  mutate(id = str_split_fixed(str_split_fixed(attributes, ";Parent",2)[1], "ID=",2)[2]) %>%
  dplyr::filter(id == "cds1962")

#................................plot
pdf(here("figures/single_reads_pilin_hvo.pdf"),
    width = 16, height = 10, paper = "special",onefile=FALSE)
ggplot(data = plot_group) +
  geom_segment(aes(x = start, xend = end, yend = lfd_n, y = lfd_n), size = 0.7) +
  theme_Publication_white() +
  geom_vline(xintercept = c(gff_annotation$end,gff_annotation$start, 
                            operon_annotation$median_utr5, operon_annotation$median_utr3),
                            linetype = "dashed", alpha = 0.5) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_reverse() +
  xlab("") +
  ylab("") +
  scale_x_continuous(limits = c(-my.ggp.xrange[1], -my.ggp.xrange[2]), trans = "reverse", expand = c(0,0))
dev.off()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# U ENRICHMENT PLOT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#................................get sequence
fasta <- readDNAStringSet(here("data/genome_Data/hvo.fasta"))
interesting_positions <- seq(from = -my.ggp.xrange[2], to = -my.ggp.xrange[1], by =  1)
int_seq <- reverseComplement(fasta$`NC_013967.1 Haloferax volcanii DS2, complete genome`[interesting_positions])
int_seq <- RNAString(gsub(int_seq, pattern = "T", replacement = "U"))

#...............................set colors for ATCG
cols=pal_npg()(10)[c(1,2,10,3)]

helper <- Biostrings::as.data.frame(int_seq) %>%
  rownames_to_column("position") %>%
  as_tibble() %>%
  mutate(position = as.numeric(position),
         "A" = ifelse(x == "A", 1, 0),
         "U" = ifelse(x == "U", 1, 0),
         "G" = ifelse(x == "G", 1, 0),
         "C" = ifelse(x == "C", 1, 0)) %>%
  dplyr::select(position, A,U,G,C) %>%
  gather(position) %>%
  dplyr::rename(base = position) %>%
  mutate(position = rep(length(int_seq):1, 4)) %>%
  mutate(position = position + min(interesting_positions))

pdf(here("figures/single_reads_pilin_hvo_gc.pdf"),
    width = 26, height = 6, paper = "special",onefile=FALSE)
ggplot() +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_col(data = helper, aes(x = position, y = value, color = base, fill = base)) +
  theme_Publication_white() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(-my.ggp.xrange[1], -my.ggp.xrange[2]), trans = "reverse", expand = c(0,0)) +
  xlab("") +
  ylab("")
dev.off()

