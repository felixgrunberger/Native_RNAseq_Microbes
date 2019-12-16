###########################################################################
###########################################################################
###
### SINGLE-READ ANALYSIS - ALBA PFU
###
###########################################################################
###########################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES AND PLOTTING FUNCTION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(here)
source(here("Rscripts/load_libraries.R"))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#................................load alba
gene_of_interest <- "PFDSM3638_09460"

#................................load tidy data 
load(here("data/tidy_data/pfu_tex_id_table"))

#................................filter for those reads
full_id_table_spliced <- full_id_table %>%
  dplyr::filter(gene %in% c(paste(gene_of_interest, ".p01", sep = ""))) %>%
  distinct(minion_read_name, .keep_all = TRUE)

#................................read in bam file again to get cigar string - information is stored in "N" tag (skipped region from reference)
input_bam_file <- here("data/mapped_data/pfu_tex.bam")
allReads <- readGAlignments(input_bam_file, use.names = T, param = ScanBamParam(tag=c("NM"), what="mapq"))
allReads_table <- GenomicAlignments::as.data.frame(allReads) %>%
  mutate(minion_read_name = names(allReads)) 

#................................filter reads
allReads_table_filtered <- allReads_table %>%
  dplyr::filter(minion_read_name %in% full_id_table_spliced$minion_read_name)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# preapare single-read plot
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#................................get genome data
gff <- read.gff(here("data/genome_data/pfu.gff"))
gff_region <-  gff %>%
  mutate(id_name = str_split_fixed(str_split_fixed(attributes, "ID=",2)[,2], ";Parent", 2)[,1]) %>%
  dplyr::select(id_name, start, end, strand) 

#................................get TSS/TTS
operon_annotation <- fread(here("data/operon_data/pfu_tex_reads_for_operons.tsv")) %>%
  dplyr::filter(gene %in% c(paste(gene_of_interest, ".p01", sep = "")))


#................................get lfdn
utr_long <- allReads_table_filtered %>%
  dplyr::filter(end %in% c((operon_annotation$median_utr5[1]-50):(operon_annotation$median_utr5[1]+50)),start > 1699400) %>%
  #sample_n(200) %>%
  arrange(desc(width)) %>%
  mutate(lfd_n = 1:n()) %>%
  head(200)

#................................plot
pdf(here("/figures/single_reads_alba_pfu.pdf"), 
    width = 16, height = 10, paper = "special",onefile=FALSE)
gg_area_plot <- ggplot(data = utr_long) +
  geom_segment(aes(x = start, xend = end, yend = lfd_n, y = lfd_n), size = 0.2) +
  theme_Publication_white() +
  geom_vline(xintercept = c(gff_region$start[gff_region$id_name == "PFDSM3638_09460.p01"],gff_region$end[gff_region$id_name == "PFDSM3638_09460.p01"]),color = "grey50", linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(gff_region$start[gff_region$id_name == "PFDSM3638_09455.p01"],gff_region$end[gff_region$id_name == "PFDSM3638_09455.p01"]),color = "blue", linetype = "dashed", alpha = 0.5, size = 1) +
  geom_vline(xintercept = c(operon_annotation$median_utr3[1], operon_annotation$median_utr5[1]),linetype = "dashed", color = "red1", alpha = 0.5) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_blank()) +
  scale_x_reverse() +
  scale_y_reverse()
gg_area_plot
dev.off()

#...................................get ranges
ggp <- ggplot_build(gg_area_plot)
my.ggp.yrange <- ggp$layout$panel_scales_y[[1]]$range$range  # data range!
my.ggp.xrange <- ggp$layout$panel_scales_x[[1]]$range$range  # data range!


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# COVERAGE CALCULATED PLOT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
interesting_positions <- seq(from = -my.ggp.xrange[2], to = -my.ggp.xrange[1], by =  1)

#...................................get coverage files
coverage_tex   <- fread(here("data/coverage_data/pfu_tex_minus_depth.txt.gz")) %>%
  dplyr::filter(V2 %in% interesting_positions) %>%
  dplyr::rename(position = V2, coverage = V3) %>%
  mutate(dataset = "TEX")

coverage_notex <- fread(here("data/coverage_data/pfu_notex_minus_depth.txt.gz")) %>%
  dplyr::filter(V2 %in% interesting_positions) %>%
  dplyr::rename(position = V2, coverage = V3) %>%
  mutate(dataset = "NOTEX")

#...................................rescale
rescale_factor <- max(coverage_tex$coverage)/max(coverage_notex$coverage)
coverage_notex$coverage <- coverage_notex$coverage * rescale_factor

#...................................change colors
heat_color_npg <- rev(c(pal_npg()(10)[4],
                        pal_npg()(10)[1]))

#...................................position specific coverage area plot
pdf(here("/figures/single_reads_alba_pfu_area.pdf"), 
    width = 26, height = 6, paper = "special",onefile=FALSE)
ggplot(data = coverage_tex, aes(x = position, y = coverage)) +
  geom_line(color = pal_npg()(10)[4],  size = .5) +
  geom_area(width = 1, fill = pal_npg()(10)[4], color = NA, alpha = 0.7) +
  geom_area(data = coverage_notex,width = 1, fill = pal_npg()(10)[1], alpha = 0.7, color = pal_npg()(10)[1]) +
  geom_line(data = coverage_notex,color = pal_npg()(10)[1],  size = .5) +
  theme_Publication_white() +
  theme(panel.grid.major = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), trans = "reverse") +
  geom_rect(aes(xmin = gff_region$start[gff_region$id_name == "PFDSM3638_09460.p01"], xmax = gff_region$end[gff_region$id_name == "PFDSM3638_09460.p01"], ymin = -10, ymax = -5), color = "black", fill = "black", alpha = 0.5) +
  geom_rect(aes(xmin = gff_region$start[gff_region$id_name == "PFDSM3638_09455.p01"], xmax = gff_region$end[gff_region$id_name == "PFDSM3638_09455.p01"], ymin = -10, ymax = -5), color = "grey50", fill = "black", alpha = 0.5) +
  xlab("") +
  ylab("coverage") +
  guides(color = F) +
  geom_vline(xintercept = c(operon_annotation$median_utr5[1] + 12, operon_annotation$median_utr3[1]))
dev.off()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# T ENRICHMENT PLOT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# get sequence
pfu_fasta <- readDNAStringSet(here("data/genome_data/pfu.fasta"))
int_seq <- reverseComplement(pfu_fasta$CP023154[interesting_positions])
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

pdf(here("/figures/single_reads_alba_pfu_gc.pdf"), 
    width = 26, height = 6, paper = "special",onefile=FALSE)
ggplot() +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_col(data = helper, aes(x = position, y = value, color = base, fill = base)) +
  theme_Publication_white() +
  scale_x_continuous(expand = c(0,0), trans = "reverse") +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
dev.off()


