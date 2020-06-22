###########################################################################
###########################################################################
###
### TRANSCRIPTION TERMINATION SITES - SINGLE GENES
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

#...................................load bam file
get_bam_input <- function(org, treatment){
  input_bam_file <- paste(here::here("data/mapped_data/"),org,"_", treatment, ".bam", sep = "")
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

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DATA ANALYIS - PILIN HVO
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#.....................................HVO PILIN (see Supplementary Fig. 8a,b,c)
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

#................................input bam file
allReads_table_tex   <- get_bam_input("hvo", "tex")
allReads_table_notex <- get_bam_input("hvo", "notex")

#................................filter reads that fall into genomic region
allReads_table_tex_filtered <- allReads_table_tex %>%
  dplyr::filter(start %in% interesting_positions, 
                end %in% interesting_positions)
allReads_table_notex_filtered <- allReads_table_notex %>%
  dplyr::filter(start %in% interesting_positions, 
                end %in% interesting_positions)

#...................................scale the data sets
rescale_factor <- max(coverage_tex$coverage)/max(coverage_notex$coverage)
coverage_notex$coverage <- coverage_notex$coverage * rescale_factor

#...................................change colors
heat_color_npg <- rev(c(pal_npg()(10)[4],
                        pal_npg()(10)[1]))

#...................................position specific coverage area plot (Supplementary Fig. 8a)
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


#...................................single-gene track plotting
#.................................get pilin reads
load(here("data/tidy_data/hvo_tex_id_table"))

cds1962_table <- full_id_table %>%
  dplyr::filter(gene == "cds1962") %>%
  arrange(start) %>%
  ungroup() %>%
  mutate(lfd_n = 1:n()) 

cds1962_table <- full_id_table %>%
  dplyr::filter(gene == "cds1962") %>%
  left_join(allReads_table_tex, by = "minion_read_name") %>%
  arrange(start.x) %>%
  ungroup() %>%
  mutate(lfd_n = 1:n()) 

#.................................add annotation information
#................................get TSS/TTS
operon_annotation <- fread(here("data/operon_data/hvo_notex_reads_for_operons.tsv")) %>%
  dplyr::filter(gene %in% cds1962_table$gene)

#................................get gff info
gff_annotation <- read.gff(here("data/genome_data/hvo.gff")) %>%
  dplyr::filter(type == "CDS") %>%
  rowwise() %>%
  mutate(id = str_split_fixed(str_split_fixed(attributes, ";Parent",2)[1], "ID=",2)[2]) %>%
  dplyr::filter(id == "cds1962")

#................................get ranges of previous plot
ggp <- ggplot_build(gg_area_plot)
my.ggp.yrange <- ggp$layout$panel_scales_y[[1]]$range$range  # data range y!
my.ggp.xrange <- ggp$layout$panel_scales_x[[1]]$range$range  # data range x!


#................................plot
pdf(here("figures/single_reads_pilin_hvo.pdf"),
    width = 16, height = 10, paper = "special",onefile=FALSE)
ggplot(data = cds1962_table) +
  geom_segment(aes(x = start.x, xend = end.x, yend = lfd_n, y = lfd_n), size = 0.7) +
  theme_Publication_white() +
  geom_vline(xintercept = c(gff_annotation$end,gff_annotation$start, 
                            operon_annotation$median_utr5, operon_annotation$median_utr3),
             linetype = "dashed", alpha = 0.5) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_reverse() +
  xlab("") +
  ylab("") 
  scale_x_continuous(limits = c(-my.ggp.xrange[1], -my.ggp.xrange[2]), trans = "reverse", expand = c(0,0))
dev.off()



# > same method used to look at alba/pilin transcripts in Pfu
