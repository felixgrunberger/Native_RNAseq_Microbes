###########################################################################
###########################################################################
###
### TRANSCRIPTION START SITES
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

#...................................violin plot
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#...................................is a gene start/end or rest of an operon
type_of_gene_in_operon <- function(genes_in_operon_table, type = c("first", "last")){
  plus  <- str_split(genes_in_operon_table$genes_in_operon[genes_in_operon_table$strand_operon == "+"], ",")
  minus <- str_split(genes_in_operon_table$genes_in_operon[genes_in_operon_table$strand_operon == "-"], ",")
  if(type == "first"){
    plus_gene  <- sapply(plus, head, 1)
    minus_gene <- sapply(minus, tail, 1)
    gene <- as.data.frame(c(plus_gene, minus_gene)) %>%
      mutate(start_operon = T) %>%
      dplyr::rename(gene = 1)
  }else{
    plus_gene  <- sapply(plus, tail, 1)
    minus_gene <- sapply(minus, head, 1)
    gene <- as.data.frame(c(plus_gene, minus_gene)) %>%
      mutate(end_operon = T) %>%
      dplyr::rename(gene = 1)
  }
  return(gene)
}

#...................................for 5´UTR violin plots
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#...................................detect TSS and calculate utr5 length
detect_tss <- function(filtered_operon_table, collapsed_operon_table){
  
  # check if a gene starts an operon or ends it
  first_gene <- type_of_gene_in_operon(collapsed_operon_table, "first")
  last_gene  <- type_of_gene_in_operon(collapsed_operon_table, "last")
  
  # calculate 5´utr length (gene_start - TSS)
  utr5_operon <- filtered_operon_table %>%
    dplyr::filter(mapped_type == "CDS") %>%
    left_join(first_gene, by = "gene") %>%
    left_join(last_gene, by = "gene") %>%
    group_by(gene) %>%
    mutate(counts = n(),
           start_operon = as.character(start_operon),
           utr5_length = ifelse(strand == "+", start_gene - median_utr5, median_utr5 - end_gene),
           group = ifelse(is.na(start_operon) == T, "rest", "start")) %>%
    dplyr::select(counts,median_utr5, gene, utr5_length,  strand, group, end_operon) %>%
    dplyr::filter(!is.na(median_utr5)) %>%
    distinct()
}

#...................................for density plots
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


#...................................plot DNA MEME motif
motif_plotter_dna <- function(motif_set){
  motif_set <-   c(as.character(droplevels(motif_set$V1)))
  ggplot() + 
    geom_logo(motif_set, font = "helvetica_bold", col_scheme = color_scale, seq_type = "dna") + 
    theme_logo() +
    theme_Publication_white() +
    theme(panel.grid.major = element_line(colour = NA),
          axis.ticks.x = element_line(colour = NA), 
          axis.text.x = element_text(size = 0)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(limits = c(0,2), expand = c(0,0))
}

#...................................plot correlation ILLLUMINA/ONT
gg_utr5_corr <- function(input_set, min_v, max_v){
  ggplot(data = input_set, aes(x = utr5_length, y = as.numeric(UTRlength), color = density, shape = counts > 5)) +
    geom_abline(intercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
    geom_abline(intercept = 12, color = "black", linetype = "dashed", alpha = 0.8) +
    geom_point(alpha = 0.7, size = 3.5, stroke = 1) +
    scale_x_continuous(limits = c(min_v, max_v)) +
    scale_y_continuous(limits = c(min_v, max_v)) +
    xlab("5´UTR - Native RNAseq (nt)") +
    ylab("5´UTR - Illumina dRNAseq (nt)") +
    scale_shape_manual(values=c(21, 15))+
    ggtitle("") +
    theme_Publication_white() +
    scale_color_gradientn(colours = heat_color_npg) +
    guides(fill = F, alpha = F) + 
    coord_equal() +
    stat_cor(method = "spearman", color = "black") +
    guides(fill = guide_legend(title = "", override.aes = list(alpha=1)), 
           color = guide_legend(title = ""))
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD TSS DATA FROM PUBLICATIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load published TSS data 

#.................................ESCHERICHIA COLI (combine with genome data to extract names)
ecoli_gff <- read.gff(file = here("data/genome_data/ecoli.gff")) %>%
  dplyr::filter(type == "CDS") %>%
  mutate(id_name = str_split_fixed(str_split_fixed(attributes, "ID=",2)[,2], ";Parent", 2)[,1],
         parent = str_split_fixed(str_split_fixed(attributes, "Parent=gene-",2)[,2], ";Dbxref", 2)[,1]) %>%
  dplyr::select(id_name, parent)

# tss from refernce data sets commented in the manuscript
ecoli_tss <- readxl::read_xlsx(path = here("data/tss_data/ecoli_tss.xlsx"), sheet = "TSS Map MasterTable", skip = 2) %>%
  as_tibble() %>%
  dplyr::select(Pos, Strand, Condition, detected, enriched, Locus_tag, UTRlength, GeneLength, Primary) %>%
  dplyr::filter(Primary == 1,
                Condition == "LB_0.4") %>%
  left_join(ecoli_gff, by = c("Locus_tag" = "parent"))

#.................................PYROCOCCUS FURIOSUS
# tss from refernce data sets commented in the manuscript
pfu_tss <- fread(here("data/tss_data/pfu_tss.tsv")) %>%
  dplyr::filter(Primary == 1) %>%
  dplyr::select(Locus_tag, UTRlength, Pos, Strand) 

#.................................HALOFERAX VOLCANII (combine with genome data to extract names)
hvo_gff <- here("data/genome_data/hvo.gff")

hvo_gff_table_gene <- read.gff(hvo_gff) %>%
  as_tibble() %>%
  dplyr::filter(type %in% "gene") %>%
  mutate(old_locus_tag = str_split_fixed(attributes, "old_locus_tag=", 2)[,2],
         GeneID = str_split_fixed(str_split_fixed(attributes, "GeneID:", 2)[,2], ";Name=",2)[,1]) %>%
  dplyr::rename(strand_gene = strand)

hvo_gff_table_cds <- read.gff(hvo_gff) %>%
  as_tibble() %>%
  dplyr::filter(type %in% "CDS") %>%
  mutate(ID = str_split_fixed(str_split_fixed(attributes, "ID=",2)[,2], ";Parent", 2)[,1],
         GeneID = str_split_fixed(str_split_fixed(attributes, "GeneID:", 2)[,2], ";Name=",2)[,1]) %>%
  dplyr::rename(strand_gene = strand)

hvo_gff_table_combined <- left_join(hvo_gff_table_cds, hvo_gff_table_gene, by = "GeneID") %>%
  mutate(strand_gene = strand_gene.x,
         start_gene = start.x,
         end_gene = end.x) %>%
  dplyr::select(start_gene, end_gene, strand_gene, ID, old_locus_tag)

# tss from refernce data sets commented in the manuscript
hvo_tss <- readxl::read_xlsx(path = here("data/tss_data/hvo_tss.xlsx"), sheet = "Suppl. Table S1") %>%
  as_tibble() %>%
  dplyr::filter(is.na(Novel)) %>%
  dplyr::select(`leaderless start distance`, `UTR Length`, Strand, `HVO (Correct or downstream feature)`) %>%
  mutate(UTRlength = ifelse(is.na(`leaderless start distance`), `UTR Length`, `leaderless start distance`)) %>%
  dplyr::rename(Gene = `HVO (Correct or downstream feature)`) %>%
  dplyr::filter(!is.na(Gene)) %>%
  dplyr::select(UTRlength, Gene) %>%
  left_join(hvo_gff_table_combined, by = c("Gene" = "old_locus_tag"))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD ONT PREDICTED DATA and combine with ILLUMINA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#...................................load transcriptional unit annotation
filtered_ids  <- paste(here("data/operon_data/"), list.files(here("data/operon_data/"),pattern = "for_operons.tsv"), sep = "")
collapsed_ids <- paste(here("data/operon_data/"), list.files(here("data/operon_data/"),pattern = "tex_operons.tsv"), sep = "")

#...................................combine ONT and ILLUMINA data
#.................................ESCHERICHIA COLI
ecoli_tss_table <- detect_tss(filtered_operon_table  = read_tsv(filtered_ids[1]), 
                              collapsed_operon_table = read_tsv(collapsed_ids[1])) %>%
  full_join(ecoli_tss,by = c("gene" = "id_name"))

ecoli_tss_table_ONT <- detect_tss(filtered_operon_table  = read_tsv(filtered_ids[1]), 
                                  collapsed_operon_table = read_tsv(collapsed_ids[1 ])) 


#.................................PYROCOCCUS FURIOSUS
pfu_tss_table <- detect_tss(filtered_operon_table  = read_tsv(filtered_ids[4]), 
                            collapsed_operon_table = read_tsv(collapsed_ids[4])) %>%
  mutate(new_gene = substr(gene, 1,15)) %>%
  full_join(pfu_tss,by = c("new_gene" = "Locus_tag"))

pfu_tss_table_ONT <- detect_tss(filtered_operon_table  = read_tsv(filtered_ids[4]), 
                            collapsed_operon_table = read_tsv(collapsed_ids[4])) %>%
  mutate(new_gene = substr(gene, 1,15)) 

#.................................HALOFERAX VOLCANII
hvo_tss_table <- detect_tss(filtered_operon_table  = read_tsv(filtered_ids[3]), 
                            collapsed_operon_table = read_tsv(collapsed_ids[3])) %>%
  full_join(hvo_tss,by = c("gene" = "ID"))

hvo_tss_table_ONT <- detect_tss(filtered_operon_table  = read_tsv(filtered_ids[3]), 
                                collapsed_operon_table = read_tsv(collapsed_ids[3])) 

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# UPSET-R GROUP COMPARISON
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................calculate intersects (off all possible variations) and save in new list
upsetr_table <- c(`PFU_ONT` = nrow(pfu_tss_table) - nrow(pfu_tss), 
                  `PFU_ILLUMINA` = nrow(pfu_tss_table) - nrow(pfu_tss_table_ONT), 
                  `PFU_ONT&PFU_ILLUMINA` = nrow(pfu_tss) + nrow(pfu_tss_table_ONT) - nrow(pfu_tss_table),
                  `ECOLI_ONT` = nrow(ecoli_tss_table) - nrow(ecoli_tss), 
                  `ECOLI_ILLUMINA` = nrow(ecoli_tss_table) - nrow(ecoli_tss_table_ONT), 
                  `ECOLI_ONT&ECOLI_ILLUMINA` = nrow(ecoli_tss) + nrow(ecoli_tss_table_ONT) - nrow(ecoli_tss_table),
                  `HVO_ONT` = nrow(hvo_tss_table) - nrow(hvo_tss), 
                  `HVO_ILLUMINA` = nrow(hvo_tss_table) - nrow(hvo_tss_table_ONT), 
                  `HVO_ONT&HVO_ILLUMINA` = nrow(hvo_tss) + nrow(hvo_tss_table_ONT) - nrow(hvo_tss_table))

#...................................group comparison (see Supplementary Fig. 6a)
pdf(here("figures/tss_intersections.pdf"),
    width = 14, height = 14, paper = "special",onefile=FALSE)
upset(data = fromExpression(upsetr_table), 
      sets = rev(c("ECOLI_ONT", "ECOLI_ILLUMINA", "PFU_ONT", "PFU_ILLUMINA", "HVO_ONT", "HVO_ILLUMINA")),
      keep.order = TRUE, order.by = c("degree"),
      sets.x.label = "Number of TSS in category", mainbar.y.label = "Intersections", 
      line.size = 3, point.size = 11, 
      text.scale = c(1.6, 1.2, 1.2, 1.2, 2, 1.6), 
      main.bar.color = pal_npg()(10)[c(4,6,5,4,4,6,6,5,5)], 
      sets.bar.color = "gray60",
      matrix.color = pal_npg()(10)[4])
dev.off()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 5´UTR comparison
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................modify input tables from ILLUMINA AND ONT (filter for utrlength)

#.................................ESCHERICHIA COLI
tss_ecoli <- inner_join(ecoli_tss, subset(ecoli_tss_table_ONT, counts >= 1), by = c("id_name" = "gene")) %>%
  distinct(id_name, .keep_all = TRUE) %>%
  mutate(nt_change = ifelse(Strand == "+", Pos - median_utr5, median_utr5 - Pos),
         sequencing_set = "ecoli") %>%
  mutate(UTRlength = as.numeric(UTRlength))

#.................................PYROCOCUCCUS FURIOSUS
tss_pfu <- inner_join(pfu_tss, subset(pfu_tss_table_ONT, counts >= 1), by = c("Locus_tag" = "new_gene")) %>%
  distinct(Locus_tag, .keep_all = TRUE) %>%
  mutate(nt_change = ifelse(Strand == "+", Pos - median_utr5, median_utr5 - Pos),
         sequencing_set = "pfu") 

#...................................HALOFERAX VOLCANII
tss_hvo <- inner_join(hvo_tss, subset(hvo_tss_table_ONT, counts >= 1), by = c("ID" = "gene")) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  mutate(Pos = ifelse(strand_gene == "+", start_gene - UTRlength, end_gene + UTRlength)) %>%
  mutate(nt_change = ifelse(strand_gene == "+", Pos - median_utr5, median_utr5 - Pos),
         sequencing_set = "hvo") 


#...................................merge all
tss_all <- bind_rows(tss_ecoli %>%
                       dplyr::select(nt_change, sequencing_set), 
                     tss_pfu %>%
                       dplyr::select(nt_change, sequencing_set), 
                     tss_hvo %>%
                       dplyr::select(nt_change, sequencing_set))

#...................................plot distance between ONT uncorrected 5´end and Illumina ends (Fig. 2a)
pdf(here("figures/utr5_uncorrected_distance.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggplot(tss_all, aes(x = nt_change)) +
  facet_grid(rows = vars(sequencing_set)) +
  geom_vline(xintercept = -12) +
  geom_histogram(binwidth = 1, color = "white", lwd = .3) +
  scale_x_continuous(limits = c(-40,40)) +
  theme_Publication_white()
dev.off()  


#.................................ESCHERICHIA COLI
ecoli_tss_ILLUMINA_filtered <- ecoli_tss %>%
  mutate(organism = "ecoli", set = "illumina") %>%
  dplyr::rename(utrlength = UTRlength) %>%
  dplyr::select(utrlength,organism, set)

ecoli_tss_ONT_filtered <- ecoli_tss_table_ONT %>%
  dplyr::filter(counts > 2) %>%
  ungroup() %>%
  mutate(organism = "ecoli", set = "nanopore", utr5_length = as.numeric(utr5_length)) %>%
  dplyr::rename(utrlength = utr5_length) %>%
  dplyr::select(utrlength,organism, set)

#.................................PYROCOCCUS FURIOSUS
pfu_tss_ILLUMINA_filtered <- pfu_tss %>%
  mutate(organism = "pfu", set = "illumina") %>%
  dplyr::rename(utrlength = UTRlength) %>%
  dplyr::select(utrlength,organism, set)

pfu_tss_ONT_filtered <- pfu_tss_table_ONT %>%
  dplyr::filter(counts > 2) %>%
  ungroup() %>%
  mutate(organism = "pfu", set = "nanopore", utr5_length = as.numeric(utr5_length)) %>%
  dplyr::rename(utrlength = utr5_length) %>%
  dplyr::select(utrlength,organism, set)

#.................................HALOFERAX VOLCANII
hvo_tss_ILLUMINA_filtered <- hvo_tss %>%
  mutate(organism = "hvo", set = "illumina") %>%
  dplyr::rename(utrlength = UTRlength) %>%
  dplyr::select(utrlength,organism, set)

hvo_tss_ONT_filtered <- hvo_tss_table_ONT %>%
  dplyr::filter(counts > 2) %>%
  ungroup() %>%
  mutate(organism = "hvo", set = "nanopore", utr5_length = as.numeric(utr5_length)) %>%
  dplyr::rename(utrlength = utr5_length) %>%
  dplyr::select(utrlength,organism, set)

#...................................merge all tables (uncorrected TSS for ONT data)
tss_all_uncorrected <- rbind(ecoli_tss_ILLUMINA_filtered, ecoli_tss_ONT_filtered, 
                             pfu_tss_ILLUMINA_filtered, pfu_tss_ONT_filtered, 
                             hvo_tss_ILLUMINA_filtered, hvo_tss_ONT_filtered) %>%
  mutate(utrlength = as.numeric(utrlength))

#...................................merge all tables (+12 corrected TSS for ONT data, filtering out many outliers based on too low coverage)
tss_all_corrected <- rbind(ecoli_tss_ILLUMINA_filtered, ecoli_tss_ONT_filtered, 
                           pfu_tss_ILLUMINA_filtered, pfu_tss_ONT_filtered, 
                           hvo_tss_ILLUMINA_filtered, hvo_tss_ONT_filtered) %>%
  mutate(utrlength = as.numeric(utrlength)) %>%
  dplyr::filter(utrlength <=300 & utrlength >= -20) %>%
  mutate(utrlength = ifelse(set == "nanopore", utrlength + 12, utrlength)) 

#...................................reorder groups
tss_all_uncorrected$organism <-  factor(tss_all_uncorrected$organism, 
                            levels = rev(c("ecoli", "pfu", "hvo")))

tss_all_corrected$organism <-  factor(tss_all_corrected$organism, 
                                        levels = rev(c("ecoli", "pfu", "hvo")))


#...................................set colors and plot using geom_densities
heat_color_npg <- rev(c(pal_npg()(10)[4],
                        pal_npg()(10)[7]))

#...................................plot corrected 5´UTR of three datasets (each TEX treated) in comparison to reference sets (Fig. 2b)
gg_utr_corrected <- ggplot(data = tss_all_corrected, aes(y = utrlength, x = organism, fill = set, color = set)) +
  geom_split_violin(scale = "width", trim = F, alpha = 0.5, size = 1) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = position_dodge(),
               geom = "crossbar", width = 0.2, color = "black") +
  scale_y_continuous(limits = c(-20,300), breaks = c(0,50,100,200, 300), expand = c(0,0)) +
  coord_flip() +
  theme_Publication_white() +
  xlab("5´ UTR length (nt)") +
  xlab("") +
  scale_fill_manual(values = heat_color_npg) +
  scale_color_manual(values = heat_color_npg)



pdf(here("figures/utr5_corrected.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_utr_corrected
dev.off()

#...................................median values
tss_all_corrected %>%
  group_by(organism, set) %>%
  summarise(median_length = median(utrlength))

#...................................filter out outliers for correct color coding
tss_ecoli_plot <- tss_ecoli %>%
  dplyr::filter(utr5_length > -21 & utr5_length < 301)
tss_pfu_plot <- tss_pfu %>%
  dplyr::filter(utr5_length > -21 & utr5_length < 301)
tss_hvo_plot <- tss_hvo %>%
  dplyr::filter(utr5_length > -21 & utr5_length < 301)

#...................................coloring
heat_color_npg <- c(pal_npg()(10)[4],
                    pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    pal_npg()(10)[5],
                    pal_npg()(10)[1],
                    pal_npg()(10)[8])

#...................................how many data for plotting? n?
nrow(tss_ecoli_plot)
nrow(tss_pfu_plot)
nrow(tss_hvo_plot)

#...................................calculate density
tss_ecoli_plot$density <- get_density(tss_ecoli_plot$utr5_length, as.numeric(tss_ecoli_plot$UTRlength))
tss_pfu_plot$density <- get_density(tss_pfu_plot$utr5_length, as.numeric(tss_pfu_plot$UTRlength))
tss_hvo_plot$density <- get_density(tss_hvo_plot$utr5_length, as.numeric(tss_hvo_plot$UTRlength))

#...................................plot ecoli full (Supplementary Fig. 6b left)
pdf(here("figures/utr5_ecoli_full.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_utr5_corr(tss_ecoli_plot, min_v = -20, max_v = 300)
dev.off()


#...................................plot ecoli small (Supplementary Fig. 6b right)
tss_ecoli_plot_filtered <- tss_ecoli_plot %>%
  dplyr::filter(utr5_length > -21 & utr5_length < 51)
tss_ecoli_plot_filtered$density <- get_density(tss_ecoli_plot_filtered$utr5_length, as.numeric(tss_ecoli_plot_filtered$UTRlength))

pdf(here("figures/utr5_ecoli_small.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_utr5_corr(tss_ecoli_plot_filtered, min_v = -20, max_v = 50)
dev.off()

#...................................plot pfu (Supplementary Fig. 6c left)
pdf(here("figures/utr5_pfu_full.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_utr5_corr(tss_pfu_plot, min_v = -20, max_v = 300)
dev.off()

#...................................plot pfu small (Supplementary Fig. 6c right)
tss_pfu_plot_filtered <- tss_pfu_plot %>%
  dplyr::filter(utr5_length > -21 & utr5_length < 51)
tss_pfu_plot_filtered$density <- get_density(tss_pfu_plot_filtered$utr5_length, as.numeric(tss_pfu_plot_filtered$UTRlength))

pdf(here("figures/utr5_pfu_small.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_utr5_corr(tss_pfu_plot_filtered, min_v = -20, max_v = 50)
dev.off()

#...................................plot hvo (Supplementary Fig. 6d left)
pdf(here("figures/utr5_hvo_full.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_utr5_corr(tss_hvo_plot, min_v = -20, max_v = 300)
dev.off()

#...................................plot hvo small (Supplementary Fig. 6d right)
tss_hvo_plot_filtered <- tss_hvo_plot %>%
  dplyr::filter(utr5_length > -21 & utr5_length < 51)
tss_hvo_plot_filtered$density <- get_density(tss_hvo_plot_filtered$utr5_length, as.numeric(tss_hvo_plot_filtered$UTRlength))

pdf(here("figures/utr5_hvo_small.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_utr5_corr(tss_hvo_plot_filtered, min_v = -20, max_v = 50)
dev.off()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SEQUENCING DEPTH INFLUENCE
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................calculate correlations
pear_cor_counts <- function(dataset, organism){
  counts <- list()
  cor <- list()
  cor_set <- NULL
  for(i in 1:100){
    subset <- dataset %>% mutate(utr5_length = utr5_length + 12) %>%
      dplyr::filter(UTRlength < 300, UTRlength >= 0,utr5_length <300, utr5_length >= 0, counts >= i)
    counts[i] <- as.integer(i)
    cor[i] <- cor(subset$UTRlength,
                  subset$utr5_length, method = "pearson")
  }
  cor_set <- data.table(counts = unlist(counts), cor = unlist(cor)) %>%
    as_tibble() %>% mutate(group = organism)
}

cor_pfu   <- pear_cor_counts(tss_pfu, "pfu")
cor_hvo   <- pear_cor_counts(tss_hvo, "hvo")
cor_ecoli <- pear_cor_counts(tss_ecoli, "ecoli")
cor_all   <- rbind(cor_pfu, cor_hvo, cor_ecoli)

#...................................plot correlations
ggplot(data = cor_pfu, aes(x = counts, y = cor, group = group, color = group, fill = group)) +
  geom_smooth(se = F, span = 0.1) +
  geom_ribbon(alpha = 0.2,aes(ymin = 0,ymax = predict(loess(cor ~ counts, span = 0.1)))) +
  scale_fill_npg() +
  scale_color_npg() +
  theme_Publication_white() +
  scale_x_continuous(limits = c(1,100),expand = expand_scale(add = c(0,0))) +
  geom_vline(xintercept = 5, linetype = "dashed")



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TEX vs NOTEX HVO INFLUENCE
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load data for TEX set
hvo_tex_tss_table <- detect_tss(filtered_operon_table  = read_tsv(filtered_ids[3]), 
                            collapsed_operon_table = read_tsv(collapsed_ids[3])) 

#...................................load data for NOTEX set
hvo_notex_tss_table <- detect_tss(filtered_operon_table  = read_tsv(filtered_ids[2]), 
                                  collapsed_operon_table = read_tsv(collapsed_ids[2])) 

#...................................combine data all
hvo_tex_vs_notex <- inner_join(hvo_tex_tss_table %>% mutate(tex_utr = utr5_length), hvo_notex_tss_table %>% mutate(notex_utr = utr5_length), by = c("gene")) %>%
  dplyr::filter(tex_utr > -21 & tex_utr < 301 & notex_utr > -21 & notex_utr < 301)
hvo_tex_vs_notex$density <- get_density(hvo_tex_vs_notex$tex_utr, hvo_tex_vs_notex$notex_utr)

#...................................plot TEX vs NOTEX full (Supplementary Fig. 6e left)
hvo_tex_vs_notex_plot_full <- ggplot(data = hvo_tex_vs_notex, aes(x = tex_utr, notex_utr, color = density, shape = counts.x > 5)) +
  geom_abline(intercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
  geom_point(alpha = 0.7, size = 3.5, stroke = 1) +
  scale_x_continuous(limits = c(-20, 300)) +
  scale_y_continuous(limits = c(-20, 300)) +
  xlab("5`UTR length MinION [nt]") +
  scale_shape_manual(values=c(21, 15))+
  ggtitle("") +
  theme_Publication_white() +
  scale_color_gradientn(colours = heat_color_npg) +
  guides(fill = F, alpha = F) + 
  coord_equal() +
  stat_cor(method = "spearman", color = "black") +
  guides(fill = guide_legend(title = "", override.aes = list(alpha=1)), 
         color = guide_legend(title = ""))

pdf(here("figures/utr5_hvo_tex_vs_notex_full.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
hvo_tex_vs_notex_plot_full
dev.off()


#...................................plot TEX vs NOTEX full (Supplementary Fig. 6e right)
hvo_tex_vs_notex <- inner_join(hvo_tex_tss %>% mutate(tex_utr = utr5_length), hvo_notex_tss %>% mutate(notex_utr = utr5_length), by = c("gene")) %>%
  dplyr::filter(tex_utr > -21 & tex_utr < 21 & notex_utr > -21 & notex_utr < 21)
hvo_tex_vs_notex$density <- get_density(hvo_tex_vs_notex$tex_utr, hvo_tex_vs_notex$notex_utr)

hvo_tex_vs_notex_plot_small <- ggplot(data = hvo_tex_vs_notex, aes(x = tex_utr, notex_utr, color = density, shape = counts.x > 5)) +
  geom_abline(intercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
  geom_point(alpha = 0.7, size = 3.5, stroke = 1) +
  scale_x_continuous(limits = c(-20, 20)) +
  scale_y_continuous(limits = c(-20, 20)) +
  xlab("5`UTR length MinION [nt]") +
  scale_shape_manual(values=c(21, 15))+
  ggtitle("") +
  theme_Publication_white() +
  scale_color_gradientn(colours = heat_color_npg) +
  guides(fill = F, alpha = F) + 
  coord_equal() +
  stat_cor(method = "spearman", color = "black") +
  guides(fill = guide_legend(title = "", override.aes = list(alpha=1)), 
         color = guide_legend(title = ""))

pdf(here("figures/utr5_hvo_tex_vs_notex_small.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
hvo_tex_vs_notex_plot_small
dev.off()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MOTIF DETECTION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................get sequences for MEME search (based on ONT table)
#.................................is a gene start or end of an transcriptional unit?
first_gene_ecoli <- type_of_gene_in_operon(fread(collapsed_ids[1]), "first")
last_gene_ecoli  <- type_of_gene_in_operon(fread(collapsed_ids[1]), "last")
first_gene_pfu   <- type_of_gene_in_operon(fread(collapsed_ids[4]), "first")
last_gene_pfu    <- type_of_gene_in_operon(fread(collapsed_ids[4]), "last")
first_gene_hvo   <- type_of_gene_in_operon(fread(collapsed_ids[3]), "first")
last_gene_hvo    <- type_of_gene_in_operon(fread(collapsed_ids[3]), "last")

#.................................load fasta files
fasta_file  <- paste(here("data/genome_data/"), list.files(here("data/genome_data/"),pattern = ".fasta"), sep = "")
ecoli_fasta <- readDNAStringSet(fasta_file[1])
pfu_fasta   <- readDNAStringSet(fasta_file[5])
hvo_fasta   <- readDNAStringSet(fasta_file[3])

#.................................filter single-read table for protein coding gene mapping reads | add position info
#..retrieve sequence in a window from -58 to -11

#...............................E. COLI
utr5_operon_ecoli <- fread(filtered_ids[1]) %>%
  dplyr::filter(mapped_type == "CDS") %>%
  left_join(first_gene_ecoli, by = "gene") %>%
  left_join(last_gene_ecoli, by = "gene") %>%
  group_by(gene) %>%
  mutate(start_operon = as.character(start_operon),
         utr5_length = ifelse(strand == "+", start_gene - median_utr5, median_utr5 - end_gene),
         group = ifelse(is.na(start_operon) == T, "rest", "start")) %>%
  dplyr::select(median_utr5, gene, utr5_length,  strand, group, end_operon) %>%
  dplyr::filter(!is.na(median_utr5)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  rowwise() %>%
  mutate(promoter_sequence = ifelse(strand == "+", as.character(ecoli_fasta$`U00096.2 Escherichia coli str. K-12 substr. MG1655, complete genome`[(median_utr5 - 46 - 12):(median_utr5 + 1 - 12)]),
                                                                    as.character(reverseComplement(ecoli_fasta$`U00096.2 Escherichia coli str. K-12 substr. MG1655, complete genome`[(median_utr5 - 1 + 12):(median_utr5 + 46 + 12)])))) %>%
  dplyr::filter(group == "start")

#...............................P. FURIOSUS
utr5_operon_pfu <- fread(filtered_ids[4]) %>%
  dplyr::filter(mapped_type == "CDS") %>%
  left_join(first_gene_pfu, by = "gene") %>%
  left_join(last_gene_pfu, by = "gene") %>%
  group_by(gene) %>%
  mutate(start_operon = as.character(start_operon),
         utr5_length = ifelse(strand == "+", start_gene - median_utr5, median_utr5 - end_gene),
         group = ifelse(is.na(start_operon) == T, "rest", "start")) %>%
  dplyr::select(median_utr5, gene, utr5_length,  strand, group, end_operon) %>%
  dplyr::filter(!is.na(median_utr5)) %>%
  distinct(gene, .keep_all = T) %>%
  rowwise() %>%
  mutate(promoter_sequence = ifelse(strand == "+", as.character(pfu_fasta$CP023154[(median_utr5 - 46 - 12):(median_utr5 + 1 - 12)]),
                                    as.character(reverseComplement(pfu_fasta$CP023154[(median_utr5 - 1 + 12):(median_utr5 + 46 + 12)])))) %>%
  mutate(new_gene = substr(gene, 1,15)) %>%
  dplyr::filter(group == "start") 

#...............................H. VOLCANII
hvo_gff_table_cds_chr <- hvo_gff_table_cds %>%
  dplyr::select(seqid, ID)
chr_lengths <- data.table()
chr_lengths$seqid <- levels(as.factor(hvo_gff_table_cds_chr$seqid))
chr_lengths$chr_size <- c(length(hvo_fasta$`NC_013964.1 Haloferax volcanii DS2 plasmid pHV3, complete sequence`),
                          length(hvo_fasta$`NC_013965.1 Haloferax volcanii DS2 plasmid pHV2, complete sequence`),
                          length(hvo_fasta$`NC_013966.1 Haloferax volcanii DS2 plasmid pHV4, complete sequence`),
                          length(hvo_fasta$`NC_013967.1 Haloferax volcanii DS2, complete genome`),
                          length(hvo_fasta$`NC_013968.1 Haloferax volcanii DS2 plasmid pHV1, complete sequence`))
chr_lengths <- as_tibble(chr_lengths)
hvo_gff_table_cds_chr <- hvo_gff_table_cds_chr %>%
  left_join(chr_lengths, by = "seqid")


utr5_operon_hvo <-  fread(filtered_ids[3]) %>%
  left_join(hvo_gff_table_cds_chr, by = c("gene" = "ID")) %>%
  dplyr::filter(mapped_type == "CDS") %>%
  left_join(first_gene_hvo, by = "gene") %>%
  left_join(last_gene_hvo, by = "gene") %>%
  group_by(gene) %>%
  mutate(start_operon = as.character(start_operon),
         utr5_length = ifelse(strand == "+", start_gene - median_utr5, median_utr5 - end_gene),
         group = ifelse(is.na(start_operon) == T, "rest", "start")) %>%
  dplyr::select(median_utr5, gene, utr5_length,  strand, group, end_operon, seqid, chr_size) %>%
  dplyr::filter(!is.na(median_utr5)) %>%
  distinct(gene, .keep_all = T) %>%
  group_by(seqid) %>%
  dplyr::filter(median_utr5 > 58, median_utr5 < chr_size) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(promoter_sequence = ifelse(strand == "+" & seqid == "NC_013964.1", as.character(hvo_fasta$`NC_013964.1 Haloferax volcanii DS2 plasmid pHV3, complete sequence`[(median_utr5 - 46 - 12):(median_utr5 + 1 - 12)]),
                                     ifelse(strand == "+" & seqid == "NC_013965.1", as.character(hvo_fasta$`NC_013965.1 Haloferax volcanii DS2 plasmid pHV2, complete sequence`[(median_utr5 - 46 - 12):(median_utr5 + 1 - 12)]),
                                           ifelse(strand == "+" & seqid == "NC_013966.1", as.character(hvo_fasta$`NC_013966.1 Haloferax volcanii DS2 plasmid pHV4, complete sequence`[(median_utr5 - 46 - 12):(median_utr5 + 1 - 12)]),
                                                  ifelse(strand == "+" & seqid == "NC_013967.1", as.character(hvo_fasta$`NC_013967.1 Haloferax volcanii DS2, complete genome`[(median_utr5 - 46 - 12):(median_utr5 + 1 - 12)]),
                                                         ifelse(strand == "+" & seqid == "NC_013968.1", as.character(hvo_fasta$`NC_013968.1 Haloferax volcanii DS2 plasmid pHV1, complete sequence`[(median_utr5 - 46 - 12):(median_utr5 + 1 - 12)]),
                                                                ifelse(strand == "-" & seqid == "NC_013968.1", as.character(reverseComplement(hvo_fasta$`NC_013968.1 Haloferax volcanii DS2 plasmid pHV1, complete sequence`[(median_utr5 - 1 + 12):(median_utr5 + 46 + 12)])),
                                                                       ifelse(strand == "-" & seqid == "NC_013967.1", as.character(reverseComplement(hvo_fasta$`NC_013967.1 Haloferax volcanii DS2, complete genome`[(median_utr5 - 1 + 12):(median_utr5 + 46 + 12)])),
                                                                              ifelse(strand == "-" & seqid == "NC_013966.1", as.character(reverseComplement(hvo_fasta$`NC_013966.1 Haloferax volcanii DS2 plasmid pHV4, complete sequence`[(median_utr5 - 1 + 12):(median_utr5 + 46 + 12)])),
                                                                                      ifelse(strand == "-" & seqid == "NC_013965.1", as.character(reverseComplement(hvo_fasta$`NC_013965.1 Haloferax volcanii DS2 plasmid pHV2, complete sequence`[(median_utr5 - 1 + 12):(median_utr5 + 46 + 12)])),
                                                                                            ifelse(strand == "-" & seqid == "NC_013964.1", as.character(reverseComplement(hvo_fasta$`NC_013964.1 Haloferax volcanii DS2 plasmid pHV3, complete sequence`[(median_utr5 - 1 + 12):(median_utr5 + 46 + 12)]))))))))))))) %>%
                                    
  mutate(new_gene = substr(gene, 1,15)) %>%
  dplyr::filter(group == "start") 


#.................................write to fasta
write.fasta(as.list(utr5_operon_ecoli$promoter_sequence), 
            as.list(utr5_operon_ecoli$gene), 
            file = here("data/meme_data/promoter_sequences_ecoli_operon_starting.fasta")) 

write.fasta(as.list(utr5_operon_pfu$promoter_sequence), 
            as.list(utr5_operon_pfu$gene), 
            file = here("data/meme_data/promoter_sequences_pfu_operon_starting.fasta")) 

write.fasta(as.list(utr5_operon_hvo$promoter_sequence), 
            as.list(utr5_operon_hvo$gene), 
            file = here("data/meme_data/promoter_sequences_hvo_operon_starting.fasta")) 


#.................................MEME command (outside R in console)
# meme >>INPUTFILE<< -dna -oc >>OUTPUTDIRECTORY<< -nostatus -time 18000 -maxsize 600000 -mod zoops -nmotifs 3 -minw 3 -maxw 16 -p 2 -brief 2000 -bfile >>CUSTOM_BGFILE_INTGERGENICREGIONS<<

#.................................save found motifs to fasta and plot motif in R

#...............................set colors for ATCG
color_scale = make_col_scheme(chars=c('A', 'T', 'C', 'G'), 
                              cols=pal_npg()(10)[c(1,3,2,10)])

#...............................load output from MEME and plot
motif_ecoli <- read.table(here("data/meme_data/ecoli_motif.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(!grepl("offset", V1)) 
motif_pfu <- read.table(here("data/meme_data/pfu_motif.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(!grepl("offset", V1)) 
motif_hvo <- read.table(here("data/meme_data/hvo_motif.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(!grepl("offset", V1)) 


pdf(here("figures/promoter_ecoli.pdf"),
    width = 6, height = 4, paper = "special",onefile=FALSE)
motif_plotter_dna(motif_ecoli)
dev.off()

pdf(here("figures/promoter_pfu.pdf"),
    width = 6, height = 4, paper = "special",onefile=FALSE)
motif_plotter_dna(motif_pfu)
dev.off()

pdf(here("figures/promoter_hvo.pdf"),
    width = 6, height = 4, paper = "special",onefile=FALSE)
motif_plotter_dna(motif_hvo)
dev.off()



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# WRITE TSS OUTPUT TO TABLE (CORRECTED POSITIONS -12)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...............................ecoli
#...................load annotation data
ecoli_annotation <- fread(here("data/genome_data/ecoli_annotation.tsv"))

#...................write table
ecoli_export <- ecoli_tss_table_ONT %>%
  mutate(median_utr5 = ifelse(strand == "+", median_utr5 - 12, median_utr5 +12),
         utr5_length = utr5_length +12) %>%
  left_join(ecoli_annotation, by = c("gene" = "id")) %>%
  dplyr::rename(GeneID = gene, TSS = median_utr5, number_of_reads = counts, annotation = name) %>%
  dplyr::select(GeneID, TSS, number_of_reads, strand, utr5_length, locus_tag, annotation) 

writexl::write_xlsx(x = ecoli_export, path = here("tables/tss_tables/tss_ecoli.xlsx"))

#...............................pfu
#...................load annotation data
pfu_annotation <- fread(here("data/genome_data/pfu_annotation.tsv"))

#...................write table
pfu_export <- pfu_tss_table_ONT %>%
  mutate(median_utr5 = ifelse(strand == "+", median_utr5 - 12, median_utr5 +12),
         utr5_length = utr5_length +12) %>%
  left_join(pfu_annotation, by = c("new_gene" = "gene")) %>%
  dplyr::rename(GeneID = gene, TSS = median_utr5, number_of_reads = counts, annotation = name) %>%
  dplyr::select(GeneID, TSS, number_of_reads, strand, utr5_length, old_name, annotation) 

writexl::write_xlsx(x = pfu_export, path = here("tables/tss_tables/tss_pfu.xlsx"))

#...............................hvo
#...................load annotation data
hvo_annotation <- fread(here("data/genome_data/hvo_annotation.tsv")) %>%
  rowwise() %>%
  mutate(name = str_split_fixed(name, ";",2)[1])

#...................write table
hvo_export <- hvo_tss_table_ONT %>%
  mutate(median_utr5 = ifelse(strand == "+", median_utr5 - 12, median_utr5 +12),
         utr5_length = utr5_length +12) %>%
  left_join(hvo_annotation, by = c("gene" = "id_name")) %>%
  dplyr::rename(GeneID = gene, TSS = median_utr5, number_of_reads = counts, annotation = name) %>%
  dplyr::select(GeneID, TSS, number_of_reads, strand, utr5_length, old_name, annotation) 

writexl::write_xlsx(x = hvo_export, path = here("tables/tss_tables/tss_hvo.xlsx"))


