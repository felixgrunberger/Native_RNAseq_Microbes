###########################################################################
###########################################################################
###
### TRANSCRIPT ABUNDANCE COMPARISON ONT/ILLUMINA
###
###########################################################################
###########################################################################


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
packages <- c("DESeq2", "ggpubr", "rtracklayer", "ggeconodist", 
              "tidyverse", "here", "ggthemes", "data.table", "ggExtra", 
              "Rsamtools", "GenomicAlignments", "seqTools", "Rsubread", 
              "ape", "DT", "ggpubr", "ggridges", "ggsci", "readxl")
invisible(lapply(packages, require, character.only = TRUE))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................for density plots
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#...............................DESeq2 transcript abundance transformation
deseq2_estimation_and_rlog <- function(input_counts_table){
  colnames_data <- c("sequencing")
  colData <- matrix(ncol = 1, nrow = 2)
  colnames(colData) <- colnames_data 
  colData <- data.table(colData) %>%
    mutate(sequencing = as.factor(c("nanopore", "illumina"))) %>%
    as_tibble()
  
  # > make DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = input_counts_table[,1:2], 
                                        colData = colData,
                                        design = ~sequencing)
  
  # > prefilter data set: remove rows (genes) with no counts
  dds     <- dds[ rowSums(counts(dds)) > 1, ]
  # > rld transformation
  rld      <- DESeq2::rlog( dds ,blind = T)
  # > plot object
  plot_rld <- as_tibble(assay(rld)[, 1:2]) 
}

#...............................plit transcript abundace datas for Illumina and ONT
transcript_abundace_ggcorr <- function(input_count_table){
  ggplot(data = input_count_table, aes(x = ont_counts, y = illumina_counts, color = density)) +
    geom_abline(alpha = 0.5, linetype = "dashed", color = "black", size = 1) +
    scale_fill_gradientn(colours = heat_color_npg) +
    scale_color_gradientn(colours = heat_color_npg) +
    geom_point(alpha = 0.3, size = 2) +
    theme_Publication_white() +
    guides(alpha = F) +
    guides(fill = guide_colorbar(barwidth = 5, barheight = 0.5, ticks = T, label = T),
           color = guide_colorbar(barwidth = 5, barheight = 0.5, ticks = T, label = T)) +
    ylab("Expression level of Illumina RNAseq \n(rlog counts)") +
    xlab("Expression level of native RNAseq \n(rlog counts)") +
    ggtitle("") +
    scale_x_continuous(limits = c(0,17)) + 
    scale_y_continuous(limits = c(0,17)) + 
    coord_equal() +
    stat_cor(method = "spearman", label.x = 0, label.y = 16,color = "black") +
    guides(color = guide_colorbar(title = "counts",barwidth = 15, barheight = 0.5, ticks = T, label = T)) 
}


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
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD AND TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................ONT transcript abundance data for TEX samples 
counts_files_tex <- paste(here("data/tidy_data/"), list.files(here("data/tidy_data/"),pattern = "_tex_gene_table"), sep = "")

#...................................get sample names
sample_names <- unlist(lapply(counts_files_tex, FUN=function(x){str_split_fixed(str_split_fixed(x, "_gene", 2)[1],"tidy_data/",2)[2]}))

#...................................write to tsv files
for (i in seq_along(sample_names)){
  load(counts_files_tex[i])
  full_gene_table <- full_gene_table %>%
    dplyr::select(GeneID, counts) %>%
    arrange(desc(counts))
  write_tsv(full_gene_table, path = paste(here("data/counts_data/"),sample_names[i],"_ONT.tsv", sep = ""))
}

#...................................get Illumina comparison data as described in the M&M section
# > Pyrococcus furiosus mixed RNA-seq
# > Haloferax volcanii mixed RNA-seq
# > Escherichia coli 

# > features were counted from bam files using featurecounts method (BAM files too large to upload)
# > stored as count tables with gene identifiers and rpm counts



#...................................combine ONT and ILLUMINA DATA FOR COMPARISON

#.................................ESCHERICHIA COLI

#...............................read in files
ecoli_ILLUMINA <- read_tsv(here("data/counts_data/ecoli_rpm_illumina.tsv")) 
ecoli_ONT      <- read_tsv(here("data/counts_data/ecoli_tex_ONT.tsv")) 

#...............................combine by gene_id and rename counts for two methods
ecoli_COUNTS   <- ecoli_ONT %>%
  mutate(gene = substr(GeneID,1,15)) %>%
  left_join(ecoli_ILLUMINA, by = "gene") %>%
  mutate(illumina_counts = as.integer(rpm),
         ont_counts = as.integer(counts)) %>%
  dplyr::filter(!is.na(ont_counts) == T & !is.na(illumina_counts) == T) %>%
  dplyr::select(ont_counts, illumina_counts,GeneID)

#...............................normalize using DESeq2 rlog
ecoli_abundace_correlation <- deseq2_estimation_and_rlog(ecoli_COUNTS)


#.................................PYROCOCCUS FURIOSUS

#...............................read in files
pfu_ILLUMINA <- read_tsv(here("data/counts_data/pfu_rpm_illumina.tsv")) 
pfu_ONT      <- read_tsv(here("data/counts_data/pfu_tex_ONT.tsv")) 

#...............................combine by gene_id and rename counts for two methods
pfu_COUNTS   <- pfu_ONT %>%
  mutate(gene = substr(GeneID,1,15)) %>%
  left_join(pfu_ILLUMINA, by = "gene") %>%
  mutate(illumina_counts = as.integer(rpm),
         ont_counts = as.integer(counts)) %>%
  dplyr::filter(!is.na(ont_counts) == T & !is.na(illumina_counts) == T) %>%
  dplyr::select(ont_counts, illumina_counts,GeneID)

#...............................normalize using DESeq2 rlog
pfu_abundace_correlation <- deseq2_estimation_and_rlog(pfu_COUNTS)


#.................................HALOFERAX VOLCANII

#...............................read in files
hvo_ILLUMINA <- read_tsv(here("data/counts_data/hvo_rpm_illumina.tsv")) 
hvo_ONT      <- read_tsv(here("data/counts_data/hvo_tex_ONT.tsv")) 

#...............................combine by gene_id and rename counts for two methods
hvo_COUNTS   <- hvo_ONT %>%
  mutate(gene = substr(GeneID,1,15)) %>%
  left_join(hvo_ILLUMINA, by = "gene") %>%
  mutate(illumina_counts = as.integer(rpm),
         ont_counts = as.integer(counts)) %>%
  dplyr::filter(!is.na(ont_counts) == T & !is.na(illumina_counts) == T) %>%
  dplyr::select(ont_counts, illumina_counts,GeneID)

#...............................normalize using DESeq2 rlog
hvo_abundace_correlation <- deseq2_estimation_and_rlog(hvo_COUNTS)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOT USING GGPLOT2
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................coloring for density
heat_color_npg <- c(pal_npg()(10)[4],
                    pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    pal_npg()(10)[5],
                    pal_npg()(10)[1],
                    pal_npg()(10)[8])

#...................................calculate densities
ecoli_abundace_correlation$density <- get_density(ecoli_abundace_correlation$ont_counts, ecoli_abundace_correlation$illumina_counts)
pfu_abundace_correlation$density <- get_density(pfu_abundace_correlation$ont_counts, pfu_abundace_correlation$illumina_counts)
hvo_abundace_correlation$density <- get_density(hvo_abundace_correlation$ont_counts, hvo_abundace_correlation$illumina_counts)

#...................................plot
correlation_expression_ecoli <- transcript_abundace_ggcorr(ecoli_abundace_correlation)
correlation_expression_pfu   <- transcript_abundace_ggcorr(pfu_abundace_correlation)
correlation_expression_hvo   <- transcript_abundace_ggcorr(hvo_abundace_correlation)

#...................................save plots
#.................................Supplementary Fig. 5a
pdf(here("figures/ONT_illumina_ecoli.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
correlation_expression_ecoli
dev.off()
#.................................Supplementary Fig. 5b
pdf(here("figures/ONT_illumina_pfu.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
correlation_expression_pfu
dev.off()
#.................................Supplementary Fig. 5c
pdf(here("figures/ONT_illumina_hvo.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
correlation_expression_hvo
dev.off()

