###########################################################################
###########################################################################
###
### TRANSCRIPT ABUNDANCE COMPARISON ONT/ILLUMINA
###
###########################################################################
###########################################################################

## @knitr polya_tail_lengths

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
packages <- c("DESeq2", "ggpubr", "rtracklayer", "ggeconodist", "tidyverse", "here", "ggthemes", "data.table", "ggExtra", "Rsamtools", "GenomicAlignments", "seqTools", "Rsubread", "ape", "DT", "ggpubr", "ggridges", "ggsci")
invisible(lapply(packages, require, character.only = TRUE))

library(viridis)
library(tidyverse)
library(here)
library(ggthemes)
library(data.table)
library(ggExtra)
library(seqTools)
library(Rsamtools)
library(GenomicAlignments)
library(grid)
library(Rsubread)
library(ape)
library(DT)
library(ggsci)
library(ggridges)
library(ggpubr)
library(DESeq2)
library(heatmaply)
library(rtracklayer)
library(readxl)

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


#...................................publication_theme_white
theme_Publication_white <- function(base_size=14) {
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
# PFU
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#...................................for nanopore sequencing << #
load("/Volumes/Lacie/direct_rna_seq/data/shinynano_data/181219_pfu_all/Files/full_gene_table181219_pfu_all")

#...................................from illumina sequencing << #
gtf_file <- paste("/Volumes/Lacie/direct_rna_seq/data/genome_data/CP023154.gtf", sep = "")

#bamFiles_il <- "/Users/felixgrunberger/Documents/R/differential_0739/data/mapping_data/star/S305_01_52_1_S1_L008_R1_001_rna_trimmed/S305_01_52_1_S1_L008_R1_001_rna_trimmed_Aligned.sortedByCoord.out.bam"
bamFiles_il <- "/Volumes/TOSHIBA/backup_lacie_190320/annogesic_linux/ANNOgesic/input/BAMs/BAMs_map_reference_genomes/fragment/FRAG_trim.bam"
names(bamFiles_il) <- "illumina"
illumina_rnaseq_counts <- featureCounts(bamFiles_il,verbose = F,annot.ext = gtf_file, isGTFAnnotationFile = T, nthreads = 4,GTF.featureType = "CDS", allowMultiOverlap = F, isLongRead = F)$counts
illumina_rnaseq_counts_table <- illumina_rnaseq_counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  dplyr::rename(illumina_counts = 2) 

#......annotate with old names
old_new <- read_excel("/Users/felixgrunberger/Downloads/Table_2_Next Generation DNA-Seq and Differential RNA-Seq Allow Re-annotation of the Pyrococcus furiosus DSM 3638 Genome and Provide Insights Into Arch.XLSX", skip = 3, col_names = c("new", "start", "end", "old_r", "old_rs", "old_rse", "old", "new2", "name", "x")) %>%
  dplyr::select(new, start, end, old, name)

illumina_rnaseq_counts_table_rpm <- illumina_rnaseq_counts_table %>%
  mutate(total_reads = sum(illumina_counts)/1000000, rpm = round(illumina_counts / total_reads, digits = 1)) %>%
  left_join(old_new, by = c("gene" = "new")) %>%
  mutate(old_name = old, annotation = name) %>%
  dplyr::select(gene, rpm, old_name, annotation) %>%
  arrange(desc(rpm))

write_tsv(illumina_rnaseq_counts_table_rpm, path = "/Users/felixgrunberger/Desktop/rpm_pyrocoocus.tsv")

#...................................combine << #
nano_illu_counts_pfu <- full_gene_table %>%
  mutate(gene = substr(GeneID,1,15)) %>%
  left_join(illumina_rnaseq_counts_table, by = "gene") %>%
  dplyr::filter(!is.na(counts) == T & !is.na(illumina_counts) == T) %>%
  dplyr::select(counts, illumina_counts,GeneID)

export_table <- nano_illu_counts_pfu %>%
  dplyr::select(counts, GeneID) %>%
  mutate(GeneID = substr(GeneID,1,15)) %>%
  arrange(desc(counts)) %>%
  left_join(old_new, by = c("GeneID" = "new")) %>%
  mutate(Nanopore_counts = counts, old_name = old, annotation = name) %>%
  dplyr::select(Nanopore_counts, old_name, annotation, GeneID)

#...................................Calculate transcript abundances
colnames_data <- c("sequencing")
colData <- matrix(ncol = 1, nrow = 2)
colnames(colData) <- colnames_data 
colData <- data.table(colData) %>%
  mutate(sequencing = as.factor(c("nanopore", "illumina"))) %>%
  as_tibble()


dds_pfu <- DESeq2::DESeqDataSetFromMatrix(countData = nano_illu_counts_pfu[,1:2], 
                                          colData = colData,
                                          design = ~sequencing)

rownames(dds_pfu) <- nano_illu_counts_pfu$GeneID

# prefilter data set: remove rows (genes) with no counts
dds_pfu <- dds_pfu[ rowSums(counts(dds_pfu)) > 1, ]
rld <- rlogTransformation( dds_pfu ,blind = T)
helper_pfu <- as_tibble(assay(rld)[, 1:2]) 
vsd <- varianceStabilizingTransformation( dds_pfu ,blind = T)
helper_pfu2 <- as_tibble(assay(vsd)[, 1:2]) 

#...................................coloring
heat_color_npg <- c(pal_npg()(10)[4],
                    pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    pal_npg()(10)[5],
                    pal_npg()(10)[1],
                    pal_npg()(10)[8])

helper_pfu$density <- get_density(helper_pfu$counts, helper_pfu$illumina_counts)
helper_pfu2$density <- get_density(helper_pfu2$counts, helper_pfu2$illumina_counts)
#...................................plot
correlation_expression_pfu <- ggplot(data = helper_pfu2, aes(x = counts, y = illumina_counts, color = density)) +
  geom_abline(alpha = 0.5, linetype = "dashed", color = "black", size = 1) +
  scale_fill_gradientn(colours = heat_color_npg) +
  scale_color_gradientn(colours = heat_color_npg) +
  geom_point(alpha = 0.3, size = 2) +
  #stat_density2d(aes(alpha=..level.., color = ..level..),geom = "polygon", fill = NA, size = 0.3) + 
  theme_Publication_white() +
  guides(alpha = F) +
  guides(fill = guide_colorbar(barwidth = 5, barheight = 0.5, ticks = T, label = T),
         color = guide_colorbar(barwidth = 5, barheight = 0.5, ticks = T, label = T)) +
  xlab("MinION\n[rlog transformed expression]") +
  ylab("Illumina\n[rlog transformed expression]") +
  ggtitle("") +
  scale_x_continuous(limits = c(0,17)) + 
  scale_y_continuous(limits = c(0,17)) + 
  coord_equal() +
  stat_cor(method = "spearman", label.x = 0, label.y = 14,color = "black") 

pdf("/Users/felixgrunberger/Desktop/nanopore_native_rna/figures/raw_figures/190802_illumina_native_pfu.pdf", 
    width = 7, height = 7, paper = "special",onefile=FALSE)
correlation_expression_pfu
dev.off()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# HVO
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................for nanopore sequencing << #
load("/Volumes/Lacie/direct_rna_seq/data/shinynano_data/190219_hvo_wt_tex/Files/full_gene_table190219_hvo_wt_tex")

#...................................from illumina sequencing << #
gtf_file <- paste("/Volumes/Lacie/direct_rna_seq/data/genome_data/GCF_000025685.1_ASM2568v1_genomic.gff", sep = "")


bamFiles_il <- "/Volumes/Lacie/sra_re_analysis/SRP160422_hvo/data/mapped_data/SRR7811297_1_trimmed/SRR7811297_1_trimmed_sorted.bam"
names(bamFiles_il) <- "illumina"
illumina_rnaseq_counts <- featureCounts(bamFiles_il,verbose = F,annot.ext = gtf_file, isGTFAnnotationFile = T, nthreads = 4,GTF.featureType = "CDS",GTF.attrType = "ID", allowMultiOverlap = F, isLongRead = F)$counts
illumina_rnaseq_counts_table <- illumina_rnaseq_counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  dplyr::rename(illumina_counts = 2) 
glimpse(illumina_rnaseq_counts_table)

glimpse(full_gene_table)
#...................................combine << #
nano_illu_counts_hvo <- full_gene_table %>%
  left_join(illumina_rnaseq_counts_table, by = c("GeneID" = "gene")) %>%
  dplyr::filter(!is.na(counts) == T & !is.na(illumina_counts) == T) %>%
  dplyr::select(counts, illumina_counts,GeneID)

#...................................Calculate transcript abundances
colnames_data <- c("sequencing")
colData <- matrix(ncol = 1, nrow = 2)
colnames(colData) <- colnames_data 
colData <- data.table(colData) %>%
  mutate(sequencing = as.factor(c("nanopore", "illumina"))) %>%
  as_tibble()


dds_hvo <- DESeq2::DESeqDataSetFromMatrix(countData = nano_illu_counts_hvo[,1:2], 
                                          colData = colData,
                                          design = ~sequencing)

rownames(dds_hvo) <- nano_illu_counts_hvo$GeneID

# prefilter data set: remove rows (genes) with no counts
dds_hvo <- dds_hvo[ rowSums(counts(dds_hvo)) > 1, ]
rld <- rlogTransformation( dds_hvo ,blind = T)
helper_hvo <- as_tibble(assay(rld)[, 1:2]) 
helper_hvo2 <- as_tibble(assay(vsd)[, 1:2]) 

helper_hvo$GeneID <- rownames(assay(rld))
helper_hvo2$GeneID <- rownames(assay(vsd))

helper_hvo <- helper_hvo %>%
  left_join(full_gene_table, by = "GeneID")
helper_hvo2 <- helper_hvo2 %>%
  left_join(full_gene_table, by = "GeneID")

#...................................coloring
heat_color_npg <- c(pal_npg()(10)[4],
                    pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    pal_npg()(10)[5],
                    pal_npg()(10)[1],
                    pal_npg()(10)[8])

helper_hvo$density <- get_density(helper_hvo$counts.x, helper_hvo$illumina_counts)

#...................................plot
#correlation_expression_hvo <- 
correlation_expression_hvo <- ggplot(data = helper_hvo, aes(x = counts.x, y = illumina_counts, color = density, text = sprintf("locus: %s", locus_name))) +
  geom_abline(alpha = 0.5, linetype = "dashed", color = "black", size = 1) +
  scale_fill_gradientn(colours = heat_color_npg) +
  scale_color_gradientn(colours = heat_color_npg) +
  geom_point(alpha = 0.3, size = 2) +
  #stat_density2d(aes(alpha=..level.., color = ..level..),geom = "polygon", fill = NA, size = 0.3) + 
  theme_Publication_white() +
  guides(alpha = F) +
  guides(fill = guide_colorbar(barwidth = 5, barheight = 0.5, ticks = T, label = T),
         color = guide_colorbar(barwidth = 5, barheight = 0.5, ticks = T, label = T)) +
  xlab("MinION\n[rlog transformed expression]") +
  ylab("Illumina\n[rlog transformed expression]") +
  ggtitle("") +
  scale_x_continuous(limits = c(0,17)) + 
  scale_y_continuous(limits = c(0,17)) + 
  coord_equal() +
  stat_cor(method = "spearman", label.x = 0, label.y = 14,color = "black") 

pdf("/Users/felix/Desktop/nanopore_native_rna/figures/raw_figures/190802_illumina_native_hvo.pdf", 
    width = 7, height = 7, paper = "special",onefile=FALSE)
correlation_expression_hvo
dev.off()

ggplotly(correlation_expression_hvo, tooltip = "text")
ggplot(data = helper_hvo2, aes(x = counts.x, y = illumina_counts, color = abs(counts.x - illumina_counts) > 2 , fill =abs(counts.x - illumina_counts) > 2)) +
  geom_abline(alpha = 1, linetype = "dashed", color = "black", size = 1) +
  #scale_fill_gradientn(colours = heat_color_npg) +
  #scale_color_gradientn(colours = heat_color_npg) +
  geom_point(alpha = 1, size = 2) +
  #stat_density2d(aes(alpha=..level.., color = ..level..),geom = "polygon", fill = NA, size = 0.3) + 
  theme_Publication_white() +
  guides(alpha = F) +
  guides(fill = guide_colorbar(barwidth = 5, barheight = 0.5, ticks = T, label = T),
         color = guide_colorbar(barwidth = 5, barheight = 0.5, ticks = T, label = T)) +
  xlab("MinION\n[rlog transformed expression]") +
  ylab("Illumina\n[rlog transformed expression]") +
  ggtitle("") +
  scale_x_continuous(limits = c(0,17)) + 
  scale_y_continuous(limits = c(0,17)) + 
  coord_equal() +
  stat_cor(method = "spearman",color = "black") 

helper_hvo_biggesterror <- helper_hvo2 %>%
  dplyr::mutate(errorgroup = ifelse(abs(counts.x - illumina_counts) > 5, "big", "small")) %>%
  dplyr::filter(errorgroup == "big") %>%
  dplyr::select(locus_name, GeneID)
helper_hvo_biggesterror

library(ggsignif)



ggplot(data = helper_hvo_biggesterror, aes(x = errorgroup, y = counts.x)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("big", "small")), 
              map_signif_level=TRUE)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ECOLI
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................for nanopore sequencing << #
load("/Volumes/Lacie/direct_rna_seq/data/shinynano_data/190123_ecoli_all/Files/full_gene_table190123_ecoli_all")

#...................................from illumina sequencing << #
gtf_file <- paste("/Volumes/Lacie/direct_rna_seq/data/genome_data/e_coli_k12.gff", sep = "")


bamFiles_il <- "/Volumes/Lacie/sra_re_analysis/SRP056485_ecoli/data/mapped_data/SRR1927169_trimmed/SRR1927169_trimmed_sorted.bam"
names(bamFiles_il) <- "illumina"
illumina_rnaseq_counts <- featureCounts(bamFiles_il,verbose = F,annot.ext = gtf_file, isGTFAnnotationFile = T, nthreads = 4,GTF.featureType = "CDS",GTF.attrType = "ID", allowMultiOverlap = F, isLongRead = F)$counts
illumina_rnaseq_counts_table <- illumina_rnaseq_counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  dplyr::rename(illumina_counts = 2) 


#...................................combine << #
nano_illu_counts_ecoli <- full_gene_table %>%
  left_join(illumina_rnaseq_counts_table, by = c("GeneID" = "gene")) %>%
  dplyr::filter(!is.na(counts) == T & !is.na(illumina_counts) == T) %>%
  dplyr::select(counts, illumina_counts,GeneID)

#...................................Calculate transcript abundances
colnames_data <- c("sequencing")
colData <- matrix(ncol = 1, nrow = 2)
colnames(colData) <- colnames_data 
colData <- data.table(colData) %>%
  mutate(sequencing = as.factor(c("nanopore", "illumina"))) %>%
  as_tibble()


dds_ecoli <- DESeq2::DESeqDataSetFromMatrix(countData = nano_illu_counts_ecoli[,1:2], 
                                            colData = colData,
                                            design = ~sequencing)

rownames(dds_ecoli) <- nano_illu_counts_ecoli$GeneID

# prefilter data set: remove rows (genes) with no counts
dds_ecoli <- dds_ecoli[ rowSums(counts(dds_ecoli)) > 1, ]
rld <- rlogTransformation( dds_ecoli ,blind = T)
helper_ecoli <- as_tibble(assay(rld)[, 1:2]) 

#...................................coloring
heat_color_npg <- c(pal_npg()(10)[4],
                    pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    pal_npg()(10)[5],
                    pal_npg()(10)[1],
                    pal_npg()(10)[8])
helper_ecoli$density <- get_density(helper_ecoli$counts, helper_ecoli$illumina_counts)

#...................................plot
correlation_expression_ecoli <- ggplot(data = helper_ecoli, aes(x = counts, y = illumina_counts, color = density)) +
  geom_abline(alpha = 0.5, linetype = "dashed", color = "black", size = 1) +
  scale_fill_gradientn(colours = heat_color_npg) +
  scale_color_gradientn(colours = heat_color_npg) +
  geom_point(alpha = 0.3, size = 2) +
  #stat_density2d(aes(alpha=..level.., color = ..level..),geom = "polygon", fill = NA, size = 0.3) + 
  theme_Publication_white() +
  guides(alpha = F) +
  guides(fill = guide_colorbar(barwidth = 5, barheight = 0.5, ticks = T, label = T),
         color = guide_colorbar(barwidth = 5, barheight = 0.5, ticks = T, label = T)) +
  xlab("MinION\n[rlog transformed expression]") +
  ylab("Illumina\n[rlog transformed expression]") +
  ggtitle("") +
  scale_x_continuous(limits = c(0,17)) + 
  scale_y_continuous(limits = c(0,17)) + 
  coord_equal() +
  stat_cor(method = "spearman", label.x = 0, label.y = 14,color = "black") 

pdf("/Users/felix/Desktop/nanopore_native_rna/figures/raw_figures/190802_illumina_native_ecoli.pdf", 
    width = 7, height = 7, paper = "special",onefile=FALSE)
correlation_expression_ecoli
dev.off()



















how_many <- 20
helper$gene <- rownames(rld)

# NEW HEATMAP PLOT
table_expression_tex <- helper %>%
  arrange(desc(counts)) %>%
  head(n = how_many) %>%
  dplyr::select(gene, counts) %>%
  rename(counts = counts) %>%
  mutate(set = "nanopore (TEX)", gene = as.factor(gene))

table_expression_illumina <- helper %>%
  arrange(desc(counts)) %>%
  head(n = how_many) %>%
  dplyr::select(gene, illumina_counts) %>%
  rename(counts = illumina_counts) %>%
  mutate(set = "illumina")

table_expression_full <- rbind(table_expression_tex, table_expression_illumina) %>%
  mutate(gene = substr(gene,1,15))


gff_file <- "/Users/felixgrunberger/Documents/R/differential_0739/data/genome_data/CP023154.gff"
gff <- readGFF(file = gff_file) %>%
  as.tibble() %>%
  dplyr::filter(type == "CDS") %>%
  mutate(ID = substr(ID,1,15)) %>%
  dplyr::select(ID, product)

list_names <- list()
for (i in 1:length(table_expression_full$counts)){
  list_names <- c(list_names, gff$product[gff$ID == table_expression_full$gene[i]])
}
table_expression_full$gene <- paste(table_expression_full$gene, list_names, sep = ",")


wanted_order <- table_expression_full %>%
  dplyr::filter(set == "illumina") %>%
  arrange(desc(counts)) %>%
  dplyr::select(gene)

table_expression_full_new <- table_expression_full
table_expression_full_new$gene <-  factor(table_expression_full_new$gene, 
                                          levels = rev(unlist(wanted_order)))



heatmap_counts <- ggplot(table_expression_full_new, aes(x = set, y = gene)) + 
  geom_tile(aes(fill = counts), colour = "white", size = 0.1) + 
  scale_fill_gradientn(colours = heat_color_npg) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  coord_equal() +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 330, hjust = 0, colour = "grey50"))

heatmap_counts

pdf(here("/figures/190529_heatmap_genes_pfu.pdf"), 
    width = 14, height = 8, paper = "special",onefile=FALSE)
heatmap_counts
dev.off()
