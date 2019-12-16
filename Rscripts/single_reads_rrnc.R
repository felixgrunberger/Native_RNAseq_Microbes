###########################################################################
###########################################################################
###
### SINGLE-READ ANALYSIS - RRNAC ECOLI
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
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#................................load tidy data 
load(here("data/tidy_data/ecoli_tex_id_table"))

#................................filter for those reads mapping to rrna-c
full_id_table_spliced <- full_id_table %>%
  dplyr::filter(gene %in% c("rna-b3756", "rna-b3757", "rna-b3758", "rna-b3759", "rna-b3760", "rna-b3761"),start > 3939000 & end < 3946000) %>%
  distinct(minion_read_name, .keep_all = TRUE)

#................................read in bam file again to get cigar string - information is stored in "N" tag (skipped region from reference)
input_bam_file <- here("data/mapped_data/ecoli_tex.bam")
allReads <- readGAlignments(input_bam_file, use.names = T, param = ScanBamParam(tag=c("NM"), what="mapq"))
allReads_table <- GenomicAlignments::as.data.frame(allReads) %>%
  mutate(minion_read_name = names(allReads)) 

#................................filter reads
allReads_table_filtered <- allReads_table %>%
  dplyr::filter(minion_read_name %in% full_id_table_spliced$minion_read_name)

#................................calculate number of skipped reads / set keep if length > 0 / calculate total skipped
allReads_table_filtered$skipped_reads <- explodeCigarOpLengths(allReads_table_filtered$cigar, ops = "N")
allReads_table_filtered$keep <- lapply(allReads_table_filtered$skipped_reads, function(x) length(x)) > 0
allReads_table_filtered$total_skipped <- lapply(allReads_table_filtered$skipped_reads, function(x) sum(x))

#................................select for reads in the genomic region
allReads_table_filtered_region <- allReads_table_filtered %>%
  dplyr::filter(start > 3939000 & end < 3946000)

#................................extra set for spliced reads
allReads_table_filtered_skipped <- allReads_table_filtered_region %>%
  unnest(skipped_reads)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CATEGORIES (for rRNAc) - 40 longest reads in each category for plotting
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ecoli_gff <- read.gff(here("data/genome_data/ecoli.gff"))

ecoli_gff_region <-  ecoli_gff %>%
  dplyr::filter(type == "rRNA") %>%
  mutate(id_name = str_split_fixed(str_split_fixed(attributes, "ID=",2)[,2], ";Parent", 2)[,1]) %>%
  dplyr::select(id_name, start, end, strand) %>%
  dplyr::filter(id_name %in% c("rna-b3756", "rna-b3757", "rna-b3758", "rna-b3759", "rna-b3760", "rna-b3761"))


#................................MATURE 16S rRNA
category1 <- allReads_table_filtered_region %>%
  dplyr::filter(start > (ecoli_gff_region$start[1] + 10) & start < (ecoli_gff_region$start[1] + 25) & 
                  end < ecoli_gff_region$end[1] + 0, 
                keep == FALSE) %>%
  mutate(group = "category1", group_screen = "mature") %>%
  arrange(desc(width)) %>%
  head(40)

#................................MATURE 23S rRNA
category2 <- allReads_table_filtered_region %>%
  dplyr::filter(start > (ecoli_gff_region$start[2] + 10) & start < (ecoli_gff_region$start[2] + 20) & 
                  end < ecoli_gff_region$end[2] + 0, 
                keep == FALSE) %>%
  mutate(group = "category2", group_screen = "mature") %>%
  arrange(desc(width)) %>%
  head(40)

#................................MATURE  5S rRNA
category3 <- allReads_table_filtered_region %>%
  dplyr::filter(start > (ecoli_gff_region$start[3] + 10) & start < (ecoli_gff_region$start[3] + 20) & 
                  end < ecoli_gff_region$end[3] + 0 & end > ecoli_gff_region$end[3] - 20, 
                keep == FALSE) %>%
  mutate(group = "category3", group_screen = "mature") %>%
  arrange(desc(width)) %>%
  head(40)

#................................TSS TO AS LONG AS IT GETS
category4 <- allReads_table_filtered_region %>%
  dplyr::filter(start > (ecoli_gff_region$start[1] - 12 - 300 ) & start < (ecoli_gff_region$start[1]), 
                keep == FALSE) %>%
  mutate(group = "category4", group_screen = "tss_long") %>%
  arrange(desc(width)) %>%
  arrange(desc(width)) %>%
  head(40)

#................................BOTH RNASE-III CUTTING SITES - 16S
category5 <- allReads_table_filtered_region %>%
  dplyr::filter(start < ecoli_gff_region$start[1] - 12 & start > ecoli_gff_region$start[1] - 120, 
                end < ecoli_gff_region$end[1] + 40 & end > ecoli_gff_region$end[1] + 20, 
                keep == FALSE) %>%
  mutate(group = "category5", group_screen = "rnaseIII_16") %>%
  arrange(desc(width)) %>%
  dplyr::filter(width > 1650) %>%
  mutate(distance = end-start) %>%
  arrange(desc(width)) %>%
  head(40)

#................................BOTH RNASE-III CUTTING SITES - 23S
category6 <- allReads_table_filtered_region %>%
  dplyr::filter(start > ecoli_gff_region$start[2] - 0 & start < ecoli_gff_region$start[2] +6, 
                end < ecoli_gff_region$end[2] + 10 & end > ecoli_gff_region$end[2]+2, 
                keep == FALSE) %>%
  mutate(group = "category6", group_screen = "rnaseIII_23") %>%
  arrange(desc(width)) %>%
  head(40)


#................................combine datasets
all_cats <- rbindlist(list(category1, category2, category3, category4, category5, category6), fill = T) %>%
  arrange(desc(start)) %>%
  group_by(group) %>%
  arrange(desc(width)) %>%
  ungroup() %>%
  mutate(lfd_n = 1:n()) 

#....plot
pdf("/Users/felixgrunberger/Documents/screen/nanopore_native_rna/figures/raw_figures/191205_singlereads_ecoli.pdf",
    width = 16, height = 10, paper = "special",onefile=FALSE)
ggplot(data = all_cats) +
  geom_segment(aes(x = start - 12, xend = end, yend = lfd_n, y = lfd_n, color = group_screen), size = 0.7) +
  theme_Publication_white() +
  geom_vline(xintercept = c(ecoli_gff_region$start[1],ecoli_gff_region$end[1],
                            ecoli_gff_region$start[2],ecoli_gff_region$end[2],
                            ecoli_gff_region$start[3],ecoli_gff_region$end[3]), linetype = "dashed", alpha = 0.5) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_blank()) 

dev.off()
