###########################################################################
###########################################################################
###
### PLOT POLY(A) TAIL LENGHTS (OUTPUT FROM NANOPOLISH)
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

#...................................load data and combine with tidy data single read information
load_polya_and_combine <- function(input_polya, input_id_file, sequencing_set){
  
  # > load poly(A) estimation from nanopolish
  input_polya_filtered <- fread(input_polya) %>%
    as_tibble() %>%
    dplyr::filter(qc_tag == "PASS")
  
  # > load single read table and classify reads based on feature lengths
  id_table <- full_id_table %>%
    ungroup() %>%
    distinct(minion_read_name, .keep_all = TRUE) %>%
    mutate(length_gene = abs(start_gene - end_gene),
           fraction_read = (aligned_reads/length_gene*100)) %>%
    mutate(template_length = end_gene - start_gene,
           group = ifelse(template_length < 500, "[0-500]",
                          ifelse(template_length < 1000 & template_length >= 500, "[500-1000]",
                                 ifelse(template_length < 1500 & template_length >= 1000, "[1000-1500]", 
                                        ifelse(template_length < 2000 & template_length >= 1500, "[1500-2000]", 
                                               ifelse(template_length >= 2000, "[>=2000]", NA)))))) %>%
    ungroup() %>%
    mutate(mapped_type = ifelse(mapped_type == "rRNA", locus_name, mapped_type)) %>%
    left_join(input_polya_filtered, by = c("minion_read_name" = "readname")) %>%
    rbind() %>%
    dplyr::filter(!(is.na(polya_length))) %>%
    mutate(sequencing_set = sequencing_set) %>%
    dplyr::select(group, sequencing_set, polya_length, group, mapped_type, aligned_reads, identity, fraction_read)
  
  return(id_table)
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................paths to poly(A).tsv files from nanopolish output (gzipped)
polya_files <- paste(here("data/polya_data/"), list.files(here("data/polya_data/"),pattern = ".tsv.gz"), sep = "")

#...................................input: pre-calculated single-read tidy files
id_files <- paste(here("data/tidy_data/"), list.files(here("data/tidy_data/"),pattern = "_id_table"), sep = "")
# filter out lowsalt and ksga sample
id_files <- id_files[c(1,2,4,6,7,8)]

#...................................get sample names
sample_names <- unlist(lapply(polya_files, FUN=function(x){str_split_fixed(str_split_fixed(x, "_polya", 2)[1],"polya_data/",2)[2]}))

#...................................calculate polya tables
polya_table <- data.frame()

for (i in seq_along(sample_names)){
  load(id_files[i])
  polya_table <- rbindlist(list(polya_table, 
                                load_polya_and_combine(polya_files[i], full_id_table, sample_names[i])))
}

#...................................reorder levels
polya_table$sequencing_set <-  factor(polya_table$sequencing_set, 
                                      levels = (c("ecoli_tex","ecoli_notex",
                                                     "pfu_tex", "pfu_notex",
                                                     "hvo_tex", "hvo_notex")))

#...................................reorder groups (length of groups)
polya_table$group <- factor(polya_table$group, levels(as.factor(polya_table$group))[c(2,5,3,4,1)])

#...................................filter out tRNAs
polya_table <- polya_table %>%
  dplyr::filter(mapped_type != "tRNA")

#...................................set colors
heat_color_npg <- c(pal_npg()(10)[4],
                    pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    pal_npg()(10)[1])

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOT USING GGPLOT2
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................poly(A) tails for all sequencing sets grouped by feature lenghts (Supplementary Fig. 4)
gg_polya_tail <- ggplot(data = polya_table, aes(x = polya_length, y = group, color = mapped_type, fill = mapped_type)) +
  geom_density_ridges2(aes(height =..ndensity..), scale = 0.9, alpha = 0.3, size = 1) +
  scale_fill_npg() +
  facet_grid(~sequencing_set) +
  scale_y_discrete(expand = c(0,0,0.2,0)) +
  scale_color_npg() +
  theme_Publication_white() +
  guides(fill = guide_legend(title = ""), color = guide_legend(title = "")) +
  xlab("Poly(A)-tail (nt)") +
  ylab("Feature length (nt)") +
  scale_x_continuous(limits = c(0,300), expand = c(0,0), labels = c(0,100,200), breaks = c(0,100,200))

pdf(here("figures/polya_tail_analysis.pdf"), 
    width = 14, height = 7, paper = "special",onefile=FALSE)
gg_polya_tail
dev.off()

