###########################################################################
###########################################################################
###
### COVERAGE DROP ANALYSIS 3´END
###
###########################################################################
###########################################################################


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
packages <- c("CoverageView", "tidyverse", "here", "ggthemes", 
              "data.table", "ggExtra", "Rsamtools", "GenomicAlignments", 
              "seqTools", "Rsubread", "ape", "DT", "ggpubr", "ggridges", "ggsci")
invisible(lapply(packages, require, character.only = TRUE))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#................................calculate enolase coverage
enolase_coverage_calculator <- function(input_bam_file, sequencing_set_name){
  as.data.table(coverage(input_bam_file)) %>%
    rownames_to_column("position") %>%
    dplyr::select(position, value) %>%
    mutate_at(vars(-position), funs(./sum(.))) %>%
    mutate(value = value * 10,
           position = 200 + (as.integer(position) - 1314)/10,
           group = "enolase [1314]",
           sequencing_set = sequencing_set_name) %>%
    dplyr::select(position, value, group, sequencing_set)
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

#...................................calculate coverage matrix for top_n genes in operon table (after featurecounts) using bam file
calculate_coverage_table <- function(bam_file, gene_table, operon_table, top_n, extend_size, step_size, filename, target_length_min, target_length_max){
  
  # > filter featourecounts output table for protein-coding genes, select top n abundant transcripts
  featurecounts_selected <- unlist(gene_table %>%
                                     dplyr::filter(type == "CDS", Length < target_length_max, Length > target_length_min) %>%
                                     arrange(desc(counts)) %>%
                                     head(n = top_n) %>%
                                     dplyr::select(GeneID))
  
  # > select transcriptional unit tables for top n genes
  terminating_positions <- operon_table %>%
    distinct(gene, .keep_all = T) %>%
    dplyr::filter(gene %in% featurecounts_selected)
  
  # > write bed.like file inclduding length of genes
  bedfile_info <- terminating_positions %>%
    mutate(chr = seqnames,
           start = as.integer(median_utr3),
           end = as.integer(median_utr3),
           name = gene,
           length = abs(median_utr3 - median_utr5)) %>%
    select(chr,start,end,strand,name, length_gene)
  
  # > bed.like file with standard info 
  bedfile_stop <- terminating_positions %>%
    mutate(chr = seqnames,
           start = as.integer(median_utr3),
           end = as.integer(median_utr3),
           name = gene) %>%
    select(chr,start,end,name)
  
  # > write bed.tables to files
  write.table(bedfile_info,
              file = paste(here("data/genome_data/"), filename, "_info.bed", sep = ""),row.names = F, col.names = F, quote = F, sep = "\t")  
  
  write.table(bedfile_stop,
              file = paste(here("data/genome_data/"), filename, ".bed", sep = ""),row.names = F, col.names = F, quote = F, sep = "\t")     
  
  bed_stop <- paste(here("data/genome_data/"), filename, ".bed", sep = "")
  
  # > calculate coverage using CoverageBamFile function
  trm <- CoverageBamFile(bam_file)
  
  # > generate coverage matrix extending "2000" nucleotdides on each side of the provided TTS
  coverage <- t(cov.matrix(trm,coordfile=bed_stop,extend=extend_size,num_cores=4, bin_width=step_size))
  
  # > filter coverage tables for strands - plus
  filter_plus <- coverage %>%
    as_tibble() %>%
    mutate(name = bedfile_stop$name) %>%
    left_join(bedfile_info) %>%
    filter(strand == "+") %>%
    dplyr::select(-name, -chr, -start,-end ,-strand, -length_gene) %>%
    t() %>%
    as_tibble() %>%
    set_colnames(bedfile_info$name[bedfile_info$strand == "+"]) %>%
    rownames_to_column(var = "position") %>%
    head(extend_size/step_size) %>%
    mutate_at(vars(-position), funs(./sum(.))) %>%
    mutate(position = as.numeric(position)) %>%
    gather(dataset,value, -position) 
  
  # > filter coverage tables for strands - minus
  filter_minus <- coverage %>%
    as_tibble() %>%
    mutate(name = bedfile_stop$name) %>%
    left_join(bedfile_info) %>%
    filter(strand == "-") %>%
    dplyr::select(-name, -chr, -start,-end ,-strand, -length_gene) %>%
    t() %>%
    as_data_frame() %>%
    set_colnames(bedfile_info$name[bedfile_info$strand == "-"]) %>%
    rownames_to_column(var = "position") %>%
    tail(extend_size/step_size) %>%
    mutate_at(vars(-position), funs(./sum(.))) %>%
    mutate(position = (extend_size/step_size):1) %>%
    arrange(position) %>%
    gather(dataset,value, -position) 
  
  # > combine strand specific coverage tables and group according to protein-coding gene lenghts
  filter_file <- rbind(filter_plus, filter_minus) %>%
    left_join(bedfile_info, by = c("dataset" = "name")) %>%
    mutate(group = ifelse(length_gene < 500, "[0-500]",
                          ifelse(length_gene < 1000 & length_gene >= 500, "[500-1000]",
                                 ifelse(length_gene < 1500 & length_gene >= 1000, "[1000-1500]", 
                                        ifelse(length_gene < 2000 & length_gene >= 1500, "[1500-2000]", 
                                               ifelse(length_gene >= 2000, "[>=2000]", NA))))))
  
  return(filter_file)
}



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CALCULATE COVERAGE
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................ECOLI

#................................load shinynano data
load(here("data/tidy_data/ecoli_tex_gene_table"))

#................................execute function
ecoli_coverage_table <- calculate_coverage_table(bam_file = here("data/mapped_data/ecoli_tex.bam"),
                                                 gene_table = full_gene_table,
                                                 operon_table = fread(here("data/operon_data/ecoli_tex_reads_for_operons.tsv")),
                                                 top_n  = 100, 
                                                 extend_size = 2000, 
                                                 step_size = 10, 
                                                 filename = "ecoli_top_100",
                                                 target_length_min = 0,
                                                 target_length_max = 6000)

#...................................PFU

#................................load shinynano data
load(here("data/tidy_data/pfu_tex_gene_table"))


#................................execute function
pfu_coverage_table <- calculate_coverage_table(bam_file = here("data/mapped_data/pfu_tex.bam"),
                                               gene_table = full_gene_table,
                                               operon_table = fread(here("data/operon_data/pfu_tex_reads_for_operons.tsv")),
                                               top_n  = 100, 
                                               extend_size = 2000, 
                                               step_size = 10, 
                                               filename = "pfu_top_100",
                                               target_length_min = 0,
                                               target_length_max = 6000)


#...................................HVO

#................................load shinynano data
load(here("data/tidy_data/hvo_tex_gene_table"))

#................................execute function
hvo_coverage_table <- calculate_coverage_table(bam_file = here("data/mapped_data/hvo_tex.bam"),
                                               gene_table = full_gene_table,
                                               operon_table = fread(here("data/operon_data/hvo_tex_reads_for_operons.tsv")),
                                               top_n  = 100, 
                                               extend_size = 2000, 
                                               step_size = 10, 
                                               filename = "hvo_top_100",
                                               target_length_min = 0,
                                               target_length_max = 6000)


#...................................COMBINE ALL TABLES
prokaryotes_coverage_table <- bind_rows(pfu_coverage_table %>%
                                          mutate(sequencing_set = "pfu"), 
                                        ecoli_coverage_table %>%
                                          mutate(sequencing_set = "ecoli"), 
                                        hvo_coverage_table %>%
                                          mutate(sequencing_set = "hvo"))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ADD ENOALSE DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#................................calculate
enolase_pfu  <- enolase_coverage_calculator(here("data/enolase_data/pfu_tex_enolase.bam"), "pfu")
enolase_coli <- enolase_coverage_calculator(here("data/enolase_data/ecoli_tex_enolase.bam"), "ecoli")
enolase_hvo  <- enolase_coverage_calculator(here("data/enolase_data/hvo_tex_enolase.bam"), "hvo")


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# COMBINE ALL
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#................................calculate
prokaryotes_coverage_table_plotting <- prokaryotes_coverage_table %>%
  dplyr::select(position, value, group, sequencing_set) %>%
  rbind(enolase_pfu, enolase_coli, enolase_hvo)

#................................color
heat_color_npg_control <- c("black",
                            pal_npg()(10)[4],
                            pal_npg()(10)[6],
                            pal_npg()(10)[7],
                            pal_npg()(10)[5],
                            pal_npg()(10)[1])

#................................grouping
prokaryotes_coverage_table_plotting$group <- factor(prokaryotes_coverage_table_plotting$group, 
                                                    levels(as.factor(prokaryotes_coverage_table_plotting$group))[c(6,1:5)])
prokaryotes_coverage_table_plotting$sequencing_set <- factor(prokaryotes_coverage_table_plotting$sequencing_set, 
                                                             levels(as.factor(prokaryotes_coverage_table_plotting$sequencing_set))[c(1,3,2)])


gg_decay_plots <- ggplot(data = prokaryotes_coverage_table_plotting, aes(x = position, y = value, color = group, fill = group, group = group)) +
  geom_smooth(span = 0.05, na.rm = TRUE, se = TRUE) +
  facet_grid(~sequencing_set) +
  theme_Publication_white() +
  scale_fill_manual(values = heat_color_npg_control) +
  scale_color_manual(values = heat_color_npg_control) +
  xlab("Postion to TTS [nt]") +
  ylab("normalized coverage") +
  scale_x_continuous(limits = c(0,200), expand = c(0,0), breaks = c(0,50, 150,200),labels = c("","-1500","-500", "TTS")) +
  theme(panel.grid.major.y = element_blank()) +
  guides(color = guide_legend(title = ""),
         fill = guide_legend(title = "")) 

pdf(here("figures/coverage_drop_3_prime.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_decay_plots
dev.off()


