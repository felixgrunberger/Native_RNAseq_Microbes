###########################################################################
###########################################################################
###
### TU DETECTED COMPARISON TO DATABASE
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


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................ECOLI

#................................database
operons_publication <- readxl::read_excel(here("data/operon_data/ecoli_database_operons.xlsx")) %>%
  rowwise() %>%
  mutate(GenesInTUC = gsub(GenesInTUC,replacement = ", ", pattern = ";"))

#................................edit names 
ecoli_gff_table_annotation <- read.gff(here("data/genome_data/ecoli.gff")) %>%
  dplyr::filter(type == "CDS") %>%
  mutate(id_name = str_split_fixed(str_split_fixed(attributes, "ID=",2)[,2], ";Parent", 2)[,1],
         parent = str_split_fixed(str_split_fixed(attributes, "Parent=gene-",2)[,2], ";Dbxref", 2)[,1]) %>%
  mutate(gene = substring(id_name, 1, 12)) %>%
  dplyr::select(gene, parent)

#................................load Nanopore detected subTUs
sub_tus_ecoli   <- readxl::read_xlsx(here("tables/tu_tables/tu_ecoli.xlsx"))
ecoli_collapsed <- readxl::read_xlsx(here("tables/tu_tables/tu_ecoli.xlsx")) %>%
  dplyr::rename(operon_number = operon_id, genes = genes_in_operon) %>%
  mutate(operon_number = 1:length(sub_tus_ecoli$start_operon)) %>%
  dplyr::select(operon_number, genes) %>%
  mutate(genes = strsplit(genes, ",")) %>% 
  unnest(genes) %>%
  rowwise() %>%
  mutate(genes = str_split_fixed(genes, "\\.",2)[1]) %>%
  left_join(ecoli_gff_table_annotation, by = c("genes" = "gene")) %>%
  group_by(operon_number) %>%
  mutate(genes = toString(parent)) %>%
  distinct(operon_number, genes) %>%
  mutate(size_operon = 1 + str_count(genes, ","))

#................................combine
ecoli_full_tu <- full_join(operons_publication,ecoli_collapsed, by = c("GenesInTUC" = "genes"))

#...................................PFU

#................................database
operon_prediction_pfu <- fread(here("data/operon_data/pfu_database_operons.tsv")) %>%
  dplyr::rename(operon_number = 1, genes = 2) %>%
  mutate(genes = strsplit(as.character(genes), " ")) %>% 
  unnest(genes) %>%
  mutate(genes = str_split_fixed(genes, ":",2)[,2]) %>%
  group_by(operon_number) %>%
  mutate(genes = toString(genes)) %>%
  distinct(operon_number, genes)

#................................load Nanopore detected subTUs
sub_tus_pfu <- readxl::read_xlsx(here("tables/tu_tables/tu_pfu.xlsx"))
pfu_collapsed <- readxl::read_xlsx(here("tables/tu_tables/tu_pfu.xlsx")) %>%
  dplyr::rename(operon_number = operon_id, genes = genes_in_operon) %>%
  ungroup() %>%
  mutate(operon_number = 1:length(sub_tus_pfu$start_operon)) %>%
  dplyr::select(operon_number, genes) %>%
  mutate(genes = strsplit(genes, ",")) %>% 
  unnest(genes) %>%
  mutate(genes = str_split_fixed(genes, ".p",2)[,1]) %>%
  group_by(operon_number) %>%
  mutate(genes = toString(genes)) %>%
  distinct(operon_number, genes) %>%
  mutate(size_operon = 1 + str_count(genes, ","))

#................................combine
pfu_full_tu <- full_join(operon_prediction_pfu,pfu_collapsed, by = c("genes"))

#...................................HVO

#................................database
gff_table_combined <- fread(here("data/genome_data/hvo_annotation.tsv"))

operon_prediction_hvo <- fread(here("data/operon_data/hvo_database_operons.opr")) %>%
  dplyr::left_join(gff_table_combined, by = c("Synonym" = "old_name")) %>%
  dplyr::rename(operon_number = 1, genes = id_name) %>%
  dplyr::select(operon_number, genes) %>%
  group_by(operon_number) %>%
  mutate(genes = toString(genes)) %>%
  distinct(operon_number, genes)

#................................load Nanopore detected subTUs
sub_tus_hvo <- readxl::read_xlsx(here("tables/tu_tables/tu_hvo.xlsx"))
hvo_collapsed <- readxl::read_xlsx(here("tables/tu_tables/tu_hvo.xlsx")) %>%
  ungroup() %>%
  dplyr::rename(operon_number = operon_id, genes = genes_in_operon) %>%
  mutate(operon_number = 1:length(sub_tus_hvo$start_operon)) %>%
  dplyr::select(operon_number, genes) %>%
  mutate(genes = strsplit(genes, ",")) %>% 
  unnest(genes) %>%
  mutate(genes = str_split_fixed(genes, ".p",2)[,1]) %>%
  group_by(operon_number) %>%
  mutate(genes = toString(genes)) %>%
  distinct(operon_number, genes) %>%
  mutate(size_operon = 1 + str_count(genes, ","))

#................................combine
hvo_full_tu <- full_join(operon_prediction_hvo,hvo_collapsed, by = c("genes"))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# UPSETR GROUP COMPARISON
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................calculate intersects (off all possible variations) and save in new list
upsetr_table <- c(`ECOLI_ONT` = nrow(ecoli_full_tu) - nrow(operons_publication), 
                  `ECOLI_ILLUMINA` = nrow(ecoli_full_tu) - nrow(ecoli_collapsed), 
                  `ECOLI_ONT&ECOLI_ILLUMINA` = nrow(operons_publication) + nrow(ecoli_collapsed) - nrow(ecoli_full_tu),
                  `PFU_ONT` = nrow(pfu_full_tu) - nrow(operon_prediction_pfu), 
                  `PFU_ILLUMINA` = nrow(pfu_full_tu) - nrow(pfu_collapsed), 
                  `PFU_ONT&PFU_ILLUMINA` = nrow(operon_prediction_pfu) + nrow(pfu_collapsed) - nrow(pfu_full_tu),
                  `HVO_ONT` = nrow(hvo_full_tu) - nrow(operon_prediction_hvo), 
                  `HVO_ILLUMINA` = nrow(hvo_full_tu) - nrow(hvo_collapsed), 
                  `HVO_ONT&HVO_ILLUMINA` = nrow(operon_prediction_hvo) + nrow(hvo_collapsed) - nrow(hvo_full_tu))

#...................................plot using upsetR function (Supplementary Fig. 9a)
tu_upset <- upset(data = fromExpression(upsetr_table), 
                  sets = rev(c("ECOLI_ONT", "ECOLI_ILLUMINA", "PFU_ONT", "PFU_ILLUMINA", "HVO_ONT", "HVO_ILLUMINA")),
                  keep.order = TRUE, order.by = c("degree"),
                  sets.x.label = "Number of TUs", mainbar.y.label = "Intersections", 
                  line.size = 3, point.size = 11, 
                  text.scale = c(1.6, 1.2, 1.2, 1.2, 2, 1.6), 
                  main.bar.color = pal_npg()(10)[c(4,6,7,4,4,6,6,7,7)], 
                  sets.bar.color = "gray60",
                  matrix.color = pal_npg()(10)[4])

#...................................save
pdf(here("figures/tu_total_comparison.pdf"), 
    width = 14, height = 14, paper = "special",onefile=FALSE)
tu_upset
dev.off()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# COMPARISON OF TU IN ECOLI WITH DATABASE
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................database ecoli 
operons_publication_bigger <- operons_publication %>%
  dplyr::filter(NumberofTUInTUC > 0) %>%
  mutate(set = "database") %>%
  dplyr::rename(tu_size = NumberofTUInTUC) %>%
  dplyr::select(set, tu_size)

#...................................Nanopore detected TU ecoli
operons_coli_shiny <- ecoli_collapsed %>%
  ungroup() %>%
  dplyr::filter(size_operon > 0) %>%
  mutate(set = "nanopore") %>%
  dplyr::rename(tu_size = size_operon) %>%
  dplyr::select(set, tu_size)

#...................................combine tables
operons_comparison <- rbind(operons_publication_bigger, operons_coli_shiny) %>%
  dplyr::filter(!is.na(tu_size)) %>%
  mutate(tu_size = (tu_size)) %>%
  group_by(tu_size, set) %>%
  mutate(n_operons = n()) %>%
  distinct(set, tu_size, n_operons)

#...................................set colors
color_npg2 <- c(pal_npg()(10)[7],
                pal_npg()(10)[4]) 

#...................................plot distribution (Supplementary Fig. 9b)
gg_tu_ecoli_comparison <- ggplot(data = operons_comparison, aes(x = n_operons, y = as.factor(tu_size), color = set)) +
  geom_point(size = 5, alpha = 0.7) +
  theme_Publication_white() +
  scale_color_manual(values = color_npg2) +
  theme(panel.grid.major.x = element_blank()) +
  ylab("Number of genes in TU") +
  xlab("Counts") +
  guides(fill = guide_legend(title = ""), color = guide_legend(title = ""))

pdf(here("figures/tu_ecoli_comparison.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_tu_ecoli_comparison
dev.off()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# plot number of genes in TUC in all three
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................combine 
all_operons_collapsed <- bind_rows(ecoli_collapsed %>%
                                     ungroup() %>%
                                     dplyr::select(size_operon) %>%
                                     mutate(sequencing_set = "Ecoli\n(TEX)"),
                                   pfu_collapsed %>%
                                     ungroup() %>%
                                     dplyr::select(size_operon) %>%
                                     mutate(sequencing_set = "Pfu\n(TEX)"),
                                   hvo_collapsed %>%
                                     ungroup() %>%
                                     dplyr::select(size_operon) %>%
                                     mutate(sequencing_set = "Hvo\n(TEX)"))

#...................................filter out non existing operons
all_operons_collapsed_multi <- all_operons_collapsed %>%
  dplyr::filter(size_operon > 0)

#...................................group by size of the operon and the sample
all_operons_comparison_grouped <- all_operons_collapsed_multi %>%
  mutate(size_operon = as.factor(size_operon)) %>%
  group_by(size_operon, sequencing_set) %>%
  mutate(n_operons = n()) %>%
  distinct(sequencing_set, size_operon, n_operons)

#...................................set colors
color_npg3 <- c(pal_npg()(10)[4],
                pal_npg()(10)[7],
                pal_npg()(10)[1]) 

#...................................plot counts of detected TU sizes (Supplementary Fig. 9c)
gg_operons_samples <- ggplot(data = all_operons_comparison_grouped, aes(x = n_operons, y = size_operon, color = sequencing_set, fill = sequencing_set)) +
  geom_point(size = 5, alpha = 0.8, shape = 21) +
  theme_Publication_white() +
  scale_color_manual(values = color_npg3) +
  scale_fill_manual(values = color_npg3) +
  theme(panel.grid.major.x = element_blank()) +
  ylab("Number of genes in TU") +
  xlab("counts") +
  xlab("n") +
  guides(fill = guide_legend(title = ""), color = guide_legend(title = "")) 

pdf(here("figures/tu_comparison_samples.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_operons_samples
dev.off()




