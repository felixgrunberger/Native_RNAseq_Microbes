###########################################################################
###########################################################################
###
### TRANSCRIPTION TERMINATION SITES
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

#...................................coloring heatmaps
heat_color_npg <- c(pal_npg()(10)[4],
                    pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    pal_npg()(10)[5],
                    pal_npg()(10)[1],
                    pal_npg()(10)[8])

#...................................subsampling intergenic regions
sampleString = function(string) {
  nStart = sample(1:(nchar(string) - 91),1)
  substr(string, nStart, nStart + 90)
}

#...................................plot rna terminator motif
plot_rna_motif <- function(input_motif){
  ggplot() + 
    geom_logo(input_motif, font = "helvetica_bold", col_scheme = color_scale, seq_type = "rna") + 
    theme_logo() +
    theme_Publication_white() +
    theme(panel.grid.major = element_line(colour = NA),
          axis.ticks.x = element_line(colour = NA), 
          axis.text.x = element_text(size = 0)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(limits = c(0,2), expand = c(0,0))
}

#...................................plot rna terminator heatmap
plot_rna_position <- function(position_table){
  numbers <- list()
  for (i in 1:length(position_table$V1)){
    numbers[i] <- list(position_table$start_motif_from_TTS[i]:position_table$end_motif_from_TTS[i])
  }
  
  heat_table <- as.data.frame(table(unlist(numbers))) %>%
    mutate(coordinate = as.numeric(as.character(Var1)))
  
  ggplot(heat_table, aes(x = coordinate, color = Freq, fill = Freq, y = as.factor(1))) +
    geom_tile(size = 0.5) +
    scale_fill_gradientn(colours = heat_color_npg) +
    scale_color_gradientn(colours = heat_color_npg) +
    scale_x_continuous(limits = c(-45,45),expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_Publication_white() +
    xlab("Position from TTS (nt)") +
    ylab("") +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "white") +
    guides(color = guide_colorbar(title = "counts",barwidth = 15, barheight = 0.5, ticks = T, label = T)) +
    guides(fill = F)
}

#...................................make consensus matrix
make_matrix_c <- function(input_term_sequence){
  consensusMatrix(input_term_sequence, as.prob = T) %>%
  t() %>%
  as_tibble() %>%
  mutate(position = c(-45:45)) %>%
  gather(key = position) %>%
  dplyr::rename(base = 1) %>%
  mutate(position = rep(c(-45:45),4))
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

#..............................color scale
color_npg4 <- c(pal_npg()(10)[1],
                pal_npg()(10)[2],
                pal_npg()(10)[10],
                pal_npg()(10)[3])

#...................................nucleotide enrichment plotting
nucleotide_enrichment_plotting <- function(input_matrix){
  ggplot(data = input_matrix, aes(x = position, y = log_value, color = base)) +
    geom_hline(yintercept = 0, alpha = 0.5, linetype = "dashed") +
    geom_line(size = 2.5, alpha = 1) +
    xlab("Position relative to 3´end (nt)") +
    scale_y_continuous(limits = c(-1.3,1.3),breaks = c(-0.5,0,0.5)) +
    ylab("Nucleotide enrichment \n(log2-fold enrichment)") +
    theme_Publication_white() +
    scale_color_manual(values = color_npg4)
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD & TIDY DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................load transcriptional unit annotation
filtered_ids  <- paste(here("data/operon_data/"), list.files(here("data/operon_data/"),pattern = "for_operons.tsv"), sep = "")
collapsed_ids <- paste(here("data/operon_data/"), list.files(here("data/operon_data/"),pattern = "tex_operons.tsv"), sep = "")

#.................................is a gene start or end of an transcriptional unit?
first_gene_ecoli <- type_of_gene_in_operon(fread(collapsed_ids[1]), "first")
last_gene_ecoli  <- type_of_gene_in_operon(fread(collapsed_ids[1]), "last")
first_gene_pfu   <- type_of_gene_in_operon(fread(collapsed_ids[4]), "first")
last_gene_pfu    <- type_of_gene_in_operon(fread(collapsed_ids[4]), "last")
first_gene_hvo   <- type_of_gene_in_operon(fread(collapsed_ids[3]), "first")
last_gene_hvo    <- type_of_gene_in_operon(fread(collapsed_ids[3]), "last")

#.................................load fasta files
fasta_file  <- paste(here("data/genome_data/"), list.files(here("data/genome_data/"),pattern = ".fasta"), sep = "")
ecoli_fasta <- readDNAStringSet(fasta_file[2])
pfu_fasta   <- readDNAStringSet(fasta_file[9])
hvo_fasta   <- readDNAStringSet(fasta_file[5])

#.................................get TTS and get sequence from -45 to +45

#...............................E. COLI
ecoli_utr3 <- fread(filtered_ids[1]) %>%
  dplyr::filter(mapped_type == "CDS") %>%
  left_join(first_gene_ecoli, by = "gene") %>%
  left_join(last_gene_ecoli, by = "gene") %>%
  group_by(gene) %>%
  mutate(end_operon = as.character(end_operon),
         utr3_length= ifelse(strand == "+", median_utr3 - end_gene, start_gene - median_utr3),
         group = ifelse(is.na(end_operon) == T, "rest", "end"))  %>%
  dplyr::select(median_utr3, gene, strand, utr3_length, start_operon, end_operon, group) %>%
  dplyr::filter(!is.na(median_utr3), group == "end") %>%
  distinct(gene, .keep_all = T) %>%
  rowwise() %>%
  mutate(terminator_sequence = ifelse(strand == "+", as.character(ecoli_fasta$`U00096.2 Escherichia coli str. K-12 substr. MG1655, complete genome`[(median_utr3 - 45):(median_utr3 + 45)]),
                                      as.character(reverseComplement(ecoli_fasta$`U00096.2 Escherichia coli str. K-12 substr. MG1655, complete genome`[(median_utr3 - 45):(median_utr3 + 45)])))) 

#...............................P. FURIOSUS
pfu_utr3 <- fread(filtered_ids[4]) %>%
  dplyr::filter(mapped_type == "CDS") %>%
  left_join(first_gene_pfu, by = "gene") %>%
  left_join(last_gene_pfu, by = "gene") %>%
  group_by(gene) %>%
  mutate(end_operon = as.character(end_operon),
         utr3_length= ifelse(strand == "+", median_utr3 - end_gene, start_gene - median_utr3),
         group = ifelse(is.na(end_operon) == T, "rest", "end"))  %>%
  dplyr::select(median_utr3, gene, strand, utr3_length, start_operon, end_operon, group) %>%
  dplyr::filter(!is.na(median_utr3), group == "end") %>%
  distinct(gene, .keep_all = T) %>%
  rowwise() %>%
  mutate(terminator_sequence = ifelse(strand == "+", as.character(pfu_fasta$CP023154[(median_utr3 - 45):(median_utr3 + 45)]),
                                      as.character(reverseComplement(pfu_fasta$CP023154[(median_utr3 - 45):(median_utr3 + 45)])))) 

#...............................H. VOLCANII
#.............................include genome information to include chr information
hvo_gff_table_cds_chr <- read.gff(here("data/genome_data/hvo.gff")) %>%
  as_tibble() %>%
  dplyr::filter(type %in% "CDS") %>%
  mutate(ID = str_split_fixed(str_split_fixed(attributes, "ID=",2)[,2], ";Parent", 2)[,1],
         GeneID = str_split_fixed(str_split_fixed(attributes, "GeneID:", 2)[,2], ";Name=",2)[,1]) %>%
  dplyr::rename(strand_gene = strand) %>%
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

hvo_utr3 <- fread(filtered_ids[3]) %>%
  dplyr::filter(mapped_type == "CDS") %>%
  left_join(hvo_gff_table_cds_chr, by = c("gene" = "ID")) %>%
  left_join(first_gene_hvo, by = "gene") %>%
  left_join(last_gene_hvo, by = "gene") %>%
  group_by(gene) %>%
  mutate(end_operon = as.character(end_operon),
         utr3_length= ifelse(strand == "+", median_utr3 - end_gene, start_gene - median_utr3),
         group = ifelse(is.na(end_operon) == T, "rest", "end"))  %>%
  dplyr::select(median_utr3, gene, strand, utr3_length, start_operon, end_operon, group, seqid, chr_size) %>%
  dplyr::filter(!is.na(median_utr3), group == "end") %>%
  distinct(gene, .keep_all = T) %>%
  rowwise() %>%
  dplyr::filter(median_utr3 > 45, median_utr3 < chr_size) %>%
  mutate(terminator_sequence = ifelse(strand == "+" & seqid == "NC_013964.1", as.character(hvo_fasta$`NC_013964.1 Haloferax volcanii DS2 plasmid pHV3, complete sequence`[(median_utr3 - 45):(median_utr3 + 45)]),
                                  ifelse(strand == "+" & seqid == "NC_013965.1", as.character(hvo_fasta$`NC_013965.1 Haloferax volcanii DS2 plasmid pHV2, complete sequence`[(median_utr3 - 45):(median_utr3 + 45)]),
                                    ifelse(strand == "+" & seqid == "NC_013966.1", as.character(hvo_fasta$`NC_013966.1 Haloferax volcanii DS2 plasmid pHV4, complete sequence`[(median_utr3 - 45):(median_utr3 + 45)]),
                                           ifelse(strand == "+" & seqid == "NC_013967.1", as.character(hvo_fasta$`NC_013967.1 Haloferax volcanii DS2, complete genome`[(median_utr3 - 45):(median_utr3 + 45)]),
                                                  ifelse(strand == "+" & seqid == "NC_013968.1", as.character(hvo_fasta$`NC_013968.1 Haloferax volcanii DS2 plasmid pHV1, complete sequence`[(median_utr3 - 45):(median_utr3 + 45)]),
                                                         ifelse(strand == "-" & seqid == "NC_013968.1", as.character(reverseComplement(hvo_fasta$`NC_013968.1 Haloferax volcanii DS2 plasmid pHV1, complete sequence`[(median_utr3 - 45):(median_utr3 + 45)])),
                                                                ifelse(strand == "-" & seqid == "NC_013967.1", as.character(reverseComplement(hvo_fasta$`NC_013967.1 Haloferax volcanii DS2, complete genome`[(median_utr3 - 45):(median_utr3 + 45)])),
                                                                       ifelse(strand == "-" & seqid == "NC_013966.1", as.character(reverseComplement(hvo_fasta$`NC_013966.1 Haloferax volcanii DS2 plasmid pHV4, complete sequence`[(median_utr3 - 45):(median_utr3 + 45)])),
                                                                            ifelse(strand == "-" & seqid == "NC_013965.1", as.character(reverseComplement(hvo_fasta$`NC_013965.1 Haloferax volcanii DS2 plasmid pHV2, complete sequence`[(median_utr3 - 45):(median_utr3 + 45)])),
                                                                              ifelse(strand == "-" & seqid == "NC_013964.1", as.character(reverseComplement(hvo_fasta$`NC_013964.1 Haloferax volcanii DS2 plasmid pHV3, complete sequence`[(median_utr3 - 45):(median_utr3 + 45)]))))))))))))) %>%
  dplyr::select(median_utr3, gene, strand, utr3_length, start_operon, end_operon, group, terminator_sequence)


#.................................write terminator sequences from -45 to +45 to fasta file
write.fasta(as.list(ecoli_utr3$terminator_sequence), 
            as.list(ecoli_utr3$gene), 
            file = here("data/meme_data/terminator_sequences_ecoli_operon_ending.fasta")) 

write.fasta(as.list(pfu_utr3$terminator_sequence), 
            as.list(pfu_utr3$gene), 
            file = here("data/meme_data/terminator_sequences_pfu_operon_ending.fasta")) 

write.fasta(as.list(hvo_utr3$terminator_sequence), 
            as.list(hvo_utr3$gene), 
            file = here("data/meme_data/terminator_sequences_hvo_operon_ending.fasta")) 

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 3´UTR ANALYSIS (lengths)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#.................................comparison of ecoli data to term seq set (combine with gff)
termseq_ecoli <- read_excel(<path to term seq study in e coli>,skip = 11) %>%
  mutate(short_name = substr(`Locus tag`, 9, 12), start = `gene fr`) %>%
  mutate(CDS_end = ifelse(`Gene strand` == "+", `gene to`, `gene fr`)) %>%
  mutate(utr3_length =    ifelse(`Gene strand` == "+", `primary 3' end position` - CDS_end, CDS_end - `primary 3' end position`)) %>%
  mutate(sequencing_set = "ecoli", seq = "Illumina") %>%
  dplyr::select(sequencing_set, utr3_length, seq)

#.................................combine all
utr3_all <- rbind(ecoli_utr3 %>%
                    mutate(sequencing_set = "ecoli", seq = "ONT"), 
                  pfu_utr3 %>%
                    mutate(sequencing_set = "pfu", seq = "ONT"), 
                  hvo_utr3 %>%
                    mutate(sequencing_set = "hvo", seq = "ONT")) %>%
  dplyr::select(utr3_length, sequencing_set, seq) %>%
  rbind(termseq_ecoli)

#.................................reorder levels
utr3_all$sequencing_set <-  factor(utr3_all$sequencing_set, 
                                   levels = rev(c("ecoli", "pfu", "hvo")))

#.................................set colors
color_npg2 <- c(pal_npg()(10)[7],
                pal_npg()(10)[4])

#.................................plot 3´UTR lengths (ONT only & comparison to termseq for ecoli, Fig. 2c)
gg_utr3 <- ggplot(data = utr3_all, aes(x = utr3_length, y = sequencing_set, fill = seq, color = seq)) +
  geom_density_ridges2(alpha = 0.4, size = 1, scale = 1, color = NA) +
  geom_density_ridges2(alpha = 0.8, size = 0.2, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 100, color = "black") +
  theme_Publication_white() +
  scale_color_manual(values = color_npg2) +
  scale_fill_manual(values = color_npg2) +
  scale_x_continuous(limits = c(-50,500)) +
  scale_y_discrete(expand = c(0,0)) +
  ylab("") +
  xlab("3´ UTR length (nt)")

pdf(here("figures/utr3_lengths.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_utr3
dev.off()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MOTIF ANALYSIS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#.................................MEME command (outside R in console)
# meme >>INPUTFILE<< -dna -oc >>OUTPUTDIRECTORY<< -nostatus -nostatus -time 18000 -maxsize 600000 -mod zoops -nmotifs 5 -minw 3 -maxw 50 -bfile >>CUSTOM_BGFILE_INTGERGENICREGIONS<<

#.................................save found motifs to fasta (from html file) and plot motif in R

#...............................set colors for ATCG
color_scale = make_col_scheme(chars=c('A', 'U', 'C', 'G'), 
                              cols=pal_npg()(10)[c(1,3,2,10)])

#...............................load sequence output from MEME 
motif_ecoli1 <- read.table(here("data/meme_data/ecoli_motif1_terminator.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(!grepl("offset", V1)) %>%
  mutate(V1 = gsub("T", "U", V1))

motif_ecoli2 <- read.table(here("data/meme_data/ecoli_motif2_terminator.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(!grepl("offset", V1)) %>%
  mutate(V1 = gsub("T", "U", V1))

motif_pfu <- read.table(here("data/meme_data/pfu_motif_terminator.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(!grepl("offset", V1)) %>%
  mutate(V1 = gsub("T", "U", V1))

motif_hvo <- read.table(here("data/meme_data/hvo_motif_terminator.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(!grepl("offset", V1)) %>%
  mutate(V1 = gsub("T", "U", V1))

#...............................load sequence position output from MEME 
position_ecoli1 <- read.table(here("data/meme_data/ecoli_motif1_terminator.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(grepl("offset", V1)) %>%
  mutate(start_motif = as.numeric(str_split_fixed(V1, "offset=",2)[,2]) + 1,
         end_motif = start_motif + 30,
         start_motif_from_TTS = -47 + start_motif,
         end_motif_from_TTS = -47 + end_motif,
         gene = str_split_fixed(str_split_fixed(V1, ">", 2)[,2], "_site_1", 2)[,1],
         group = "motif1") 

position_ecoli2 <- read.table(here("data/meme_data/ecoli_motif2_terminator.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(grepl("offset", V1)) %>%
  mutate(start_motif = as.numeric(str_split_fixed(V1, "offset=",2)[,2]) + 1,
         end_motif = start_motif + 21,
         start_motif_from_TTS = -47 + start_motif,
         end_motif_from_TTS = -47 + end_motif,
         gene = str_split_fixed(str_split_fixed(V1, ">", 2)[,2], "_site_1", 2)[,1],
         group = "motif1") 

position_pfu <- read.table(here("data/meme_data/pfu_motif_terminator.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(grepl("offset", V1)) %>%
  mutate(start_motif = as.numeric(str_split_fixed(V1, "offset=",2)[,2]) + 1,
         end_motif = start_motif + 21,
         start_motif_from_TTS = -47 + start_motif,
         end_motif_from_TTS = -47 + end_motif,
         gene = str_split_fixed(str_split_fixed(V1, ">", 2)[,2], "_site_1", 2)[,1],
         group = "motif1") 

position_hvo <- read.table(here("data/meme_data/hvo_motif_terminator.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(grepl("offset", V1)) %>%
  mutate(start_motif = as.numeric(str_split_fixed(V1, "offset=",2)[,2]) + 1,
         end_motif = start_motif + 6,
         start_motif_from_TTS = -46 + start_motif,
         end_motif_from_TTS = start_motif_from_TTS + 6,
         gene = str_split_fixed(str_split_fixed(V1, ">", 2)[,2], "_site_1", 2)[,1],
         group = "motif1") 


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MOTIF PLOTTING
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...............................plot best motifs
gg_motif_term_ecoli1 <- ggarrange(plot_rna_motif(motif_ecoli1), plot_rna_position(position_ecoli1), nrow = 2, heights = c(0.7,0.3))
gg_motif_term_ecoli2 <- ggarrange(plot_rna_motif(motif_ecoli2), plot_rna_position(position_ecoli2), nrow = 2, heights = c(0.7,0.3))
gg_motif_term_pfu    <- ggarrange(plot_rna_motif(motif_pfu), plot_rna_position(position_pfu), nrow = 2, heights = c(0.7,0.3))
gg_motif_term_hvo    <- ggarrange(plot_rna_motif(motif_hvo), plot_rna_position(position_hvo), nrow = 2, heights = c(0.7,0.3))

#...............................save plots (Supplementary Fig. 7a,c,d)
pdf(here("figures/motif1_term_ecoli.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_motif_term_ecoli1
dev.off()

pdf(here("figures/motif2_term_ecoli.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_motif_term_ecoli2
dev.off()

pdf(here("figures/motif_term_pfu.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_motif_term_pfu
dev.off()

pdf(here("figures/motif_term_hvo.pdf"),
    width = 7, height = 7, paper = "special",onefile=FALSE)
gg_motif_term_hvo
dev.off()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# NUCLEOTIDE ENRICHMENT PLOT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...............................filter for two motifs in ecoli
motif_ecoli1_genelist <- read.table(here("data/meme_data/ecoli_motif1_terminator.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(grepl("_site", V1)) %>%
  mutate(gene = str_split_fixed(str_split_fixed(V1, "_site_1",2)[,1], ">",2)[,2]) %>%
  dplyr::select(gene)

ecoli_utr3_motif1 <- ecoli_utr3 %>%
  dplyr::filter(gene %in% motif_ecoli1_genelist$gene)

motif_ecoli2_genelist <- read.table(here("data/meme_data/ecoli_motif2_terminator.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(grepl("_site", V1)) %>%
  mutate(gene = str_split_fixed(str_split_fixed(V1, "_site_1",2)[,1], ">",2)[,2]) %>%
  dplyr::select(gene)

ecoli_utr3_motif2 <- ecoli_utr3 %>%
  dplyr::filter(gene %in% motif_ecoli2_genelist$gene)

#..............................terminator sequence
ecoli_utr3$terminator_sequence <- gsub(ecoli_utr3$terminator_sequence, pattern = "T", replacement = "U")
ecoli_utr3_motif1$terminator_sequence <- gsub(ecoli_utr3_motif1$terminator_sequence, pattern = "T", replacement = "U")
ecoli_utr3_motif2$terminator_sequence <- gsub(ecoli_utr3_motif2$terminator_sequence, pattern = "T", replacement = "U")
pfu_utr3$terminator_sequence <- gsub(pfu_utr3$terminator_sequence, pattern = "T", replacement = "U")
hvo_utr3$terminator_sequence <- gsub(hvo_utr3$terminator_sequence, pattern = "T", replacement = "U")

#..............................background model
ecoli_fasta_bg <- as.character(unlist(readDNAStringSet(here("data/genome_data/ecoli_intergenic.fasta"))))
pfu_fasta_bg <- as.character(unlist(readDNAStringSet(here("data/genome_data/pfu_intergenic.fasta"))))
hvo_fasta_bg <- as.character(unlist(readDNAStringSet(here("data/genome_data/hvo_intergenic.fasta"))))

#..............................set seed
set.seed(87)

#..............................make subsamples (run 10000 times)
substrings_ecoli <- replicate(10000, sampleString(ecoli_fasta_bg))
substrings_ecoli <- gsub(substrings_ecoli, pattern = "T", replacement = "U")
substrings_pfu <- replicate(10000, sampleString(pfu_fasta_bg))
substrings_pfu <- gsub(substrings_pfu, pattern = "T", replacement = "U")
substrings_hvo <- replicate(10000, sampleString(hvo_fasta_bg))
substrings_hvo <- gsub(substrings_hvo, pattern = "T", replacement = "U")

#..............................make matrix for data and bg model and join tables
ecoli_terminator_matrix1 <- make_matrix_c(ecoli_utr3_motif1$terminator_sequence)
ecoli_background_matrix <- make_matrix_c(substrings_ecoli)
ecoli_compare_matrix1    <- ecoli_terminator_matrix1 %>%
  dplyr::rename(value_terminator = 2) %>%
  left_join(ecoli_background_matrix) %>%
  mutate(log_value = log2(value_terminator/value)) 

ecoli_terminator_matrix2 <- make_matrix_c(ecoli_utr3_motif2$terminator_sequence)
ecoli_background_matrix <- make_matrix_c(substrings_ecoli)
ecoli_compare_matrix2    <- ecoli_terminator_matrix2 %>%
  dplyr::rename(value_terminator = 2) %>%
  left_join(ecoli_background_matrix) %>%
  mutate(log_value = log2(value_terminator/value)) 

ecoli_terminator_matrix <- make_matrix_c(ecoli_utr3$terminator_sequence)
ecoli_background_matrix <- make_matrix_c(substrings_ecoli)
ecoli_compare_matrix    <- ecoli_terminator_matrix %>%
  dplyr::rename(value_terminator = 2) %>%
  left_join(ecoli_background_matrix) %>%
  mutate(log_value = log2(value_terminator/value)) 

pfu_terminator_matrix <- make_matrix_c(pfu_utr3$terminator_sequence)
pfu_background_matrix <- make_matrix_c(substrings_pfu)
pfu_compare_matrix    <- pfu_terminator_matrix %>%
  dplyr::rename(value_terminator = 2) %>%
  left_join(pfu_background_matrix) %>%
  mutate(log_value = log2(value_terminator/value)) 

hvo_terminator_matrix <- make_matrix_c(hvo_utr3$terminator_sequence)
hvo_background_matrix <- make_matrix_c(substrings_hvo)
hvo_compare_matrix    <- hvo_terminator_matrix %>%
  dplyr::rename(value_terminator = 2) %>%
  left_join(hvo_background_matrix) %>%
  mutate(log_value = log2(value_terminator/value)) 

#..............................plot nucleotide enrichment plots
nucleotide_enrich_motif1_ecoli <- nucleotide_enrichment_plotting(ecoli_compare_matrix1)
nucleotide_enrich_motif2_ecoli <- nucleotide_enrichment_plotting(ecoli_compare_matrix2)
nucleotide_enrich_all_ecoli    <- nucleotide_enrichment_plotting(ecoli_compare_matrix)
nucleotide_enrich_all_pfu      <- nucleotide_enrichment_plotting(pfu_compare_matrix)
nucleotide_enrich_all_hvo      <- nucleotide_enrichment_plotting(hvo_compare_matrix)

#..............................save nucleotide enrichment plots (Supplementary Fig. 7b)
pdf(here("figures/nucleotide_enrichment_motif1_ecoli.pdf"),
    width = 14, height = 7, paper = "special",onefile=FALSE)
nucleotide_enrich_motif1_ecoli
dev.off()

pdf(here("figures/nucleotide_enrichment_motif2_ecoli.pdf"),
    width = 14, height = 7, paper = "special",onefile=FALSE)
nucleotide_enrich_motif2_ecoli
dev.off()

#..............................save nucleotide enrichment plots (Fig. 2e)
pdf(here("figures/nucleotide_enrichment_all_ecoli.pdf"),
    width = 14, height = 7, paper = "special",onefile=FALSE)
nucleotide_enrich_all_ecoli
dev.off()

pdf(here("figures/nucleotide_enrichment_all_pfu.pdf"),
    width = 14, height = 7, paper = "special",onefile=FALSE)
nucleotide_enrich_all_pfu
dev.off()

pdf(here("figures/nucleotide_enrichment_all_hvo.pdf"),
    width = 14, height = 7, paper = "special",onefile=FALSE)
nucleotide_enrich_all_hvo
dev.off()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# WRITE TTS OUTPUT TO TABLE
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...............................ecoli
#...................load annotation data
ecoli_annotation <- fread(here("data/genome_data/ecoli_annotation.tsv"))

#...................write table
ecoli_export <- ecoli_utr3 %>%
  left_join(ecoli_annotation, by = c("gene" = "id")) %>%
  dplyr::rename(GeneID = gene, TTS = median_utr3, annotation = name) %>%
  dplyr::select(GeneID, TTS, strand, utr3_length, locus_tag, annotation) 

writexl::write_xlsx(x = ecoli_export, path = here("tables/tts_tables/tts_ecoli.xlsx"))

#...............................pfu
#...................load annotation data
pfu_annotation <- fread(here("data/genome_data/pfu_annotation.tsv"))
pfu_utr3$gene
#...................write table
pfu_export <- pfu_utr3 %>%
  rowwise() %>%
  mutate(gene = str_split_fixed(gene, ".p01", 2)[1]) %>%
  left_join(pfu_annotation, by = c("gene" = "gene")) %>%
  dplyr::rename(GeneID = gene, TTS = median_utr3, annotation = name) %>%
  dplyr::select(GeneID, TTS, strand, utr3_length, old_name, annotation) 

writexl::write_xlsx(x = pfu_export, path = here("tables/tts_tables/tts_pfu.xlsx"))

#...............................hvo
#...................load annotation data
hvo_annotation <- fread(here("data/genome_data/hvo_annotation.tsv"))

#...................write table
hvo_export <- hvo_utr3 %>%
  left_join(hvo_annotation, by = c("gene" = "id_name")) %>%
  dplyr::rename(GeneID = gene, TTS = median_utr3, annotation = name) %>%
  dplyr::select(GeneID, TTS, strand, utr3_length, old_name, annotation) 

writexl::write_xlsx(x = pfu_export, path = here("tables/tts_tables/tts_hvo.xlsx"))


