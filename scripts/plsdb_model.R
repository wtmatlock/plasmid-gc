### wrangle data ###

library(tidyverse)
library(Biostrings)
library(janitor)

nuccore <- read_csv("plsdb/nuccore.csv", show_col_types = FALSE) %>%
  mutate(NUCCORE_UID = as.character(NUCCORE_UID))

assembly <- read_csv("plsdb/assembly.csv", show_col_types = FALSE) %>%
  mutate(
    NUCCORE_UID_split = str_remove_all(NUCCORE_UID, "\\[|\\]")
  ) %>%
  separate_rows(NUCCORE_UID_split, sep = ",\\s*") %>%
  mutate(NUCCORE_UID_split = as.character(NUCCORE_UID_split))

biosample <- read_csv("plsdb/biosample.csv", show_col_types = FALSE) %>%
  mutate(BIOSAMPLE_taxonomy = as.character(BIOSAMPLE_taxonomy))

taxonomy <- read_csv("plsdb/taxonomy.csv", show_col_types = FALSE) %>%
  mutate(TAXONOMY_UID = as.character(TAXONOMY_UID))

plsdb_df <- nuccore %>%
  left_join(assembly, by = c("NUCCORE_UID" = "NUCCORE_UID_split")) %>%
  left_join(biosample, by = "BIOSAMPLE_UID") %>%
  left_join(taxonomy, by = c("BIOSAMPLE_taxonomy" = "TAXONOMY_UID"))

assemblies <- assembly %>%
  distinct(ASSEMBLY_ACC) %>%
  filter(!is.na(ASSEMBLY_ACC))

#writeLines(assemblies$ASSEMBLY_ACC, "accession_list.txt")

### NCBI datasets download ###

# alternatively, run this in the terminal
system("
source ~/miniconda3/etc/profile.d/conda.sh &&
conda activate datasets_env &&
datasets download genome accession \
  --inputfile accession_list.txt \
  --include genome \
  --dehydrated \
  --filename genomes.zip \
  --api-key your_api_key &&
unzip genomes.zip &&
datasets rehydrate --directory .
")

### organise downloaded data ###

seq_report_files <- list.files(
  path = "ncbi_dataset/data",
  pattern = "genomic\\.fna$",
  full.names = TRUE,
  recursive = TRUE
)

get_assembly_id <- function(path) {
  basename(dirname(path))
}

file_df <- tibble(fasta = seq_report_files) %>%
  mutate(ASSEMBLY_ACC = map_chr(fasta, get_assembly_id))

complete_assemblies <- plsdb_df %>%
  filter(
    ASSEMBLY_Status == "Complete Genome",
    NUCCORE_Topology == "circular"
  ) %>%
  pull(ASSEMBLY_ACC) %>%
  unique()

files_to_process <- file_df %>%
  filter(ASSEMBLY_ACC %in% complete_assemblies)

### calculate chromosomal GC-content ###

calc_chrom_stats <- function(fasta_path) {
  dna_seqs <- readDNAStringSet(fasta_path)
  
  is_plasmid <- grepl("plasmid", names(dna_seqs), ignore.case = TRUE)
  putative_chroms <- dna_seqs[!is_plasmid]
  
  if(length(putative_chroms) == 0) putative_chroms <- dna_seqs
  
  main_chrom <- putative_chroms[which.max(width(putative_chroms))]
  
  freqs <- letterFrequency(main_chrom, c("G","C","A","T"), OR = 0)
  gc <- sum(freqs[, c("G","C")])
  total_gcat <- sum(freqs)
  
  chrom_len <- unname(width(main_chrom))
  
  return(tibble(
    chromosomal_gc = 100 * gc / total_gcat,
    chromosome_length = chrom_len
  ))
}

gc_df <- files_to_process %>%
  mutate(
    stats = map(fasta, calc_chrom_stats)
  ) %>%
  unnest(stats)

### combine ###

plot_df <- plsdb_df %>%
  filter(
    ASSEMBLY_Status == "Complete Genome",
    NUCCORE_Topology == "circular"
  ) %>%
  inner_join(gc_df, by = "ASSEMBLY_ACC") %>%
  transmute(
    nuccore_uid = NUCCORE_UID.x,
    plasmid_gc = NUCCORE_GC * 100,
    chromosome_gc = chromosomal_gc,
    chromosome_length = chromosome_length,
    species = TAXONOMY_species
  ) %>%
  drop_na()

cor.test(plot_df$chromosome_gc, plot_df$plasmid_gc, method = "pearson")

### figure ###

top10 <- plot_df %>%
  count(species, sort = TRUE) %>%
  head(10) %>%
  pull(species)

plot_df <- plot_df %>%
  mutate(
    species_clean = str_replace_all(species, "_", " "),
    species_group = ifelse(species %in% top10, species_clean, "Other")
  ) %>%
  mutate(
    species_group = factor(
      species_group, 
      levels = c(str_replace_all(top10, "_", " "), "Other")
    )
  )

species_colours <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442", 
  "#0072B2",
  "#D55E00",
  "#CC79A7",
  "#332288",
  "#117733",
  "#882255",
  "grey70" # other
)

names(species_colours) <- levels(plot_df$species_group)

format_sp_labels <- function(breaks) {
  lapply(breaks, function(x) {
    if (is.na(x)) return(x)
    if (x == "Other") "Other" else bquote(italic(.(x)))
  })
}

plot_df

ggplot(plot_df, aes(x = chromosome_gc, y = plasmid_gc, colour = species_group)) +
  geom_point(data = filter(plot_df, species_group == "Other"), alpha = 0.3, size = 0.7) +
  geom_point(data = filter(plot_df, species_group != "Other"), alpha = 0.3, size = 0.7) +
  geom_abline(slope = 1, intercept = 0, colour = "grey10", linetype="dashed") +
  scale_x_continuous(limits=c(0,100), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
  scale_colour_manual(
    name = "Species",
    values = species_colours,
    labels = format_sp_labels
  ) +
  labs(
    x = "Host chromosome GC-content",
    y = "Plasmid GC-content"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white", face = "bold"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )
