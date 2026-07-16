### wrangle data ###

library(readr)
library(dplyr)
library(ape)
library(phytools)
library(brms)
library(rstan)
library(loo)
library(ggplot2)
library(patchwork)
library(tidyr)
library(forcats)
library(ggExtra)

# fasta headers and contig names sometimes differ
# the sense checks at the end ensure everything merged correctly

data_dir = "~/Desktop/2026_plasmid_gc/enterobacterales_analysis/data"
setwd("~")
setwd(data_dir)

sample_sheet_file <- file.path("sample.csv")

sample_sheet <- read_csv(sample_sheet_file, show_col_types = FALSE) %>%
  dplyr::rename("genome_id" = genome, "contig_id" = contig, "contig_type" = type) %>%
  mutate(genome_id = gsub("_contig", "", genome_id)) %>%
  mutate(contig_id = gsub("_contig", "", contig_id)) %>%
  mutate(contig_id = ifelse(contig_type =="chromosome", gsub("_1", "_chromo", contig_id), contig_id))

# all chromosomal contig_id's end in _chromo
# plasmid chromosomal contigs can end in _1, _2, ...

results_dir = "~/Desktop/2026_plasmid_gc/enterobacterales_analysis/results"
setwd("~")
setwd(results_dir)

contig_stats_file <- file.path("contig_statistics.tsv")
pling_file <- file.path("pling_results/dcj_thresh_4_graph/objects/typing.tsv")
pling_dist_file <- file.path("pling_results/all_plasmids_distances.tsv")
pling_hubs_file <- file.path("pling_results/dcj_thresh_4_graph/objects/hub_plasmids.csv")
plasmidfinder_file <- file.path("plasmidfinder.tsv")
mobtyper_file <- file.path("mobtyper.tsv")
mlst_file <- file.path("mlst.tsv")

contig_stats <- read_tsv(contig_stats_file, show_col_types = FALSE) %>%
  mutate(contig_id_clean = gsub("_contig", "", contig_id)) %>%
  mutate(contig_id_clean = ifelse(length>3480000, gsub("_1", "_chromo", contig_id_clean), contig_id_clean)) %>%
  mutate(contig_id_clean = case_when(
    contig_id_clean == "27e2513d-6427-411b-ba7c-0f2fce69a57b_NODE_edge_2+_length_11432_cov_84" ~ 
      "27e2513d-6427-411b-ba7c-0f2fce69a57b_2",
    contig_id_clean == "27e2513d-6427-411b-ba7c-0f2fce69a57b_NODE_edge_chromo+_length_5092231_cov_447" ~
      "27e2513d-6427-411b-ba7c-0f2fce69a57b_chromo",
    contig_id_clean == "1" ~
      "73d0f482-64ad-42de-a177-9f8ab30ea9e1_chromo",
    contig_id_clean == "2" ~
      "73d0f482-64ad-42de-a177-9f8ab30ea9e1_2",
    contig_id_clean == "3" ~
      "73d0f482-64ad-42de-a177-9f8ab30ea9e1_3",
    contig_id_clean == "2233c4fd-c922-4fa8-8b2a-e47069550774_2" ~
      "2233c4fd-c922-4fa8-8b2a-e47069550774_chromo",
    contig_id_clean == "22a9cac9-e769-4e6b-bed8-01f835c0e170_2" ~
      "22a9cac9-e769-4e6b-bed8-01f835c0e170_chromo",
    contig_id_clean == "22a9cac9-e769-4e6b-bed8-01f835c0e170_2" ~
      "22a9cac9-e769-4e6b-bed8-01f835c0e170_chromo",
    contig_id_clean == "719706c9-d81a-499b-bec6-25781f8ff82d_2" ~
      "2719706c9-d81a-499b-bec6-25781f8ff82d_chromo",
    TRUE ~ contig_id_clean))

pling <- read_tsv(pling_file, show_col_types = FALSE) %>%
  dplyr::rename("contig_id" = plasmid, "subcommunity" = type) %>%
  group_by(subcommunity) %>%
  mutate(subcommunity_size = n())

hub_plasmids <- read_csv(pling_hubs_file) %>% pull(hub_plasmids)

plasmidfinder <- read_tsv(plasmidfinder_file, show_col_types = FALSE) %>%
  select(SEQUENCE, GENE) %>%
  group_by(SEQUENCE) %>%
  summarise(
    replicon_types = paste(sort(GENE), collapse = ", "),
    .groups = "drop"
  ) %>%
  dplyr::rename("contig_id" = SEQUENCE)

mobtyper <- read_tsv(mobtyper_file, show_col_types = FALSE) %>%
  select(sample_id, predicted_mobility) %>%
  mutate(contig_id = sub(" .*", "", sample_id)) %>%
  select(contig_id, predicted_mobility)

mlst <- read_tsv(mlst_file, show_col_types = FALSE, col_names = FALSE) %>%
  mutate(contig_id = sub("\\.fasta$", "", basename(X1))) %>%
  select(contig_id, X2) %>%
  mutate(
    mlst = case_when(
      X2 == "ecoli_achtman_4"     ~ "escherichia_coli",
      X2 == "koxytoca"            ~ "klebsiella_oxytoca",
      X2 == "klebsiella"          ~ "klebsiella_pneumoniae",
      X2 == "ecloacae"            ~ "enterobacter_cloacae",
      X2 == "proteus"             ~ "proteus_spp",
      X2 == "cfreundii"           ~ "citrobacter_freundii",
      X2 == "kaerogenes"          ~ "klebsiella_aerogenes",
      X2 == "aeromonas"           ~ "aeromonas_spp",
      X2 == "senterica_achtman_2" ~ "salmonella_enterica",
      X2 == "-"                   ~ "untyped"
    )
  ) %>%
  select(-X2) %>%
  mutate(contig_id = gsub("_contig", "", contig_id)) %>%
  mutate(contig_id = gsub("_1", "_chromo", contig_id)) 

combined <- merge(sample_sheet, contig_stats, by.x = "contig_id", by.y = "contig_id_clean", all.x = TRUE) %>%
  merge(., pling, by.x = "contig_id.y", by.y = "contig_id", all.x = TRUE) %>%
  mutate(subcommunity = ifelse(contig_id %in% hub_plasmids, "hub", subcommunity)) %>%
  merge(., plasmidfinder, by.x = "contig_id.y", by.y = "contig_id", all.x = TRUE) %>%
  merge(., mobtyper, by.x = "contig_id.y", by.y = "contig_id", all.x = TRUE) %>%
  merge(., mlst, by.x = "contig_id", by.y = "contig_id", all.x = TRUE) %>%
  select(-contig_id.y) 

combined[combined$contig_id == "73d0f482-64ad-42de-a177-9f8ab30ea9e1_2", "subcommunity"] <- "community_0_subcommunity_695"
combined[combined$contig_id == "73d0f482-64ad-42de-a177-9f8ab30ea9e1_2", "subcommunity_size"] <- 30

combined[combined$contig_id == "27e2513d-6427-411b-ba7c-0f2fce69a57b_2", "subcommunity"] <- "community_13_subcommunity_1"
combined[combined$contig_id == "27e2513d-6427-411b-ba7c-0f2fce69a57b_2", "subcommunity_size"] <- 23

combined[combined$contig_id == "73d0f482-64ad-42de-a177-9f8ab30ea9e1_3", "subcommunity"] <- "community_0_subcommunity_709"
combined[combined$contig_id == "73d0f482-64ad-42de-a177-9f8ab30ea9e1_3", "subcommunity_size"] <- 110

combined[combined$contig_id == "719706c9-d81a-499b-bec6-25781f8ff82d_chromo", "length"] <- 5357707
combined[combined$contig_id == "719706c9-d81a-499b-bec6-25781f8ff82d_chromo", "gc_percent"] <- 57.34

# sense checks

sum(combined$contig_type == "plasmid", na.rm=TRUE) # 3682
sum(combined$contig_type == "chromosome", na.rm=TRUE) # 1454

sum(combined$contig_type == "plasmid" & !is.na(combined$length), na.rm=TRUE) # 3682
sum(combined$contig_type == "chromosome" & !is.na(combined$length), na.rm=TRUE) # 1454

sum(combined$contig_type == "plasmid" & !is.na(combined$subcommunity), na.rm=TRUE) # 3682
sum(combined$contig_type == "plasmid" & combined$subcommunity!="hub", na.rm=TRUE) # 3634
sum(combined$contig_type == "plasmid" & combined$subcommunity=="hub", na.rm=TRUE) # 48

sum(combined$contig_type == "plasmid" & !is.na(combined$replicon_types), na.rm=TRUE) # 3081

sum(combined$contig_type == "plasmid" & !is.na(combined$predicted_mobility), na.rm=TRUE) # 3682

sum(combined$contig_type == "chromosome" & !is.na(combined$mlst), na.rm=TRUE) # 1454
sum(combined$contig_type == "chromosome" & combined$mlst=="untyped", na.rm=TRUE) # 4

### make table 1 ###

combined %>%
  filter(contig_type=="chromosome") %>%
  group_by(mlst) %>%
  summarise(
    n = n(),
    mean_gc = mean(gc_percent),
    sd_gc = sd(gc_percent)
  ) %>%
  arrange(-n) 

### fit main model ###

# prepare data
chrom_gc <- combined %>% 
  filter(contig_type == "chromosome") %>% 
  select(genome_id, chrom_gc = gc_percent)

plasmids <- combined %>% 
  filter(contig_type == "plasmid") %>% 
  left_join(chrom_gc, by = "genome_id") %>%
  filter(subcommunity!="hub") %>%
  group_by(subcommunity) %>% 
  filter(n() >= 5) %>% 
  ungroup()

# read & clean tree
tree <- read.tree("chromosome_mashtree.dnd")
tree$tip.label <- sub("_.*$", "", tree$tip.label)

# keep only genomes present in plasmids and prune tree
plasmids <- plasmids %>% filter(genome_id %in% tree$tip.label)
tree <- drop.tip(tree, setdiff(tree$tip.label, plasmids$genome_id))

# midpoint root and make ultrametric
tree_mid   <- midpoint.root(tree)
tree_ultra <- chronos(tree_mid, lambda = 1, model = "correlated")

# phylogenetic covariance matrix
A <- vcv(tree_ultra, corr = FALSE)

# ensure A matches genome order
genomes <- sort(unique(plasmids$genome_id))
if (!all(genomes %in% rownames(A))) stop("Some genomes present in data are missing from A")
A <- A[genomes, genomes]
plasmids$genome_id <- factor(plasmids$genome_id, levels = genomes)

# scale/centre continuous predictors
plasmids <- plasmids %>%
  mutate(
    log10_length = log10(length),
    log10_length_c = log10_length - mean(log10_length, na.rm = TRUE),
    chrom_gc_c = chrom_gc - mean(chrom_gc, na.rm = TRUE),
    gc_diff = gc_percent - chrom_gc,
    subcommunity = factor(subcommunity),
    predicted_mobility = factor(predicted_mobility),
    community = sub("(_subcommunity_.*)$", "", subcommunity)
  ) %>%
  select(gc_percent, gc_diff, chrom_gc, chrom_gc_c, length, log10_length_c, predicted_mobility, community, subcommunity, genome_id) 

# fit model with smooths
# fit_smooth <- brm(
#   gc_diff ~ s(log10_length_c) + predicted_mobility + s(chrom_gc_c) +
#     (1 + chrom_gc_c | community/subcommunity) + (1 | gr(genome_id, cov = A)),
#   data = plasmids,
#   data2 = list(A = A),
#   family = gaussian(),
#   chains = 4,
#   iter = 4000,
#   control = list(adapt_delta = 0.95)
# )

# saveRDS(fit_smooth, "fit_smooth.rds")

fit_smooth <- readRDS("fit_smooth.rds")

# > summary(fit_smooth)
# Family: gaussian 
# Links: mu = identity 
# Formula: gc_diff ~ s(log10_length_c) + predicted_mobility + s(chrom_gc_c) + (1 + chrom_gc_c | community/subcommunity) + (1 | gr(genome_id, cov = A)) 
# Data: plasmids (Number of observations: 2305) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Smoothing Spline Hyperparameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sds(slog10_length_c_1)    17.46      5.46     9.08    30.29 1.00     1966     2518
# sds(schrom_gc_c_1)         1.92      1.84     0.05     6.59 1.00     1524     3452
# 
# Multilevel Hyperparameters:
#   ~community (Number of levels: 21) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)                 5.35      1.06     3.57     7.65 1.00     2725     3980
# sd(chrom_gc_c)                0.29      0.18     0.02     0.71 1.00     1342     1684
# cor(Intercept,chrom_gc_c)     0.02      0.39    -0.74     0.80 1.00     3888     3160
# 
# ~community:subcommunity (Number of levels: 89) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)                 2.27      0.22     1.89     2.73 1.00     1425     2667
# sd(chrom_gc_c)                0.14      0.06     0.04     0.26 1.00     1239     1508
# cor(Intercept,chrom_gc_c)     0.13      0.30    -0.43     0.73 1.00     3367     2560
# 
# ~genome_id (Number of levels: 1010) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.56      0.12     0.31     0.78 1.00     1219     1766
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                            -5.40      1.26    -7.80    -2.87 1.00     2046     3282
# predicted_mobilitymobilizable         0.24      0.30    -0.34     0.83 1.00     4234     5622
# predicted_mobilitynonMmobilizable    -1.46      0.27    -1.99    -0.94 1.00     4167     5543
# slog10_length_c_1                   -66.32     11.14   -88.10   -45.25 1.00     3029     4534
# schrom_gc_c_1                       -10.94      2.87   -16.07    -4.16 1.00     3653     3482
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     1.55      0.03     1.50     1.60 1.00     3857     5114
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# both smooth terms look consistent with linearity - drop them
p_cond_length <- plot(
  conditional_effects(fit_smooth, "log10_length_c", re_formula = NA),
  ask = FALSE
)[[1]] +
  labs(
    x = expression(log[10]*"-scaled plasmid length (bp; centred)"),
    y = expression(Delta*"GC")
  ) 

p_cond_gc <- plot(
  conditional_effects(fit_smooth, "chrom_gc_c", re_formula = NA),
  ask = FALSE
)[[1]] +
  labs(
    x = "Chromosome GC content (centred)",
    y = expression(Delta*"GC")
  ) 

p_cond_length / p_cond_gc +
  plot_annotation(tag_levels = "a")

# fit_linear <- brm(
#   gc_diff ~ log10_length_c + predicted_mobility + chrom_gc_c +
#     (1 + chrom_gc_c | community/subcommunity) + (1 | gr(genome_id, cov = A)),
#   data = plasmids,
#   data2 = list(A = A),
#   family = gaussian(),
#   chains = 4,
#   iter = 4000,
#   control = list(adapt_delta = 0.95)
# )

# saveRDS(fit_linear, "fit_linear.rds")

fit_linear <- readRDS("fit_linear.rds")

# > summary(fit_linear)
# Family: gaussian 
# Links: mu = identity 
# Formula: gc_diff ~ log10_length_c + predicted_mobility + chrom_gc_c + (1 + chrom_gc_c | community/subcommunity) + (1 | gr(genome_id, cov = A)) 
# Data: plasmids (Number of observations: 2305) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
#   ~community (Number of levels: 21) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)                 5.87      1.15     4.00     8.43 1.00     2720     4163
# sd(chrom_gc_c)                0.16      0.12     0.01     0.46 1.00     1916     2600
# cor(Intercept,chrom_gc_c)    -0.11      0.48    -0.92     0.87 1.00     5366     3813
# 
# ~community:subcommunity (Number of levels: 89) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)                 2.50      0.24     2.09     3.02 1.00      974     2204
# sd(chrom_gc_c)                0.14      0.05     0.06     0.24 1.00     2269     3643
# cor(Intercept,chrom_gc_c)     0.18      0.30    -0.40     0.74 1.00     3196     4038
# 
# ~genome_id (Number of levels: 1010) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.53      0.12     0.28     0.76 1.00     1360     1984
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                            -4.10      1.37    -6.79    -1.36 1.00     1818     2852
# log10_length_c                        1.03      0.19     0.64     1.40 1.00     3595     5192
# predicted_mobilitymobilizable        -1.01      0.28    -1.55    -0.47 1.00     3241     4718
# predicted_mobilitynonMmobilizable    -2.18      0.26    -2.69    -1.66 1.00     3487     5119
# chrom_gc_c                           -1.03      0.11    -1.27    -0.84 1.00     3294     3727
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     1.62      0.03     1.57     1.68 1.00     4254     4907
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

loose_groups <- plasmids %>%
  group_by(subcommunity) %>%
  summarise(
    sd_gc = sd(gc_percent)
  ) %>%
  filter(sd_gc >= 1) %>%
  pull(subcommunity)

plasmids_loose <- plasmids %>%
  filter(subcommunity %in% loose_groups)

# 27/89 subcommunities representing 1530/2073 plasmids

# fit_loose <- brm(
#   gc_diff ~ log10_length_c + predicted_mobility + chrom_gc_c +
#     (1 + chrom_gc_c | community/subcommunity) + (1 | gr(genome_id, cov = A)),
#   data = plasmids_loose,
#   data2 = list(A = A),
#   family = gaussian(),
#   chains = 4,
#   iter = 4000,
#   control = list(adapt_delta = 0.95)
# )

# saveRDS(fit_loose, "fit_loose.rds")

fit_loose <- readRDS("fit_loose.rds")

# > summary(fit_loose)
# Family: gaussian 
# Links: mu = identity 
# Formula: gc_diff ~ log10_length_c + predicted_mobility + chrom_gc_c + (1 + chrom_gc_c | community/subcommunity) + (1 | gr(genome_id, cov = A)) 
# Data: plasmids_loose (Number of observations: 1530) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
#   ~community (Number of levels: 6) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)                 3.98      2.14     0.39     8.69 1.00     2123     2383
# sd(chrom_gc_c)                0.49      0.45     0.03     1.56 1.00     1915     2865
# cor(Intercept,chrom_gc_c)    -0.15      0.51    -0.95     0.85 1.00     4372     4044
# 
# ~community:subcommunity (Number of levels: 27) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)                 2.93      0.55     2.06     4.19 1.00     1331     2906
# sd(chrom_gc_c)                0.14      0.07     0.03     0.29 1.00     2738     2833
# cor(Intercept,chrom_gc_c)     0.15      0.41    -0.60     0.89 1.00     4045     3526
# 
# ~genome_id (Number of levels: 760) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.51      0.19     0.11     0.86 1.01      687      863
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                            -1.58      1.81    -5.12     2.05 1.00     3873     4397
# log10_length_c                        0.81      0.26     0.30     1.31 1.00     5816     6255
# predicted_mobilitymobilizable        -1.37      0.37    -2.10    -0.64 1.00     5271     5707
# predicted_mobilitynonMmobilizable    -2.49      0.35    -3.16    -1.80 1.00     5956     5695
# chrom_gc_c                           -1.08      0.31    -1.76    -0.56 1.00     3075     3192
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     1.96      0.04     1.89     2.05 1.00     3006     5008
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

### ICCs ###

#draws <- brms::as_draws_df(fit_linear)
draws <- brms::as_draws_df(fit_loose)

var_community    <- (draws$sd_community__Intercept)^2
var_subcommunity <- (draws$`sd_community:subcommunity__Intercept`)^2
var_phylo        <- (draws$sd_genome_id__Intercept)^2 
var_residual     <- (draws$sigma)^2

var_total <- var_community + var_subcommunity + var_phylo + var_residual

icc_plasmid_structure <- (var_community + var_subcommunity) / var_total
icc_phylo             <- var_phylo / var_total

quantile(icc_plasmid_structure, probs = c(0.025, 0.5, 0.975)) * 100
quantile(icc_phylo, probs = c(0.025, 0.5, 0.975)) * 100

### figure ###

# prepare plotting data
chrom_gc <- combined %>%
  filter(contig_type == "chromosome") %>%
  select(genome_id, chrom_gc = gc_percent, mlst)

plot_df <- combined %>%
  filter(contig_type == "plasmid") %>%
  left_join(chrom_gc, by = "genome_id") %>%
  filter(subcommunity != "hub") %>%
  group_by(subcommunity) %>%
  filter(n() >= 5) %>%
  ungroup()

# species ordering by plasmid count
species_summary <- plot_df %>%
  count(mlst.y, name = "plasmids") %>%
  left_join(
    plot_df %>%
      group_by(mlst.y) %>%
      summarise(
        subcommunities = n_distinct(subcommunity),
        .groups = "drop"
      ),
    by = "mlst.y"
  ) %>%
  arrange(desc(plasmids))

# single consistent ordering used everywhere
species_order <- species_summary$mlst.y

plot_df <- plot_df %>%
  mutate(
    mlst.y = factor(mlst.y, levels = species_order)
  )

plot_counts <- species_summary %>%
  mutate(
    mlst.y = factor(mlst.y, levels = species_order)
  ) %>%
  pivot_longer(
    c(plasmids, subcommunities),
    names_to = "metric",
    values_to = "count"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("plasmids", "subcommunities"),
      labels = c("Plasmids", "Pling subcommunities")
    )
  )

# ensure consistent where applicable between fig 1 and fig 2
species_colours <- c(
  escherichia_coli      = "#E69F00",
  klebsiella_pneumoniae = "#56B4E9",
  citrobacter_freundii  = "#882255",
  enterobacter_cloacae  = "#F0E442",
  klebsiella_oxytoca    = "#0072B2",
  proteus_spp           = "#D55E00",
  aeromonas_spp         = "#CC79A7",
  klebsiella_aerogenes  = "#332288",
  salmonella_enterica   = "#009E73",
  untyped               = "grey70"
)

species_labels <- c(
  escherichia_coli      = "italic('Escherichia coli')",
  klebsiella_pneumoniae = "italic('Klebsiella pneumoniae')",
  citrobacter_freundii  = "italic('Citrobacter freundii')",
  enterobacter_cloacae  = "italic('Enterobacter cloacae')",
  klebsiella_oxytoca    = "italic('Klebsiella oxytoca')",
  proteus_spp           = "italic('Proteus')~spp.",
  aeromonas_spp         = "italic('Aeromonas')~spp.",
  klebsiella_aerogenes  = "italic('Klebsiella aerogenes')",
  salmonella_enterica   = "italic('Salmonella enterica')",
  untyped               = "'Untyped'"
)

common_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(face = "bold", colour = "white")
  )

p_counts <-
  ggplot(plot_counts,
         aes(fct_rev(mlst.y), count, fill = mlst.y)) +
  geom_col() +
  geom_point(aes(colour = mlst.y), alpha = 0) + 
  coord_flip() +
  facet_wrap(~metric, scales = "free_x") +
  scale_x_discrete(
    labels = function(x) parse(text = species_labels[x])
  ) +
  scale_fill_manual(
    values = species_colours,
    breaks = species_order,
    labels = function(x) parse(text = species_labels[x])
  ) +
  scale_colour_manual(
    values = species_colours,
    breaks = species_order,
    labels = function(x) parse(text = species_labels[x])
  ) +
  labs(
    x = NULL,
    y = "Count",
    colour = "Species"
  ) +
  guides(
    fill = "none", 
    colour = guide_legend(override.aes = list(alpha = 0.6, size = 1))
  ) +
  common_theme 

p_gc <-
  ggplot(plot_df,
         aes(chrom_gc, gc_percent, colour = mlst.y)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey20"
  ) +
  scale_colour_manual(
    values = species_colours,
    breaks = species_order,
    labels = function(x) parse(text = species_labels[x])
  ) +
  scale_x_continuous(
    breaks = seq(40, 55, by = 5),
    limits = c(37, 60),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(35, 55, by = 5),
    limits = c(30, 59),
    expand = c(0, 0)
  ) +
  labs(
    x = "Host chromosome GC content",
    y = "Plasmid GC content"
  ) +
  common_theme

p_gc_no_legend <- p_gc + theme(legend.position = "none")

p_gc_marginal <- ggMarginal(
  p_gc_no_legend,
  type = "density",
  margins = "both",
  groupColour = FALSE,
  size = 5,
  fill = "grey85",
  colour = "grey40"
)

tag_style <- element_text(face = "bold", size = 20)

p_combined <- (
  free(p_counts) / 
    (wrap_elements(p_gc_marginal) + theme(plot.tag = tag_style))
) +
  plot_layout(
    heights = c(1.5, 2),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = 'a') &
  theme(
    legend.position = "right",
    plot.tag = tag_style
  )

p_combined
