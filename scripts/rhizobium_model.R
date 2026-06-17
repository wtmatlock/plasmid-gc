### wrangle data ###

library(readr)
library(readxl)
library(dplyr)

setwd("~")
setwd("~/Desktop/rhizobium_analysis")

gene_gc <- read_delim("rhizobium_gc.tsv", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

rhizosphere <- read_excel("pnas.2009094117.sd02.xlsx", 
                                   sheet = 5) %>%
  mutate(niche = "rhizosphere")
colonisation <- read_excel("pnas.2009094117.sd02.xlsx", 
                          sheet = 6) %>%
  mutate(niche = "colonisation")
nodule <- read_excel("pnas.2009094117.sd02.xlsx", 
                           sheet = 7) %>%
  mutate(niche = "nodule")
bacteroid <- read_excel("pnas.2009094117.sd02.xlsx", 
                     sheet = 8) %>%
  mutate(niche = "bacteroid")

gene_essentiality <- rbind(rhizosphere, colonisation, nodule, bacteroid) %>%
  group_by(No, `Gene Id`) %>%
  summarise(count = sum(grepl("ES", `State Call`))) 

data <- merge(gene_gc, gene_essentiality, by.x="old_locus_tag", by.y="Gene Id") %>%
  filter(!grepl("NC_008380", replicon))

chromosome_gc <- 61.09

data <- data %>%
  mutate(deltaGC1 = gc1 - chromosome_gc,
         deltaGC2 = gc2 - chromosome_gc,
         deltaGC3 = gc3 - chromosome_gc,
         deltaGC4 = gc4 - chromosome_gc) %>%
  na.omit() %>%
  mutate(replicon = as.factor(replicon))

### plot ###

library(ggplot2)
library(tidyr)

data_plot <- data %>%
  pivot_longer(
    cols = starts_with("deltaGC"),
    names_to = "delta_type",
    values_to = "delta_value"
  ) %>%
  mutate(
    at_rich = ifelse(delta_value >= 0.5, TRUE, FALSE),
    delta_type = recode(
      delta_type,
      "deltaGC1" = "First codon position",
      "deltaGC2" = "Second codon position",
      "deltaGC3" = "Third codon position",
      "deltaGC4" = "Four-fold degenerate sites"
    ),
    # Explicitly set the factor levels here
    delta_type = factor(delta_type, levels = c(
      "First codon position",
      "Second codon position",
      "Third codon position",
      "Four-fold degenerate sites"
    ))
  )

plot <- ggplot(data_plot, aes(x = as.factor(count), y = delta_value, color = at_rich)) +
  geom_jitter(alpha = 0.3, width = 0.25) +
  geom_boxplot(alpha = 0, fill = NA, colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ delta_type) +
  scale_color_manual(
    values = c("TRUE" = "#6FB1E4", "FALSE" = "#C76626"),
    breaks = c("TRUE", "FALSE"),
    labels = c("FALSE" = "Chromosome", "TRUE" = "Plasmid")
  ) +
  labs(
    x = "Essential niche count",
    y = "",
    color = "GC-content\nhigher on:"
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

plot

### fit model ###

library(brms)

formula <- bf(mvbind(deltaGC1, deltaGC2, deltaGC3, deltaGC4) ~ count + 
  (1 | replicon)) +
  set_rescor(TRUE)

fit <- brm(
  formula = formula,
  data = data,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  control = list(adapt_delta = 0.95)
)

data_plot %>%
  select(replicon, old_locus_tag, length, gc_content, 
         gc1, gc2, gc3, gc4, delta_type, delta_value, count) %>%
  rename("locus_tag" = old_locus_tag, "plasmid" = replicon,
         "n_essential_niches" = count) %>%
  arrange(plasmid, locus_tag, delta_type) %>%
  write.csv(., "rhizobium_output.csv", row.names = FALSE)

summary(fit)

# > summary(fit)
# Family: MV(gaussian, gaussian, gaussian, gaussian) 
# Links: mu = identity
# mu = identity
# mu = identity
# mu = identity 
# Formula: deltaGC1 ~ count + (1 | replicon) 
# deltaGC2 ~ count + (1 | replicon) 
# deltaGC3 ~ count + (1 | replicon) 
# deltaGC4 ~ count + (1 | replicon) 
# Data: data (Number of observations: 2305) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
#   ~replicon (Number of levels: 6) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(deltaGC1_Intercept)     0.60      0.38     0.08     1.51 1.00     1983     2455
# sd(deltaGC2_Intercept)     0.15      0.15     0.00     0.53 1.00     3437     4235
# sd(deltaGC3_Intercept)     4.15      1.63     2.13     8.25 1.00     3122     3947
# sd(deltaGC4_Intercept)     5.43      2.08     2.83    10.67 1.00     2658     4866
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# deltaGC1_Intercept     1.88      0.31     1.17     2.43 1.00     3054     3291
# deltaGC2_Intercept   -15.99      0.13   -16.23   -15.73 1.00     6801     4025
# deltaGC3_Intercept    11.64      1.75     8.13    15.19 1.00     2063     2985
# deltaGC4_Intercept    15.48      2.24    10.93    20.02 1.00     2469     3702
# deltaGC1_count         0.12      0.15    -0.16     0.41 1.00     9087     6155
# deltaGC2_count        -0.13      0.13    -0.39     0.12 1.00    11953     6151
# deltaGC3_count        -0.44      0.21    -0.86    -0.02 1.00     6861     6077
# deltaGC4_count        -0.40      0.26    -0.92     0.12 1.00     6793     5987
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma_deltaGC1     4.60      0.07     4.48     4.74 1.00    12243     6306
# sigma_deltaGC2     4.09      0.06     3.97     4.21 1.00    10620     5593
# sigma_deltaGC3     6.77      0.10     6.58     6.97 1.00     5767     6148
# sigma_deltaGC4     8.32      0.12     8.09     8.57 1.00     5582     5943
# 
# Residual Correlations: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# rescor(deltaGC1,deltaGC2)     0.30      0.02     0.26     0.33 1.00    11079     5943
# rescor(deltaGC1,deltaGC3)     0.01      0.02    -0.03     0.05 1.00     6244     5781
# rescor(deltaGC2,deltaGC3)    -0.10      0.02    -0.15    -0.06 1.00     6526     5828
# rescor(deltaGC1,deltaGC4)     0.08      0.02     0.04     0.12 1.00     6310     5839
# rescor(deltaGC2,deltaGC4)    -0.12      0.02    -0.16    -0.08 1.00     6279     5298
# rescor(deltaGC3,deltaGC4)     0.93      0.00     0.92     0.93 1.00     6259     6003
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
