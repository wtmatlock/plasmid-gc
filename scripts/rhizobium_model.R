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
      "deltaGC1" = "Delta*GC[1]",
      "deltaGC2" = "Delta*GC[2]",
      "deltaGC3" = "Delta*GC[3]",
      "deltaGC4" = "Delta*GC[4]"
    )
  )

plot <- ggplot(data_plot, aes(x = as.factor(count), y = delta_value, color = at_rich)) +
  geom_jitter(alpha = 0.3, width = 0.25) +
  geom_boxplot(alpha = 0, fill = NA, colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ delta_type, labeller = label_parsed) +
  scale_color_manual(
    values = c("TRUE" = "#6FB1E4", "FALSE" = "#C76626"),
    breaks = c("TRUE", "FALSE"),
    labels = c("FALSE" = "No", "TRUE" = "Yes")
  ) +
  labs(
    x = "Essential niche count",
    y = "",
    color = "AT-rich"
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
# sd(deltaGC1_Intercept)     0.60      0.40     0.07     1.59 1.00     1874     2078
# sd(deltaGC2_Intercept)     0.14      0.15     0.00     0.51 1.00     3464     3440
# sd(deltaGC3_Intercept)     4.13      1.57     2.13     8.13 1.00     3264     4317
# sd(deltaGC4_Intercept)     5.40      2.06     2.82    10.74 1.00     2918     3864
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# deltaGC1_Intercept     1.90      0.32     1.21     2.47 1.00     3220     3671
# deltaGC2_Intercept   -16.00      0.12   -16.22   -15.74 1.00     5649     4424
# deltaGC3_Intercept    11.68      1.70     8.41    15.20 1.00     2520     3195
# deltaGC4_Intercept    15.52      2.27    11.00    20.19 1.00     2240     2912
# deltaGC1_count         0.13      0.15    -0.16     0.42 1.00     8285     6574
# deltaGC2_count        -0.13      0.13    -0.38     0.12 1.00    10389     6436
# deltaGC3_count        -0.44      0.22    -0.86    -0.02 1.00     5575     5701
# deltaGC4_count        -0.40      0.26    -0.92     0.12 1.00     5572     5709
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma_deltaGC1     4.60      0.07     4.47     4.74 1.00    11610     5733
# sigma_deltaGC2     4.09      0.06     3.97     4.21 1.00    10849     6371
# sigma_deltaGC3     6.77      0.10     6.58     6.97 1.00     5122     6088
# sigma_deltaGC4     8.32      0.12     8.08     8.56 1.00     5103     5930
# 
# Residual Correlations: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# rescor(deltaGC1,deltaGC2)     0.30      0.02     0.26     0.33 1.00    11818     6096
# rescor(deltaGC1,deltaGC3)     0.01      0.02    -0.03     0.05 1.00     5846     5342
# rescor(deltaGC2,deltaGC3)    -0.10      0.02    -0.14    -0.06 1.00     5728     5765
# rescor(deltaGC1,deltaGC4)     0.08      0.02     0.04     0.12 1.00     5719     5244
# rescor(deltaGC2,deltaGC4)    -0.12      0.02    -0.16    -0.08 1.00     5364     5887
# rescor(deltaGC3,deltaGC4)     0.93      0.00     0.92     0.93 1.00     5960     6636
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

