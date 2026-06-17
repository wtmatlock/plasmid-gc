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
    # set factor levels
    delta_type = factor(delta_type, levels = c(
      "First codon position",
      "Second codon position",
      "Third codon position",
      "Four-fold degenerate sites"
    ))
  )

# make figure
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

summary(fit)
