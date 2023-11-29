##Relative abundace graph 

#required packages 
options(width = 100)
require(tidyverse)
require(phyloseq)
require(microViz)
require(rstatix)
require(ggpubr)
require(ggplot2)
require(patchwork) # for arranging groups of plots
knitr::opts_chunk$set(fig.height = 6, fig.width = 9)

pscount2

ps <- pscount2

sample_data(ps)
psBASAL <- ps_filter(ps, Treatment == "Basal Diet")
sample_data(psBASAL)
# get example phyloseq data from corncob package and tidy up
psBASAL %>%
  tax_filter(min_prevalence = 2) %>%
  tax_fix() %>%
  phyloseq_validate()


hueRank <- "Phylum"
hueRankPlural <- "Phyla"
shadeRank <- "Family"

# Sort phyloseq at lower, and then higher ranks
pseq2 <- psBASAL %>%
  tax_fix() %>% 
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# Specify number of hues and shades desired
nHues <- 3 # "Other" phyla will be shades of grey
nShades <- 4 # "Other" families will be the lightest shade of each hue

hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(pseq2)[, hueRank]),
  shade = as.vector(tt_get(pseq2)[, shadeRank]),
  counts = taxa_sums(otu_get(pseq2))
)

hierarchicalPalInfo <- hierarchicalPalInfo %>%
  dplyr::mutate(
    hue = forcats::fct_other(
      f = hue, keep = unique(hue)[seq_len(nHues)],
      other_level = paste("Other", hueRankPlural)
    ),
    nChrHue = nchar(as.character(hue)), padHue = max(nChrHue) - nChrHue
  ) %>%
  dplyr::group_by(hue) %>%
  dplyr::mutate(
    shade = forcats::fct_other(
      f = shade, keep = unique(shade)[seq_len(nShades - 1)],
      other_level = "Other"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    nChrShade = nchar(as.character(shade)), padShade = max(nChrShade) - nChrShade,
    Taxa = paste0(hue, ": ", strrep(" ", padHue), shade, strrep(" ", padShade))
  )

hierarchicalPalMatrix <- matrix(
  data = sapply(
    X = seq(from = 30, to = 75, length.out = nShades),
    FUN = function(l) scales::hue_pal(l = l, h.start = 30)(n = nHues)
  ),
  byrow = TRUE, ncol = nHues
)
hierarchicalPalMatrix <- cbind(hierarchicalPalMatrix, grey.colors(n = nShades))

hierarchicalPal <- hierarchicalPalMatrix %>%
  as.vector() %>%
  setNames(unique(hierarchicalPalInfo$Taxa))


plotBASAL <- pseq2 %>%
  ps_get() %>%
  phyloseq::merge_samples(group = "Age_Days") %>% 
  tax_mutate("Phylum: Family" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    tax_level = "Phylum: Family", n_taxa = length(hierarchicalPal),
    tax_order = "asis", palette = hierarchicalPal, bar_width = 0.975
  ) 


p1 <- plotBASAL + scale_x_discrete(limits = c ("21", "20", "19", "18", "17", "16", "15", "14", "13","12","11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1")) + coord_flip() + ggtitle("Basal Diet") +
  theme(legend.text = element_text(family = "mono")) + # for text alignment  
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 10)) +
  theme(plot.title = element_text(size = 20, face = "bold"))


p2 <- plotProb + scale_x_discrete(limits = c ("21", "20", "19", "18", "17", "16", "15", "14", "13","12","11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1")) + coord_flip() + ggtitle("Probiotic") +
  theme(legend.text = element_text(family = "mono")) + # for text alignment  
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 10)) +
  theme(plot.title = element_text(size = 20, face = "bold"))

p3 <- plotEssential + scale_x_discrete(limits = c ("21", "20", "19", "18", "17", "16", "15", "14", "13","12","11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1")) + coord_flip() + ggtitle("Essential Oils") +
  theme(legend.text = element_text(family = "mono")) + # for text alignment  
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 10)) +
  theme(plot.title = element_text(size = 20, face = "bold"))

p4 <- plotATB + scale_x_discrete(limits = c ("21", "20", "19", "18", "17", "16", "15", "14", "13","12","11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1")) + coord_flip() + ggtitle("Antibiotic") +
  theme(legend.text = element_text(family = "mono")) + # for text alignment  
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 10)) +
  theme(plot.title = element_text(size = 20, face = "bold")) 

leg <- get_legend(p4)  
legenda <- as_ggplot(leg) 
ggsave("legenda.pdf")

plot1_2 <- ggarrange(p1, p2, legend = "none") 
plot3_4 <- ggarrange(p3, p4, legend = "none")
ggsave("plot1_2.pdf")
ggsave("plot3_4.pdf")
