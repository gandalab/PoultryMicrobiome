##  Beta Analysis
## AF - Last Edited April 20-2023

# ---- setup ----
# load packages
library(BiocManager)
library(phyloseq)
library(vegan)
library(patchwork)
library(microViz)
library(ggplot2)
library(dplyr)
library(stringr)
library(pairwiseAdonis)
library(scales)

set.seed(200789)

# load phyloseq objects
# set directory 
setwd("/Users/anafonseca/Chicken_Microbiome_data/broiler_microbiome_groups/fastpfiles/asv-fastp/phyloseq_object/")

# load data
ps <- readRDS("ps-decontam-filtered-relabun.rds")

ps
ps <-  subset_samples(ps, Treatment == "Basal Diet" | Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Essential Oils")
#ps <-  subset_samples(ps, Age_Days == "1" | Age_Days == "10" | Age_Days == "21")
sample_data(ps)

# check if we have NAs
anyNA(tax_table(ps)[,"Phylum"])
anyNA(tax_table(ps)[,"Phylum"])

# tax fix our phyloseq object
ps <- tax_fix(ps)

ps <- ps %>% tax_fix(unknowns = c("Incertae Sedis"))

# put the time series in order
#sample_data(ps)$"Age_Days" <- factor(sample_data(ps)$"Age_Days", 
#                                     levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
#                                                "16", "17", "18", "19", "20", "21")) 
sample_data(ps)$"Treatment" <- factor(sample_data(ps)$"Treatment", 
                                      levels = c("Basal Diet", "Probiotic", "Essential Oils", "BMD"))


## ---- Setup Treatments to be Antibiotic = BMD ----

ps <- ps %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Antibiotic"), true = "Antibiotic", false = "BMD")
  ) 


### ---- BETA DIVERSITY ----

# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")

#ADONIS test
# age, treatment, and their interaction were included in the model


vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Treatment*sample_data(beta1)$Age_Days, permutations = 10000) 
# interaction -> p = 0.6130387 
# age -> p = 9.999e-05 ***
# treatment -> p = 0.0005999 ***

#figure out cols (4 of them) for Age_Days
classcols <- colorRampPalette(c("#000033", "white"))
plot(rep(1,30), col = classcols(30), pch = 19, cex = 3) # keep cols 2,5,9,11
#get color numbers
classcols(30)
classcols <- c("1" = "#000033",
               "2" = "#08083A",
               "3" = "#111141",
               "4" = "#1A1A48",
               "5" = "#23234F",
               "6" = "#2B2B56",
               "7" = "#34345D", 
               "8" = "#3D3D64",
               "9" = "#46466B",
               "10" = "#4F4F72", 
               "11" = "#7F7F99",
               "12" = "#8C8CA3",
               "13" = "#696987",
               "14" = "#72728E",
               "15" = "#7B7B95",
               "16" = "#83839C",
               "17" = "#8C8CA3",
               "18" = "#9595AA",
               "19" = "#9E9EB1", 
               "20" = "#A7A7B8",
               "21" = "#AFAFBF")

# transforming variable Age_Days from factor to numeric 
#https://github.com/joey711/phyloseq/issues/201
class(sample_data(ps)$Age_Days)           
sample_data(ps)$Age_Days <- as(sample_data(ps)$Age_Days, "numeric")


##  PCA plot 
p1 <- ps %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Age_Days", plot_taxa = 1:5, size = 4, shape = "Treatment") +
  scale_fill_continuous() +
  theme_classic() +
  ggtitle("PCA ") +
  #theme(legend.position = "none") +
  labs(caption = "")
p1 
ggsave(filename = "PCA-microViz-allagextreatment.pdf", dpi = 600)

p1 + scale_colour_gradient(low = "brown", high = "lightgrey")
ggsave(filename = "PCA-microViz-allagextreatment_PLOT2.pdf", dpi = 600)

p1 + scale_colour_gradient(low = "black", high = "lightgrey")
ggsave(filename = "PCA-microViz-allagextreatment_PLOT3.pdf", dpi = 600)

p1 + scale_colour_gradientn(colours = c("#000033", "#2B2B56", "#4F4F72", "#9E9EB1", "#AFAFBF"))
ggsave(filename = "PCA-microViz-allagextreatment_PLOT4.pdf", dpi = 600)
#https://r-graphics.org/recipe-colors-palette-continuous

##______________________________________________________________##

# load phyloseq objects
# load data
ps <- readRDS("ps-decontam-filtered-relabun.rds")

ps
ps <-  subset_samples(ps, Treatment == "Basal Diet" | Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Essential Oils")

###Subset Age_Days samples 
ps <-  subset_samples(ps, Age_Days == "1" | Age_Days == "2" | Age_Days == "3" | Age_Days == "5" | Age_Days == "7" | Age_Days == "10" | Age_Days == "15" |Age_Days == "21")

sample_data(ps)$"Age_Days" <- factor(sample_data(ps)$"Age_Days", 
                                     levels = c("1", "2", "3", "5", "7", "10","15","21")) 
#plot
p2 <- ps %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Age_Days", plot_taxa = 1:5, size = 2, shape = "Treatment") +
  scale_color_manual(values = c("#028571", "#669999", "#009999", "#dfc27d", "#FFCC33", "#a76119", "#CC6633", "#543005")) +
  stat_ellipse(aes(group = Age_Days, color = Age_Days)) +
  theme_classic() +
  ggtitle("B") +
  # theme(legend.position = "none") +
  labs(caption = "")
p2 

ggsave(filename = "PCA-microViz-subsetedages.pdf", dpi = 600)

# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")
vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Treatment*sample_data(beta1)$Age_Days, permutations = 999)
vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Age_Days, permutations = 999) # p=0.001 ***

# Pairwise comparison using pairwiseAdonis 
data <- phyloseq::sample_data(beta1)
data
pairwiseAdonis::pairwise.adonis(psdist, data$Age_Days)

# Pairwise comparison using pairwiseAdonis 
data <- phyloseq::sample_data(beta1)
data
pairwiseAdonis::pairwise.adonis(psdist, data$Treatment)


###Subset for 1, 10 and 21 days 

# load phyloseq objects
# load data
ps <- readRDS("ps-decontam-filtered-relabun.rds")

ps
ps <-  subset_samples(ps, Treatment == "Basal Diet" | Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Essential Oils")

###Subset Age_Days samples 
ps <-  subset_samples(ps, Age_Days == "1" | Age_Days == "10" |Age_Days == "21")

sample_data(ps)$"Age_Days" <- factor(sample_data(ps)$"Age_Days", 
                                     levels = c("1", "10","21")) 
#plot
p2 <- ps %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Age_Days", size = 5, shape = "Treatment") +
  scale_color_manual(values = c("#028571", "#CC6633", "#543005")) +
  stat_ellipse(aes(group = Age_Days, color = Age_Days)) +
  theme_classic() +
  ggtitle("B") +
  # theme(legend.position = "none") +
  labs(caption = "")
p2 

ggsave(filename = "PCA-microViz-subsetedages1-10-21.pdf", dpi = 600)

# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")
vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Treatment*sample_data(beta1)$Age_Days, permutations = 999)
#p age 0.001
vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Age_Days, permutations = 999) # p=0.001 ***

# Pairwise comparison using pairwiseAdonis 
data <- phyloseq::sample_data(beta1)
data
pairwiseAdonis::pairwise.adonis(psdist, data$Age_Days)

#      pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1 vs 10  1  2899.685 34.86673 0.3637000   0.001      0.003   *
# 1 vs 21  1  4437.285 48.75890 0.4442364   0.001      0.003   *
# 10 vs 21  1  1482.811 10.86316 0.1532977   0.001      0.003   *

ps %>% 
  tax_transform("clr", rank = "Genus") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>% 
  ord_plot(color = "Age_Days", shape = "Treatment", size = 3) +
  scale_color_viridis_d()

