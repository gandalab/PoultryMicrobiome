# Repeated Measures - Alpha Diversity Analysis
# April 20, 2023 - AF

# setup ----

# load packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(writexl)

# load phyloseq objects
# set directory 
setwd("/Users/anafonseca/Chicken_Microbiome_data/broiler_microbiome_groups/fastpfiles/asv-fastp/phyloseq_object/")

# load data
ps <- readRDS("broiler_microbiome_groups/fastpfiles/asv-fastp/phyloseq_object/ps-decontam-filtered-counts.rds")

ps
ps <-  subset_samples(ps, Treatment == "Basal Diet" | Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Essential Oils")
ps <-  subset_samples(ps, Age_Days == "1" | Age_Days == "10" | Age_Days == "21")
sample_data(ps)

##LOAD THE DIVERSITY TABLE

# create data frames ----
adivA <- data.frame(
  "Observed" = phyloseq::estimate_richness(ps, measures = "Observed"),
  "Simpson" = phyloseq::estimate_richness(ps, measures = "Simpson"),
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(ps)$Treatment,
  "Age_Days" = phyloseq::sample_data(ps)$Age_Days
)

##LOAD THE DIVERSITY TABLE
write_xlsx(adivA, "table_AlphaDiversity.xlsx")
alpha_table <- readxl::read_excel("broiler_microbiome_groups/fastpfiles/asv-fastp/phyloseq_object/table_AlphaDiversity.xlsx")
alpha_table

###NEW ALPHA DIVERSITY 
ShaDiv <- lm(Shannon ~ PROB + EO + ABX, data = alpha_table)
summary(ShaDiv)
anova(ShaDiv)

ObservedDiv <- lm(Observed ~ PROB + EO + ABX, data = alpha_table)
summary(ObservedDiv)
anova(ObservedDiv)


testAr <- aov(Shannon ~ Treatment*Age_Days, data = adivA)
summary(testAr)

# repeated measures anova ----
# https://m-clark.github.io/docs/mixedModels/anovamixed.html 

anovaA <- aov(Shannon ~ Treatment*Age_Days, data = adivA)
summary(anovaA) #  F = 1.939, P = 0.0842, Df = 6

# age came up significant so I checked just for curiositys sake
testage <- aov(Shannon ~ Age_Days, data = adivA)
summary(testage) # 9.99e-07 ***


TukeyHSD(testage)

#$Age_Days
#diff        lwr        upr     p adj
#10-1   0.2356566 -0.1669311  0.6382442 0.3479122
#21-1  -0.6942969 -1.0968846 -0.2917093 0.0002540
#21-10 -0.9299535 -1.3357237 -0.5241833 0.0000012

## box plots -----
sample_data(ps)$"Age_Days" <- factor(sample_data(ps)$"Age_Days", 
                                     levels = c("1", "10", "21")) 
sample_data(ps)$"Treatment" <- factor(sample_data(ps)$"Treatment", 
                                      levels = c("Basal Diet", "Probiotic", "Essential Oils", "BMD"))

A <- plot_richness(ps, x="Age_Days", measures=c("Shannon"), title = "A", color = "Treatment") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#FF3366","#663366","#333399","#9999CC")) + 
  theme_classic() + 
  theme(legend.position = "none")

sample_data(ps)$"Age_Days" <- factor(sample_data(ps)$"Age_Days", 
                                     levels = c("1", "10", "21"))
sample_data(ps)$"Treatment" <- factor(sample_data(ps)$"Treatment", 
                                      levels = c("Basal Diet", "Probiotic", "Essential Oils", "BMD"))

A

ggsave(filename = "alpha-1-10-21days.pdf", dpi = 600)

# lineplot ----

adivA$"Age_Days" <- factor(adivA$"Age_Days", 
                           levels = c("1", "10", "21"))

adivA$"Treatment" <- factor(adivA$"Treatment", 
                            levels = c("Basal Diet", "Probiotic", "Essential Oils", "BMD"))


adivA %>% 
  ggline(x = "Age_Days", y = "Shannon", group = "Treatment", color = "Treatment", palette = c("#FF3366","#663366","#333399","#9999CC"), add = "mean_sd", error.plot = "errorbar")

ggsave(filename = "geomline-plot-treatment*age-time.pdf", dpi = 600)

##_______________________________________________________________________________________##

###USING THE ALL DAYS ###

# load data
ps <- readRDS("ps-decontam-filtered-counts.rds")

ps
ps <-  subset_samples(ps, Treatment == "Basal Diet" | Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Essential Oils")
sample_data(ps)

class(alpha$Age_Days)

#transforme Age_Days as a numeric 
alpha$Age_Days <- as.numeric(alpha$Age_Days)

#Check if it worked 
str(alpha)

# create data frames ----
alldays <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "Observed" = phyloseq::estimate_richness(ps, measures = "Observed"),
  "Simpson" = phyloseq::estimate_richness(ps, measures = "Simpson"),
  "Treatment" = phyloseq::sample_data(ps)$Treatment,
  "Age_Days" = phyloseq::sample_data(ps)$Age_Days
)

testageall <- aov(Shannon ~ Treatment*Age_Days, data = alldays)
summary(testageall)

# repeated measures anova ----
# https://m-clark.github.io/docs/mixedModels/anovamixed.html 

class(alldays$Age_Days)
#transforme Age_Days as a numeric 
alldays$Age_Days <- as.numeric(alldays$Age_Days)

#Check if it worked 
str(alpha)

anovaalldays <- aov(Shannon ~ Treatment*alldays$Age_Days, data = alldays)
summary(anovaalldays) #  F = 1.068, P = 0.345, Df = 60
summary(lm(Shannon ~ Treatment*Age_Days, data = alldays))
matrix_coef <- summary(lm(Shannon ~ Treatment*Age_Days, data = alldays))$coefficients
matrix_coef

# age came up significant so I checked just for curiositys sake
testageall <- aov(Shannon ~ Age_Days, data = alldays)
summary(testageall) # <2e-16 ***

ageallhsd <- TukeyHSD(testageall)

library(gt)

agellhsd.df <- as.data.frame(ageallhsd$Age_Days)

age.hsd.table <- gt(agellhsd.df, rownames_to_stub = TRUE)

age.hsd.table

write.csv(age.hsd.table, "newalphadiversity_allages_hsdresults.csv", row.names=TRUE)

## box plots -----
sample_data(ps)$"Age_Days" <- factor(sample_data(ps)$"Age_Days", 
                                     levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
                                                "16", "17", "18", "19", "20", "21")) 
sample_data(ps)$"Treatment" <- factor(sample_data(ps)$"Treatment", 
                                      levels = c("Basal Diet", "Probiotic", "Essential Oils", "BMD"))

allagesplot <- plot_richness(ps, x="Age_Days", measures=c("Shannon"), title = "A", color = "Treatment") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#FF3366","#663366","#333399","#9999CC")) + 
  theme_classic() + 
  theme(legend.position = "none")

sample_data(ps)$"Age_Days" <- factor(sample_data(ps)$"Age_Days", 
                                     levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
                                                "16", "17", "18", "19", "20", "21")) 
sample_data(ps)$"Treatment" <- factor(sample_data(ps)$"Treatment", 
                                      levels = c("Basal Diet", "Probiotic", "Essential Oils", "BMD"))

allagesplot

ggsave(filename = "alpha-div-allages.pdf", dpi = 600)

# lineplot ----

alldays$"Age_Days" <- factor(alldays$"Age_Days", 
                             levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
                                        "16", "17", "18", "19", "20", "21")) 

alldays$"Treatment" <- factor(alldays$"Treatment", 
                              levels = c("Basal Diet", "Probiotic", "Essential Oils", "BMD"))


alldays %>% 
  ggline(x = "Age_Days", y = "Shannon", group = "Treatment", color = "Treatment", palette = c("#FF3366","#663366","#333399","#9999CC"), add = "mean_sd", error.plot = "errorbar")

ggsave(filename = "geomline-plot-treatment*allagedays-time.pdf", dpi = 600)


