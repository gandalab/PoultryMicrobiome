## differential rel abundance - log2 linear model 

require(microViz)
require(tidyverse)
require(phyloseq)

# set directory 
setwd("/Users/anafonseca/Chicken_Microbiome_data/broiler_microbiome_groups/fastpfiles/asv-fastp/phyloseq_object/")

# load data
ps <- readRDS("ps-decontam-filtered-counts.rds")

# check phyloseq
ps

#remove baseline 
ps <-  subset_samples(ps, Treatment == "Basal Diet" | Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Essential Oils")

# check
sample_data(ps)
ps

ps <- ps %>% 
  tax_fix()

## ---- test1: basal diet v probiotic ----

# subset 
phylo <- ps %>% 
  ps_filter(Treatment %in% c("Basal Diet", "Probiotic")) 

# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Treatment")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

# adjust for multiple comparisons (within the test)
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

## DOUBLE CHECK BY HAND TO SEE WHICH IS THE REFERENCE LEVEL
check <- phylo %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% ps_get() %>%  ps_melt()

# plot just Bifidobacterium
forplot <- check %>% 
  filter(str_detect(Genus, "Bifidobacterium")) %>% 
  filter(Treatment %in% c("Basal Diet", "Probiotic"))

ggboxplot(forplot, x = "Treatment", y = "Abundance")

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.BH.rank < 0.05])

# save for comparison
mod1 <- lm_stats

##THERE WAS A SIGNIFICANCE BTW PROBIOTIC AND BASAL DIET
#1 TreatmentProbiotic G: Bifidobacterium Genus    -1.44     0.399     -3.61 0.000354        0.0258
#[1] "G: Bifidobacterium"


## ---- test2: difference between Basal Diet v Essential Oils ----

# subset 
phylo <- ps %>% 
  ps_filter(Treatment %in% c("Basal Diet", "Essential Oils")) 

# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Treatment")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

# adjust for multiple comparisons (within the test)
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.BH.rank < 0.05])

# save for comparison
mod2 <- lm_stats

##THERE WAS NO SIGNIFICANCE BTW BASAL DIET AND ESSENTIAL OILS 

## ---- test3: difference between Basal Diet v BMD ----

# subset 
phylo <- ps %>% 
  ps_filter(Treatment %in% c("Basal Diet", "BMD")) 

# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Treatment")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

# adjust for multiple comparisons (within the test)
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.BH.rank < 0.05])

# save for comparison
mod3 <- lm_stats

##THERE WAS SIGNIFICANCE BTW BASAL DIET AND EBMD 
#1 TreatmentBMD G: Bifidobacterium   Genus   -1.58      0.311     -5.09 0.000000610     0.0000464
#2 TreatmentBMD G: Ligilactobacillus Genus   -1.60      0.378     -4.24 0.0000295       0.00112  
#3 TreatmentBMD G: Intestinimonas    Genus    0.564     0.154      3.66 0.000294        0.00744  

#"G: Bifidobacterium"   "G: Ligilactobacillus" "G: Intestinimonas"   


## ---- save tables for summary stats ----

# model 1
m1 <- mod1$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "Basal Diet versus Probiotic") %>% 
  ungroup() %>% 
  select(-c(taxon, p.adj.BH.rank, rank, term)) %>% 
  relocate(model, Genus)

write.table(m1, file = "rel-abund-basal diet_probiotic-results.txt", sep = "\t", row.names = FALSE)

# model 2
m2 <- mod2$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "Basal Diet versus Essential Oils") %>% 
  ungroup() %>% 
  select(-c(taxon, p.adj.BH.rank, rank, term)) %>% 
  relocate(model, Genus)

write.table(m2, file = "rel-abund-basal diet_essential oils-results.txt", sep = "\t", row.names = FALSE)

# model 3
m3 <- mod3$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "Basal Diet versus BMD") %>% 
  ungroup() %>% 
  select(-c(taxon, p.adj.BH.rank, rank, term)) %>% 
  relocate(model, Genus)

write.table(m3, file = "rel-abund-Basal Diet_BMD-results.txt", sep = "\t", row.names = FALSE)

## ---- test4: overall difference  ----

# subset 
phylo <- ps %>% 
  ps_filter(Treatment %in% c("Basal Diet", "Probiotic", "BMD", "Essential Oils")) 

# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Treatment")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

# adjust for multiple comparisons (within the test)
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.BH.rank < 0.05])

# save for comparison
mod4 <- lm_stats

# model 4
m4 <- mod4$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "All") %>% 
  ungroup() %>% 
  select(-c(taxon, p.adj.BH.rank, rank, term)) %>% 
  relocate(model, Genus)

write.table(m4, file = "rel-abund-all-results.txt", sep = "\t", row.names = FALSE)

------###Compare the experimental group with the BMD (positive control)###-----

# subset 
phylo <- ps %>% 
  ps_filter(Treatment %in% c("BMD", "Probiotic")) 

# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Treatment")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

# adjust for multiple comparisons (within the test)
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

## DOUBLE CHECK BY HAND TO SEE WHICH IS THE REFERENCE LEVEL
check <- phylo %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% ps_get() %>%  ps_melt()

# plot just Bifidobacterium
forplot <- check %>% 
  filter(str_detect(Genus, "Bifidobacterium")) %>% 
  filter(Treatment %in% c("BMD", "Probiotic"))

ggboxplot(forplot, x = "Treatment", y = "Abundance")

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.BH.rank < 0.05])

# save for comparison
mod5 <- lm_stats

##THERE WAS NO SIGNIFICANCE BTW PROBIOTIC AND BMD

#BMD vs. EO
# subset 
phylo <- ps %>% 
  ps_filter(Treatment %in% c("BMD", "Essential Oils")) 

# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Treatment")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

# adjust for multiple comparisons (within the test)
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

## DOUBLE CHECK BY HAND TO SEE WHICH IS THE REFERENCE LEVEL
check <- phylo %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% ps_get() %>%  ps_melt()

# plot just Bifidobacterium
forplot <- check %>% 
  filter(str_detect(Genus, "Bifidobacterium")) %>% 
  filter(Treatment %in% c("BMD", "Essential Oils"))

ggboxplot(forplot, x = "Treatment", y = "Abundance")

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.BH.rank < 0.05])

# save for comparison
mod6 <- lm_stats

##THERE WAS NO SIGNIFICANCE BTW ESSENTIAL OILS AND BMD


## ---- save tables for summary stats ----
# model 5
m5 <- mod5$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "BMDvsPB") %>% 
  ungroup() %>% 
  select(-c(taxon, p.adj.BH.rank, rank, term)) %>% 
  relocate(model, Genus)

write.table(m5, file = "rel-abund-BMDvsPB-results.txt", sep = "\t", row.names = FALSE)

# model 6
m6 <- mod6$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "BMDvsEO") %>% 
  ungroup() %>% 
  select(-c(taxon, p.adj.BH.rank, rank, term)) %>% 
  relocate(model, Genus)

write.table(m6, file = "rel-abund-BMDvsEO-results.txt", sep = "\t", row.names = FALSE)

