## Creating a phyloseq object

#Ana 07/20/22

# install (if necessary) and load packages
require(tidyverse)
#install phyloseq
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("phyloseq")
require(phyloseq)
#install DADA2
#install.packages("devtools")
library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.16")
require(dada2)
# install latest version of microViz
#devtools::install_github("david-barnett/microViz@0.7.1")
require(microViz)
require(BiocManager)
#BiocManager::install("decontam")
require(decontam)
require(ggpubr)

# set directory
setwd("/Users/anafonseca/OneDrive - The Pennsylvania State University/Documents/Chicken_Microbiome_data/broiler_microbiome_groups/fastpfiles/asv-fastp/phyloseq_object/")

### ---- phyloseq ----
# LOAD ASV table
load("asv-table.RData")

# get PCR-neg
grep("PCR-neg", rownames(seqtab.nochim))

# get row 707
neg <- seqtab.nochim[707,]

# range of taxa presence
range(neg)

### REMOVE PCR-NEG
seqtab.noneg <- seqtab.nochim[-707,]

# load tax table 
load("tax-table.RData")

# matrices in R 
dim(seqtab.nochim)
dim(tax)

# create phyloseq object 
ps <- phyloseq(otu_table(t(seqtab.noneg), taxa_are_rows = TRUE),
               tax_table(tax))

# look at the phyloseq object 
ps

# get number of taxa
ntaxa(ps)

#get taxa ranks
rank_names(ps)

# access the data "slots" with @
head(ps@tax_table)
head(ps@otu_table)

# remove all the duplicates and "play" samples
to_remove <- c("116_S116_r1_fastp.fq","310_S41_r1_fastp.fq","629_S215_r1_fastp.fq","709-1_S86_r1_fastp.fq","93-2_S162_r1_fastp.fq","38_S6_r1_fastp.fq","49_S252_r1_fastp.fq","74_S53_r1_fastp.fq",
               "109_S160_r1_fastp.fq","144_S260_r1_fastp.fq","154_S270_r1_fastp.fq","256_S18_r1_fastp.fq",
               "290_S16_r1_fastp.fq","324_S104_r1_fastp.fq","392_S149_r1_fastp.fq","426_S107_r1_fastp.fq",
               "460_S74_r1_fastp.fq","494_S63_r1_fastp.fq","528_S274_r1_fastp.fq","562_S88_r1_fastp.fq","664_S72_r1_fastp.fq")
newps <- prune_samples(!(sample_names(ps) %in% to_remove), ps)
newps
ps <- newps
ps

# fix ASV names 
### from dada2 tutorial: fix ASV names
dna <- Biostrings::DNAStringSet(taxa_names(newps))
names(dna) <- taxa_names(newps)
newDNAps <- merge_phyloseq(newps, dna)
taxa_names(newDNAps) <- paste0("ASV", seq(ntaxa(newDNAps)))
newDNAps

# get taxa names 
head(taxa_names(newDNAps))
head(sample_names(newDNAps))

# re-name
ps <- newDNAps

## ---- metadata ----
# read in metadata file 
dat <- readxl::read_xlsx(path = "/Users/anafonseca/Chicken_Microbiome_data/broiler_microbiome_groups/fastpfiles/asv-fastp/updated_chicken_microbiome_metadata.xlsx",
                         sheet = 1)

# look at data
head(dat)
str(dat)

## REMOVE PCR-neg
dat <- dat %>% filter(!Project_ID == "PCR-neg")
tail(dat)
str(dat)

# fix sample names to get ONLY the sample ID
#names <- sapply(str_split(sample_names(ps), "_"), `[`, 1) # I had an error with telling me about duplicates 
names <- str_remove_all(sample_names(ps), "_S(\\d){1,3}_r1_fastp.fq")

# figure out the duplicates 
length(unique(names))
#namedf <- data.frame(id = names)
#namedf %>% group_by(id) %>% filter(n() > 1)
#which(names == "280")
#sample_names(ps)[194]
#sample_names(ps)[195]


# change sample names to NAMES
sample_names(ps) <- names

# format our data to add to phyloseq
sampdf <- dat %>%
  column_to_rownames(var = "Project_ID")


#row.names(dat) <-dat$Sample_ID

#names <- rownames(seqtab.nochim) #we are making a list 
#gsub("_.*","", names) # removing file extension
#names <- gsub("_.*","", names)
#rownames(seqtab.nochim) <- names

# add to phyloseq
#ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows= FALSE),
#              sample_data(sampdf), 
#             tax_table(tax))

## build phyloseq
sample_data(ps) <- sampdf

# did it work?
#ps

#rownames(seqtab.nochim) %in% rownames(sample_data(ps))
## save our output
# re-name our phyloseq

psraw <- ps

# save as Rimage
save(psraw, file = "ps-raw.rds")

#Check bacillus presence
psraw
sample_data(psraw)
ps <-  subset_samples(psraw, Treatment == "Basal Diet" | Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Essential Oils")
ps

grep("Bacillus", ps@tax_table[,"Genus"])
bac <- subset_taxa(ps, Genus == "Bacillus")
bac@tax_table
unique(ps@tax_table[,"Genus"])


