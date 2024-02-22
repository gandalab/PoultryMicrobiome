##Alpha Diversity Rarefied
#Alpha Diversity Analysis
# Feb 20, 2024 - AF

# setup ----
# Set seed for reproducibility
set.seed(123)

# load packages
library(phyloseq)
library(vegan)
library(microViz)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(writexl)
library(car)
#install.packages("dunn.test")
library(dunn.test)

# load pcar# load phyloseq objects
# set directory 
setwd("/Users/aff30/OneDrive - The Pennsylvania State University/Documents/Chicken_Microbiome_data/broiler_microbiome_groups/fastpfiles/asv-fastp/phyloseq_object/")

# load data
ps <- readRDS("ps-decontam-filtered-counts.rds")

ps

#sample depth before rarefy

# Extract sample sums
sample_depths <- sample_sums(ps)

# Convert to a data frame for ggplot2
sample_depths_df <- data.frame(Sample = names(sample_depths), Depth = sample_depths)

# Plot
ggplot(sample_depths_df, aes(x = Sample, y = Depth)) +
  geom_bar(stat = "identity") +
  theme_minimal() 
  
# First, determine a suitable rarefaction depth. This is often the sample size of the smallest sample, but you can choose another threshold based on your data.

sample_sizes <- sample_sums(ps)
rarefaction_depth <- min(sample_sizes[sample_sizes > 0]) # Ensure no zero-size samples

# Rarefy to even depth
ps_rarefied <- rarefy_even_depth(ps, sample.size = rarefaction_depth, rngseed = 123, verbose = TRUE)

##check after rarefy
# Extract sample sums
sample_depths_rarefied <- sample_sums(ps_rarefied)

# Convert to a data frame for ggplot2
sample_depths_df_rarefied <- data.frame(Sample = names(sample_depths_rarefied), Depth = sample_depths_rarefied)

# Plot
ggplot(sample_depths_df_rarefied, aes(x = Sample, y = Depth)) +
  geom_bar(stat = "identity") +
  theme_minimal() 

# Calculate observed species (OTUs) before rarefaction
obs_species_before <- estimate_richness(ps, measures = "Observed")

# Calculate observed species (OTUs) after rarefaction
obs_species_after <- estimate_richness(ps_rarefied, measures = "Observed")

# Combine and prepare for plotting
obs_species_combined <- data.frame(
  Sample = rep(sample_names(ps), 2),
  Richness = c(obs_species_before$Observed, obs_species_after$Observed),
  Condition = rep(c("Before", "After"), each = nsamples(ps))
)

# Plot
ggplot(obs_species_combined, aes(x = Sample, y = Richness, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Observed Species Richness Before and After Rarefaction",
       x = "Sample",
       y = "Observed Species")


ps <-  subset_samples(ps_rarefied, Treatment == "Basal Diet" | Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Essential Oils")

#ps <-  subset_samples(ps, Age_Days == "1" | Age_Days == "10" | Age_Days == "21")
sample_data(ps)
ps

##LOAD THE DIVERSITY TABLE

# create data frames ----
adivA <- data.frame(
  "Observed" = phyloseq::estimate_richness(ps, measures = "Observed"),
  "Simpson" = phyloseq::estimate_richness(ps, measures = "Simpson"),
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "Chao1" = phyloseq::estimate_richness(ps, measures = "Chao1"),
  "Treatment" = phyloseq::sample_data(ps)$Treatment,
  "Age_Days" = phyloseq::sample_data(ps)$Age_Days
)


# check the normality of the data
#Density plot and Q-Q plot can be used to check normality visually.
#Density plot: the density plot provides a visual judgment about whether the distribution is bell shaped.

ggdensity(adivA$Shannon, 
          main = "Density plot of Shannon",
          xlab = "Shannon")

#Q-Q plot: Q-Q plot (or quantile-quantile plot) draws the correlation between a given sample and the normal distribution. A 45-degree reference line is also plotted.
ggqqplot(adivA$Shannon)


sh <- ggplot(adivA, aes(x=Shannon)) + 
  geom_histogram(fill="blue", color="black") + 
  ggtitle("Histogram of Data") +
  xlab("Shannon") +
  ylab("Frequency") + 
  theme_bw()

obs <- ggplot(adivA, aes(x=Observed)) + 
  geom_histogram(fill="blue", color="black") + 
  ggtitle("Histogram of Data") +
  xlab("Observed") +
  ylab("Frequency") 

simp <- ggplot(adivA, aes(x=Simpson)) + 
  geom_histogram(fill="blue", color="black") + 
  ggtitle("Histogram of Data") +
  xlab("Simpson") +
  ylab("Frequency") 

chao <- ggplot(adivA, aes(x=Chao1.Chao1)) + 
  geom_histogram(fill="blue", color="black") + 
  ggtitle("Histogram of Data") +
  xlab("Chao1") +
  ylab("Frequency") 

# Perform Shapiro-Wilk test
shapiro_test <- shapiro.test(adivA$Shannon)

# Show the results of the Shapiro-Wilk test
print(shapiro_test)

#From the output, the p-value < 0.05 implying that the distribution of the data are significantly different from normal distribution. In other words, we can assume the data is not normal.

#since it is not normal lets do the kruskal_test 
kruskal_test <- kruskal.test(Shannon ~ Treatment, data = adivA)
print(kruskal_test)

kruskal_test <- kruskal.test(Observed ~ Treatment, data = adivA)
print(kruskal_test)

#Check age 
kruskal_test <- kruskal.test(Shannon ~ Age_Days, data = adivA)
print(kruskal_test)
