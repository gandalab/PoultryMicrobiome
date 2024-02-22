# Repeated Measures - Alpha Diversity Analysis
# Feb 4, 2024 - AF

# setup ----
# Set seed for reproducibility
set.seed(123)

# load packages
library(phyloseq)
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
ps <-  subset_samples(ps, Treatment == "Basal Diet" | Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Essential Oils")


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

##LOAD THE DIVERSITY TABLE
write_xlsx(adivA, "table_AlphaDiversity.xlsx")
alpha_table <- readxl::read_excel("table_AlphaDiversity.xlsx")
alpha_table

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


# Perform Shapiro-Wilk test
shapiro_test_obs <- shapiro.test(adivA$Observed)

# Show the results of the Shapiro-Wilk test
print(shapiro_test_obs)


#From the output, the p-value < 0.05 implying that the distribution of the data are significantly different from normal distribution. In other words, we can assume the data is not normal.

#since it is not normal lets do the kruskal_test 
kruskal_test <- kruskal.test(Shannon ~ Treatment, data = adivA)
print(kruskal_test)

kruskal_test <- kruskal.test(Observed ~ Treatment, data = adivA)
print(kruskal_test)

#interaction
adivA$interaction <- with(adivA, paste(Treatment, Age_Days, sep = "_"))

test_inter <- kruskal.test(Shannon ~ interaction, data = adivA)
test_inter

# Post-hoc pairwise comparisons
post_hoc_result <- pairwise.wilcox.test(adivA$Shannon, adivA$interaction, p.adjust.method = "BH")
print(post_hoc_result)

#test the interaction for Observed 
test_inter_ob <- kruskal.test(Observed ~ interaction, data = adivA)
test_inter_ob

# Post-hoc pairwise comparisons for observed 
post_hoc_result_ob <- pairwise.wilcox.test(adivA$Observed, adivA$interaction, p.adjust.method = "BH")
print(post_hoc_result_ob)


## plots -----
sample_data(ps)$"Age_Days" <- factor(sample_data(ps)$"Age_Days", 
                                     levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
                                                "16", "17", "18", "19", "20", "21")) 
sample_data(ps)$"Treatment" <- factor(sample_data(ps)$"Treatment", 
                                      levels = c("Basal Diet", "Probiotic", "Essential Oils", "BMD"))


#group by treatment violin plots 
#colors
trts <- c("Basal Diet" = "#FF3366", 
          "Probiotic" = "#9999CC",
          "Essential Oils" = "#333399",
          "BMD" = "#663366")


AGE <- c("1" ="#028571", 
         "10" = "#CC6633", 
         "21" = "#543005")


shannonviol <- ggplot(adivA,#filter(treatment == "Control"),
                      aes(x = Treatment, y = Shannon, fill = Treatment)) +
  geom_boxplot(trim = FALSE, alpha = 0.8)  +
  scale_fill_manual(values = trts) +
  xlab("") +
  geom_point(position = "jitter", color = "black", size = 3) +
  theme_bw() + 
  labs(title = "A",
       x = "", 
       y = "Shannon's Index",
       shape = "Treatment",
       fill = "Treatment") + 
  #scale_x_discrete(name="", labels = c("cloaca" = "Cloaca", "sock"="Sock")) +
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none")

shannonviol 


obsviol <- 
  ggplot(adivA,#filter(treatment == "Control"),
         aes(x = Treatment, y = Observed, fill = Treatment)) +
  geom_boxplot(trim = FALSE, alpha = 0.8)  +
  scale_fill_manual(values = trts) +
  xlab("") +
  geom_point(position = "jitter", color = "black", size = 3) +
  theme_bw() + 
  labs(title = "B",
       x = "", 
       y = "Observed ASVs",
       shape = "Treatment",
       fill = "Treatment") + 
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.0)) +
  guides(fill = "none")

obsviol 



alphdiv_prelim_boxplot <- ggarrange(shannonviol, obsviol, ncol = 3, common.legend = TRUE, legend = "right" )

pdf("New_graphs/alphadiv_prelim_boxplot.pdf", width = 18, height = 12)
alphdiv_prelim_boxplot
dev.off()

alphdiv_prelim_violin_article <- ggarrange(shannonviol, obsviol, ncol = 3, common.legend = TRUE, legend = "right" )

pdf("New_graphs/alphadiv_prelim_violin_article.pdf", width = 18, height = 12)
alphdiv_prelim_violin_article
dev.off()

##----------------------------------------------------------------##
## Subsetting age 1, 10, 21 for the analysis 

# load data
ps <- readRDS("ps-decontam-filtered-counts.rds")

ps
ps <-  subset_samples(ps, Treatment == "Basal Diet" | Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Essential Oils")

ps <-  subset_samples(ps, Age_Days == "1" | Age_Days == "10" | Age_Days == "21")
sample_data(ps)

ps

# create data frames ----
adivA_age <- data.frame(
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

ggdensity(adivA_age$Shannon, 
          main = "Density plot of Shannon",
          xlab = "Shannon")

#Q-Q plot: Q-Q plot (or quantile-quantile plot) draws the correlation between a given sample and the normal distribution. A 45-degree reference line is also plotted.
ggqqplot(adivA_age$Shannon)


sh <- ggplot(adivA_age, aes(x=Shannon)) + 
  geom_histogram(fill="blue", color="black") + 
  ggtitle("Histogram of Data") +
  xlab("Shannon") +
  ylab("Frequency") + 
  theme_bw()
sh

obs <- ggplot(adivA_age, aes(x=Observed)) + 
  geom_histogram(fill="blue", color="black") + 
  ggtitle("Histogram of Data") +
  xlab("Observed") +
  ylab("Frequency") 
obs

# Perform Shapiro-Wilk test
shapiro_test_age <- shapiro.test(adivA_age$Shannon)

# Show the results of the Shapiro-Wilk test
print(shapiro_test_age)

#From the output, the p-value < 0.05 implying that the distribution of the data are significantly different from normal distribution. In other words, we can assume the data is not normal.

#since it is not normal lets do the kruskal_test 
#Check age 
kruskal_test <- kruskal.test(Shannon ~ Age_Days, data = adivA_age)
print(kruskal_test)

kruskal_testOb <- kruskal.test(Observed ~ Age_Days, data = adivA_age)
print(kruskal_testOb)

#age is different 
# Perform Dunn's test for multiple comparisons
# Adjusting the p-values using the Benjamini-Hochberg method to control the false discovery rate (FDR)
dunn_result <- dunn.test(adivA_age$Shannon, adivA_age$Age_Days, method="bh")

dunn_resultob <- dunn.test(adivA_age$Observed, adivA_age$Age_Days, method="bh")
# View the results
print(dunn_result)
#1-10 3.409840e-02 
#1-21 1.914332e-04 
# 10-21 8.097819e-08

print(dunn_resultob)
#1-10 1.338572e-02
#1-21 7.152591e-04 
# 10-21 6.547986e-08

#interaction
adivA_age$interaction <- with(adivA_age, paste(Treatment, Age_Days, sep = "_"))

test_inter <- kruskal.test(Shannon ~ interaction, data = adivA_age)
test_inter

# Post-hoc pairwise comparisons
post_hoc_result <- pairwise.wilcox.test(adivA_age$Shannon, adivA_age$interaction, p.adjust.method = "BH")
print(post_hoc_result)

#test the interaction for Observed 
test_inter_ob <- kruskal.test(Observed ~ interaction, data = adivA_age)
test_inter_ob

# Post-hoc pairwise comparisons for observed 
post_hoc_result_ob <- pairwise.wilcox.test(adivA_age$Observed, adivA_age$interaction, p.adjust.method = "BH")
print(post_hoc_result_ob)

shannonviolAGE <- ggplot(adivA_age,#filter(treatment == "Control"),
                         aes(x = Age_Days, y = Shannon, fill = Age_Days)) +
  geom_boxplot(trim = FALSE, alpha = 0.8)  +
  scale_fill_manual(values = AGE) +
  xlab("") +
  geom_point(position = "jitter", color = "black", size = 3) +
  theme_bw() + 
  labs(title = "A",
       x = "", 
       y = "Shannon's Index",
       shape = "Age_Days",
       fill = "Age_Days") + 
  #scale_x_discrete(name="", labels = c("cloaca" = "Cloaca", "sock"="Sock")) +
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.0)) +
  guides(fill = "none")

shannonviolAGE 


obsviolAGE <- 
  ggplot(adivA_age,#filter(treatment == "Control"),
         aes(x = Age_Days, y = Observed, fill = Age_Days)) +
  geom_boxplot(trim = FALSE, alpha = 0.8)  +
  scale_fill_manual(values = AGE) +
  xlab("") +
  geom_point(position = "jitter", color = "black", size = 3) +
  theme_bw() + 
  labs(title = "B",
       x = "", 
       y = "Observed ASVs",
       shape = "Age_Days",
       fill = "Age_Days") + 
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.0)) +
  guides(fill = "none")

obsviolAGE 

alphdiv_prelim_boxplotAGE_article <- ggarrange(shannonviolAGE, obsviolAGE, ncol = 3, common.legend = TRUE, legend = "right" )

pdf("New_graphs/alphadiv_prelim_boxplotAGE_article.pdf", width = 18, height = 12)
alphdiv_prelim_boxplotAGE_article
dev.off()

ps
# Calculate Shannon diversity for each sample
shannon_diversity <- estimate_richness(ps, measures = "Shannon")

# Calculate the maximum Shannon diversity (log(S))
species_counts <- taxa_sums(ps)
H_max <- log(length(species_counts[species_counts > 0]))

# Calculate Pielou's evenness for each sample
pielou_evenness <- shannon_diversity$Shannon / H_max

# pielou_evenness now contains the evenness index for each sample


###-----------------------------------------###

### subsetting 7 first days of age 

# load data
ps <- readRDS("ps-decontam-filtered-counts.rds")

ps
ps <-  subset_samples(ps, Treatment == "Basal Diet" | Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Essential Oils")

ps <-  subset_samples(ps, Age_Days == "1" | Age_Days == "2" | Age_Days == "3" | Age_Days == "4" | Age_Days == "5" | Age_Days == "6" | Age_Days == "7")

sample_data(ps)

ps

# create data frames ----
adivA_age <- data.frame(
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

ggdensity(adivA_age$Shannon, 
          main = "Density plot of Shannon",
          xlab = "Shannon")

#Q-Q plot: Q-Q plot (or quantile-quantile plot) draws the correlation between a given sample and the normal distribution. A 45-degree reference line is also plotted.
ggqqplot(adivA_age$Shannon)


sh <- ggplot(adivA_age, aes(x=Shannon)) + 
  geom_histogram(fill="blue", color="black") + 
  ggtitle("Histogram of Data") +
  xlab("Shannon") +
  ylab("Frequency") + 
  theme_bw()
sh

obs <- ggplot(adivA_age, aes(x=Observed)) + 
  geom_histogram(fill="blue", color="black") + 
  ggtitle("Histogram of Data") +
  xlab("Observed") +
  ylab("Frequency") 
obs

# Perform Shapiro-Wilk test
shapiro_test_age <- shapiro.test(adivA_age$Shannon)

# Show the results of the Shapiro-Wilk test
print(shapiro_test_age)

#From the output, the p-value < 0.05 implying that the distribution of the data are significantly different from normal distribution. In other words, we can assume the data is not normal.

#since it is not normal lets do the kruskal_test 
#Check age 
kruskal_test <- kruskal.test(Shannon ~ Age_Days, data = adivA_age)
print(kruskal_test)

#age is different 
# Perform Dunn's test for multiple comparisons
# Adjusting the p-values using the Benjamini-Hochberg method to control the false discovery rate (FDR)
dunn_result <- dunn.test(adivA_age$Shannon, adivA_age$Age_Days, method="bh")

# View the results
print(dunn_result)

##Observed analysis 
kruskal_test_ob <- kruskal.test(Observed ~ Age_Days, data = adivA_age)
print(kruskal_test_ob)

#age is different 
# Perform Dunn's test for multiple comparisons
# Adjusting the p-values using the Benjamini-Hochberg method to control the false discovery rate (FDR)
dunn_result_ob <- dunn.test(adivA_age$Observed, adivA_age$Age_Days, method="bh")

# View the results
print(dunn_result_ob)


#interaction
adivA_age$interaction <- with(adivA_age, paste(Treatment, Age_Days, sep = "_"))

test_inter <- kruskal.test(Shannon ~ interaction, data = adivA_age)
test_inter

# Post-hoc pairwise comparisons
post_hoc_result <- pairwise.wilcox.test(adivA_age$Shannon, adivA_age$interaction, p.adjust.method = "BH")
print(post_hoc_result)

#test the interaction for Observed 
test_inter_ob <- kruskal.test(Observed ~ interaction, data = adivA)
test_inter_ob

# Post-hoc pairwise comparisons for observed 
post_hoc_result_ob <- pairwise.wilcox.test(adivA$Observed, adivA$interaction, p.adjust.method = "BH")
print(post_hoc_result_ob)


##Graphs
shannonviolAGE_1_7 <- ggplot(adivA_age,#filter(treatment == "Control"),
                             aes(x = Age_Days, y = Shannon, fill = Age_Days)) +
  geom_boxplot(trim = FALSE, alpha = 0.8)  +
  scale_fill_viridis_d(option="mako") +
  xlab("") +
  geom_point(position = "jitter", color = "black", size = 3) +
  theme_bw() + 
  labs(title = "A",
       x = "", 
       y = "Shannon's Index",
       shape = "Age_Days",
       fill = "Age_Days") + 
  #scale_x_discrete(name="", labels = c("cloaca" = "Cloaca", "sock"="Sock")) +
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.0)) +
  guides(fill = "none")

shannonviolAGE_1_7


obsviolAGE_1_7 <- 
  ggplot(adivA_age,#filter(treatment == "Control"),
         aes(x = Age_Days, y = Observed, fill = Age_Days)) +
  geom_boxplot(trim = FALSE, alpha = 0.8)  +
  scale_fill_viridis_d(option="mako") +
  xlab("") +
  geom_point(position = "jitter", color = "black", size = 3) +
  theme_bw() + 
  labs(title = "B",
       x = "", 
       y = "Observed ASVs",
       shape = "Age_Days",
       fill = "Age_Days") + 
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.0)) +
  guides(fill = "none")

obsviolAGE_1_7

alphdiv_prelim_AGE_1_7_article <- ggarrange(shannonviolAGE_1_7, obsviolAGE_1_7, ncol = 3, common.legend = TRUE, legend = "right" )

pdf("New_graphs/alphadiv_prelim_AGE_1_7_article.pdf", width = 18, height = 12)
alphdiv_prelim_violinAGE_1_7_article
dev.off()


###Extracting the p-values 

# Let's say post_hoc_result is the result of your pairwise.wilcox.test
# To filter out significant interactions, look at the p-value adjustment

# Extracting the p-value matrix from the post-hoc results
p_values <- post_hoc_result$p.value

# Adjusted p-values are already in p_values if you specified p.adjust.method in pairwise.wilcox.test
# Now, filter for significant results
significant_interactions <- p_values < 0.05

# Print significant interactions
print(significant_interactions)

# To get a cleaner list of significant pairs and their p-values
significant_pairs <- which(significant_interactions, arr.ind = TRUE)

# If you want to display the significant pairs with their adjusted p-values
if (length(significant_pairs) > 0) {
  cat("Significant pairwise comparisons:\n")
  apply(significant_pairs, 1, function(idx) {
    cat(rownames(p_values)[idx[1]], "-", colnames(p_values)[idx[2]], 
        "p-value:", p_values[idx[1], idx[2]], "\n")
  })
} else {
  cat("No significant pairwise comparisons found.\n")
}

