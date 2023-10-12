#Uniformity graphs 
# packages 
library(ggplot2)
library(magrittr)
library(dplyr)
library(ggpubr)


# set directory 
setwd("/Users/anafonseca/Chicken_Microbiome_data/broiler_microbiome_groups/fastpfiles/asv-fastp/phyloseq_object/")

# Get the data 
#unidata <- readxl::read_excel("/Users/anafonseca/OneDrive - The Pennsylvania State University/Documents/Master_projects/Poultry_Science_Manuscript/BW_UNIFORMITY.xlsx")

#unidata

df1 <- data.frame(Treatment=c("Basal Diet", "Probiotic", "Essential Oils", "Antibiotic"),
                 Uniformity=c(81.2, 76.2, 86.2, 83.8),
                 Weight=c(1, 1, 1, 1))

head(df1)

df2 <- data.frame(Treatment=c("Basal Diet", "Probiotic", "Essential Oils", "Antibiotic"),
                  Uniformity=c(55.1, 67.5, 53.2, 53.8),
                  Weight=c(2, 2, 2, 2))

head(df2)

df3 <- data.frame(Treatment=c("Basal Diet", "Probiotic", "Essential Oils", "Antibiotic"),
                  Uniformity=c(54.5, 68.8, 65.8, 67.9),
                  Weight=c(3, 3, 3, 3))

head(df3)

# plots

df1$"Treatment" <- factor(df1$"Treatment", 
                            levels = c("Antibiotic", "Essential Oils", "Probiotic", "Basal Diet"))


df2$"Treatment" <- factor(df2$"Treatment", 
                          levels = c("Antibiotic", "Essential Oils", "Probiotic", "Basal Diet"))

df3$"Treatment" <- factor(df3$"Treatment", 
                          levels = c("Antibiotic", "Essential Oils", "Probiotic", "Basal Diet"))
#Weight1

p1 <-ggplot(df1, aes(x=Treatment, y=Uniformity, fill=Treatment)) +
  geom_bar(stat="identity",  width=0.7) + 
  scale_y_continuous(limits = c(0, 90)) +
  coord_flip() +  theme_minimal()
p1
ggsave(filename = "bar-plot-treatmentuniformity1.pdf", dpi = 600)

#Weight2
p2 <-ggplot(df2, aes(x=Treatment, y=Uniformity, fill=Treatment)) +
  geom_bar(stat="identity",  width=0.7) +  
  scale_y_continuous(limits = c(0, 90)) +
  coord_flip() +  theme_minimal()
p2
ggsave(filename = "bar-plot-treatmentuniformity2.pdf", dpi = 600)

#Weight3
p3 <-ggplot(df3, aes(x=Treatment, y=Uniformity, fill=Treatment)) +
  geom_bar(stat="identity",  width=0.7) +  
  scale_y_continuous(limits = c(0, 90)) +
  coord_flip() + theme_minimal()
p3
ggsave(filename = "bar-plot-treatmentuniformity3.pdf", dpi = 600)

ggarrange(p1, p2, p3, ncol = 2, nrow = 2, common.legend = TRUE)

ggsave(filename = "bar-plot-treatmentuniformityall.pdf", dpi = 600)

