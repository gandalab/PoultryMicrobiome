#Contrast Analysis 

library(afex)
library(tidyverse)
library(dplyr)
library(psych)
library(emmeans)
library(graphics)

## R contrast analysis

##load data
performance_data <- read.csv("2023_paper_performance_analysis.csv",
                             header = TRUE, colClasses = c("factor", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                           "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))

#check type of your data 
str(performance_data)
is.na(performance_data)

performance_data <- na.omit(performance_data)

## get the dummy code
c1 <- contrasts(as.factor(performance_data$TRT))

#check the levels 
#with the levels you can see how many contrasts you can have 
#if you have 4 levels k=4, you can have maximum of 3 contrasts k-1 
attributes(performance_data$TRT)

#part 2 - planned contrasts using aov_4() fuction for some reason did not work and I could not figure the error 

#using the contrast in the lm model 
#step 1: define the contrasts using the 7 rules 
ContvsPro <- c(1,-1,0,0)
ContvsEo <- c(1,0,-1,0)
ProvsAbx <- c(0,1,0,-1)

#CtlvsEoPrb <- c(-1, 0.5, 0.5, 0)
#AtxvsEoPrb <- c(0,0.5,0.5, -1)
#ProbEo <- c(0,-1,1,0)

#step 2
# Calculate the inverse of the contrast matrix
#temp was added just to make the matrix square
contr_temp <- rbind(temp=1,
                    ContvsPro,
                    ContvsEo,
                    ProvsAbx)


#contr_temp <- rbind(temp=1,
#CtlvsEoPrb,
# AtxvsEoPrb,
#ProbEo)

#Inverse of the the contr matrix using solve() function
contr <- solve(contr_temp)

#Remove the temp (1st column) since it serves no purpose
contr.notemp <- contr[,-1]

#step 3 - link the contrast to the predictors 
contrasts(performance_data$TRT) <- contr.notemp
attributes(performance_data$TRT)

#step 4 - run the linear model with planned contrasts
m1 <- lm(LIVE_WEIGHT_STARTER ~ TRT, data = performance_data)
summary(m1)

# Extract p-values
p_values1 <- summary(m1)$coefficients[, 4]

# Adjust p-values using Bonferroni correction
adjusted_p_m1 <- p.adjust(p_values1, method = "BH")

#print
print(adjusted_p_m1)
round(adjusted_p_m1, 2)

m2 <- lm(LIVE_WEIGHT_GROWER ~ TRT, data = performance_data)
summary(m2)

# Extract p-values
p_values2 <- summary(m2)$coefficients[, 4]

# Adjust p-values using Bonferroni correction
adjusted_p_m2 <- p.adjust(p_values2, method = "BH")

#print
print(adjusted_p_m2)
round(adjusted_p_m2, 2)

m3 <- lm(LIVE_WEIGHT_TOTAL ~ TRT, data = performance_data)
summary(m3)

# Extract p-values
p_values3 <- summary(m3)$coefficients[, 4]

# Adjust p-values using Bonferroni correction
adjusted_p_m3 <- p.adjust(p_values3, method = "BH")

#print
print(adjusted_p_m3)
round(adjusted_p_m3, 2)

m4 <- lm(FI_STARTER.PEN ~ TRT, data = performance_data)
summary(m4)

# Extract p-values
p_values4 <- summary(m4)$coefficients[, 4]

# Adjust p-values using Bonferroni correction
adjusted_p_m4 <- p.adjust(p_values4, method = "BH")

#print
print(adjusted_p_m4)
round(adjusted_p_m4, 2)

m5 <- lm(FI_GROWER.PEN ~ TRT, data = performance_data)
summary(m5)

# Extract p-values
p_values5 <- summary(m5)$coefficients[, 4]

# Adjust p-values using Bonferroni correction
adjusted_p_m5 <- p.adjust(p_values5, method = "BH")

#print
print(adjusted_p_m5)
round(adjusted_p_m5, 2)

m6 <- lm(FI_TOTAL.PEN ~ TRT, data = performance_data)
summary(m6)

# Extract p-values
p_values6 <- summary(m6)$coefficients[, 4]

# Adjust p-values using Bonferroni correction
adjusted_p_m6 <- p.adjust(p_values6, method = "BH")

#print
print(adjusted_p_m6)
round(adjusted_p_m6, 2)

m7 <- lm(FCR_STARTER.PEN ~ TRT, data = performance_data)
summary(m7)

# Extract p-values
p_values7 <- summary(m7)$coefficients[, 4]

# Adjust p-values using Bonferroni correction
adjusted_p_m7 <- p.adjust(p_values7, method = "BH")

#print
print(adjusted_p_m7)
round(adjusted_p_m7, 2)

m8 <- lm(FCR_GROWER.PEN ~ TRT, data = performance_data)
summary(m8)

# Extract p-values
p_values8 <- summary(m8)$coefficients[, 4]

# Adjust p-values using Bonferroni correction
adjusted_p_m8 <- p.adjust(p_values8, method = "BH")

#print
print(adjusted_p_m8)
round(adjusted_p_m8, 2)

m9 <- lm(FCR_TOTAL.PEN ~ TRT, data = performance_data)
summary(m9)

# Extract p-values
p_values9 <- summary(m9)$coefficients[, 4]

# Adjust p-values using Bonferroni correction
adjusted_p_m9 <- p.adjust(p_values9, method = "BH")

#print
print(adjusted_p_m9)
round(adjusted_p_m9, 2)


##Descriptive statistics 
##this is an overall descriptive 
descriptive_stat <- describe(performance_data1[,2:13], IQR = TRUE)

#descriptive analysis by treatment 
#Basal Diet
Basal <- performance_data1 %>% 
  filter(TRT == "Basal Diet")

descriptive_basal <- describe(Basal[,2:13], IQR = TRUE)

#Probiotic
Prob <- performance_data1 %>% 
  filter(TRT == "Probiotic")

descriptive_prob <- describe(Prob[,2:13], IQR = TRUE)

#Eo 
EO <- performance_data1 %>% 
  filter(TRT == "Essential Oils")

descriptive_eo <- describe(EO[,2:13], IQR = TRUE)

#ATX
Atx <- performance_data1 %>% 
  filter(TRT == "Antibiotic")

descriptive_atx <- describe(Atx[,2:13], IQR = TRUE)

citation(package = "emmeans")

library(gtExtras)
library(svglite)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggExtra)

##Descriptive graph 
performance_data <- na.omit(performance_data) 
performance_data = subset(performance_data, 
                          select = -c(ADG_STARTER.PEN,FI_STARTER.PEN,ADG_GROWER.PEN,FI_GROWER.PEN,ADG_TOTAL.PEN,FI_TOTAL.PEN,EO,PROB,ABX) )
discriptive_graph <- gt_plt_summary(performance_data)

ggsave("discriptive_graph_BD.pdf")

Probiotic <- performance_data %>% 
  filter(TRT == "Probiotic")

PBdescrip <- subset(PB, 
                    select = -c(ADG_STARTER.PEN,FI_STARTER.PEN,ADG_GROWER.PEN,FI_GROWER.PEN,ADG_TOTAL.PEN,FI_TOTAL.PEN,EO,PROB,ABX) )

discriptive_graph_PB <- gt_plt_summary(Probiotic)


EssentialOils <- performance_data %>% 
  filter(TRT == "Essential Oils")

discriptive_graph_EO <- gt_plt_summary(EssentialOils)

Antibiotic <- performance_data %>% 
  filter(TRT == "Antibiotic")

discriptive_graph_ATB <- gt_plt_summary(Antibiotic)

BasalDiet <- performance_data %>% 
  filter(TRT == "Basal Diet")

discriptive_graph_BD <- gt_plt_summary(BasalDiet)


# Save the scatter plot in a variable
p <- ggplot(performance_data, aes(x = TRT, y = LIVE_WEIGHT_STARTER, color = TRT)) +
  geom_point() + theme_bw()

# Marginal histograms by group
ggMarginal(p, size = 2, type = "histogram", 
           groupColour = TRUE,
           groupFill = TRUE)


