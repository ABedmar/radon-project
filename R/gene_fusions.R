setwd('/Users/alexb/OneDrive/Escritorio/IDIBAPS')

data <- read.csv(file = "gene_fusion_12June2022.csv", sep = ";", header = T)
head(data)

##################################
#T-Test

library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(reshape2)

pwc <- data %>%
  pairwise_t_test(Total_gene_fusions ~ Group, p.adjust.method = "bonferroni")
pwc

pwc_inter <- data %>%
  pairwise_t_test(Inter.Chromosomal ~ Group, p.adjust.method = "bonferroni")
pwc_inter

pwc_intra <- data %>%
  pairwise_t_test(Intra.Chromosomal ~ Group, p.adjust.method = "bonferroni")
pwc_intra

res.aov <- data %>% anova_test(Total_gene_fusions ~ Group)
res.aov


#################################
#Plots

ggplot(data, aes(x=Sample , y=Total_gene_fusions, fill=Group))+
  geom_col()+
  ggtitle("")+
  ylab("")+
  xlab("")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


data_grouped <- read.csv(file = "gene_fusion_groups_12June2022.csv", sep = ";", header = T)

ggplot(data_grouped, aes(x=Group , y=Sum, fill=Group))+
  geom_col()+
  ggtitle("Total gene fusions")+
  ylab("")+
  xlab("")+
  theme_minimal()+
  theme(legend.title=element_blank())

ggplot(data_grouped, aes(x=Group , y=Mean, fill=Group))+
  geom_col()+
  ggtitle("Mean gene fusions per group")+
  ylab("")+
  xlab("")+
  theme_minimal()+
  theme(legend.title=element_blank())


df <-data_grouped[,-c(2,3)]
df1 <- data.frame(df)
df2 <- melt(df1, id.vars='Group')
head(df2)

ggplot(df2, aes(x=Group, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge')+
  ggtitle("Mean gene fusions")+
  ylab("")+
  xlab("")+
  theme_minimal()+
  theme(legend.title=element_blank())+
  scale_fill_discrete(labels=c('Intra Chromosomal', 'Inter Chromosomal'))



