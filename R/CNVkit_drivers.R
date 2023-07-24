  setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/CNVkit/drivers/")
getwd()
library(maftools)
library(RColorBrewer)
library(stringr)
library(tibble)
library(mclust)
library(ggplot2)
library(data.table)
library(tidyr)
library(plotfunctions)
library("gplots")
library(gridExtra)
library(dplyr)
library("viridis")

data<-read.table(file="108_88.TXT", header = TRUE)

head(data)

medians<-data %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)


ggplot()+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  geom_jitter(data=data, aes(x=gene, y=log2, color=gene), show.legend = FALSE)+
  geom_line(data=medians, aes(x=gene, y=log2, group=1), show.legend = FALSE)+
  geom_hline(yintercept = median(data$log2),color = "black", linetype="dashed")+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  xlab(element_blank())+
  ggtitle("108_88")+
  ylim(-18.1368, 3.88749)



sample <- c("108_88.txt")
name<-tools::file_path_sans_ext(sample)
name
name<-paste0(name,".svg")
name

####################################################
#sample <- c("108_88.txt","118_86.txt","119_90.txt","166_86.txt","219_86.txt","269_87.txt","278_90.txt","292_87.txt","300_88.txt","309_90.txt","340_87.txt","350_88.txt","372_89.txt","374_90.txt","389_88.txt","419_87.txt","442_87.txt","443_88.txt","457_89.txt","520_89.txt","585_88.txt","603_87.txt","616_90.txt","631_88.txt","631_90.txt","635_89.txt","639_90.txt")

sample <- c("108_88.txt","118_86.txt","119_90.txt","166_86.txt","219_86.txt","269_87.txt","278_90.txt","292_87.txt","300_88.txt","309_90.txt","340_87.txt","350_88.txt","372_89.txt","374_90.txt","389_88.txt","419_87.txt","442_87.txt","443_88.txt","457_89.txt","520_89.txt","585_88.txt","603_87.txt","616_90.txt","631_88.txt","631_90.txt","635_89.txt","639_90.txt")

for(i in 1:length(sample)) {
  data<-read.table(file=sample[i], header = TRUE)
  name<-tools::file_path_sans_ext(sample[i])
  
  head(data)
  
  medians<-data %>% 
    group_by(gene) %>%  
    summarise_all(.funs = median)
  
  filename<-paste0(name,".svg")

  p<-ggplot()+
    geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
    scale_shape_identity()+
    scale_size_identity()+
    geom_jitter(data=data, aes(x=gene, y=log2, color=gene), show.legend = FALSE)+
    geom_line(data=medians, aes(x=gene, y=log2, group=1), show.legend = FALSE)+
    geom_hline(yintercept = median(data$log2),color = "black", linetype="dashed")+
    geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
    scale_shape_identity()+
    scale_size_identity()+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
    xlab(element_blank())+
    ggtitle(name)+
    ylim(-18.1368, 3.88749)
  
  
  ggsave(file=filename, plot=p)
  
}

####################################### SCATTELPLOT MERGED

data<-read.table(file="all.txt", header = TRUE)

data<-data[order(data$gene),]

medians<-data %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)


p<-ggplot()+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  geom_jitter(data=data, aes(x=gene, y=log2, color=gene), show.legend = FALSE)+
  geom_line(data=medians, aes(x=gene, y=log2, group=1), show.legend = FALSE)+
  geom_hline(yintercept = median(data$log2),color = "black", linetype="dashed")+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  xlab(element_blank())+
  ggtitle("All drivers")+
  ylim(-18.1368, 3.88749)

p

ggsave(file="All_drivers.svg", plot=p)




####################################### SCATTELPLOT MERGED BY GROUPS
# 2 GROUPS

data<-read.table(file="2GROUPS/HD/HD.txt", header = TRUE)
medians<-data %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)
hd<-ggplot()+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  geom_jitter(data=data, aes(x=gene, y=log2, color=gene), show.legend = FALSE)+
  geom_line(data=medians, aes(x=gene, y=log2, group=1), show.legend = FALSE)+
  geom_hline(yintercept = median(data$log2),color = "black", linetype="dashed")+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  xlab(element_blank())+
  ggtitle("High_dose    N=14")+
  ylim(-18.1368, 3.88749)

data<-read.table(file="2GROUPS/LD/LD.txt", header = TRUE)
medians<-data %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)
ld<-ggplot()+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  geom_jitter(data=data, aes(x=gene, y=log2, color=gene), show.legend = FALSE)+
  geom_line(data=medians, aes(x=gene, y=log2, group=1), show.legend = FALSE)+
  geom_hline(yintercept = median(data$log2),color = "black", linetype="dashed")+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  xlab(element_blank())+
  ggtitle("Low_dose    N=13")+
  ylim(-18.1368, 3.88749)



grid.arrange(hd,ld,nrow=1)

ggsave(file="2GROUPS/high_vs_low_drivers_CNA.svg", plot=grid.arrange(hd,ld,nrow=1))




# 3 GROUPS

data<-read.table(file="3GROUPS/HD/HD.txt", header = TRUE)
medians<-data %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)
hd<-ggplot()+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  geom_jitter(data=data, aes(x=gene, y=log2, color=gene), show.legend = FALSE)+
  geom_line(data=medians, aes(x=gene, y=log2, group=1), show.legend = FALSE)+
  geom_hline(yintercept = median(data$log2),color = "black", linetype="dashed")+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1.9, hjust=1))+
  xlab(element_blank())+
  ggtitle("High_dose    N=14")+
  ylim(-18.1368, 3.88749)

data<-read.table(file="3GROUPS/ID/ID.txt", header = TRUE)
medians<-data %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)
id<-ggplot()+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  geom_jitter(data=data, aes(x=gene, y=log2, color=gene), show.legend = FALSE)+
  geom_line(data=medians, aes(x=gene, y=log2, group=1), show.legend = FALSE)+
  geom_hline(yintercept = median(data$log2),color = "black", linetype="dashed")+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1.9, hjust=1))+
  xlab(element_blank())+
  ggtitle("Intermediate_dose    N=9")+
  ylim(-18.1368, 3.88749)

data<-read.table(file="3GROUPS/LD/LD.txt", header = TRUE)
medians<-data %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)
ld<-ggplot()+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  geom_jitter(data=data, aes(x=gene, y=log2, color=gene), show.legend = FALSE)+
  geom_line(data=medians, aes(x=gene, y=log2, group=1), show.legend = FALSE)+
  geom_hline(yintercept = median(data$log2),color = "black", linetype="dashed")+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1.9, hjust=1))+
  xlab(element_blank())+
  ggtitle("Low_dose    N=4")+
  ylim(-18.1368, 3.88749)



grid.arrange(hd,id,ld,nrow=1)

ggsave(file="3GROUPS/high_vs_intermediate_vs_low_drivers_CNA.svg", plot=grid.arrange(hd,id,ld,nrow=1))

################################################# HEATMAP MEDIAN

sample <- c("119_90.txt","166_86.txt","219_86.txt","269_87.txt","278_90.txt","292_87.txt","300_88.txt","309_90.txt","340_87.txt","350_88.txt","372_89.txt","374_90.txt","389_88.txt","419_87.txt","442_87.txt","443_88.txt","457_89.txt","520_89.txt","585_88.txt","603_87.txt","616_90.txt","631_88.txt","631_90.txt","635_89.txt","639_90.txt")

d1<-read.table(file="108_88.txt", header = TRUE)
d2<-read.table(file="118_86.txt", header = TRUE)

m1<-d1 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)

m2<-d2 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)

df <- merge(m1,m2,by="gene",all.x = TRUE)
df

for(i in 1:length(sample)) {
  
  data<-read.table(file=sample[i], header = TRUE)
  
  medians<-data %>% 
    group_by(gene) %>%  
    summarise_all(.funs = median)
  
  df <- merge(df,medians,by="gene",all.x =TRUE)
}

df
colnames(df) <- c("gene","108_88","118_86","119_90","166_86","219_86","269_87","278_90","292_87","300_88","309_90","340_87","350_88","372_89","374_90","389_88","419_87","442_87","443_88","457_89","520_89","585_88","603_87","616_90","631_88","631_90","635_89","639_90")
df

rownames(df) <- df[,1]
df <- df[,-1]

data<-as.matrix(df)
heatmap.2(data, scale = "none", col = bluered(100), trace = "none", density.info = "none", main = "Median log2 CNA",  lhei=c(2,10), lwid=c(2,4), keysize=0.75, key.title = "log2",Colv=FALSE)


############################HEATMAP MEAN

sample <- c("119_90.txt","166_86.txt","219_86.txt","269_87.txt","278_90.txt","292_87.txt","300_88.txt","309_90.txt","340_87.txt","350_88.txt","372_89.txt","374_90.txt","389_88.txt","419_87.txt","442_87.txt","443_88.txt","457_89.txt","520_89.txt","585_88.txt","603_87.txt","616_90.txt","631_88.txt","631_90.txt","635_89.txt","639_90.txt")

d1<-read.table(file="108_88.txt", header = TRUE)
d2<-read.table(file="118_86.txt", header = TRUE)

m1<-d1 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = mean)

m2<-d2 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = mean)

df <- merge(m1,m2,by="gene",all.x =TRUE)
df

for(i in 1:length(sample)) {
  
  data<-read.table(file=sample[i], header = TRUE)
  
  medians<-data %>% 
    group_by(gene) %>%  
    summarise_all(.funs = mean)
  
  df <- merge(df,medians,by="gene",all.x =TRUE)
}

df
colnames(df) <- c("gene","108_88","118_86","119_90","166_86","219_86","269_87","278_90","292_87","300_88","309_90","340_87","350_88","372_89","374_90","389_88","419_87","442_87","443_88","457_89","520_89","585_88","603_87","616_90","631_88","631_90","635_89","639_90")
df

rownames(df) <- df[,1]
df <- df[,-1]


data<-as.matrix(df)
heatmap.2(data, scale = "none", col = bluered(100), trace = "none", density.info = "none", main = "Mean log2 CNA",  lhei=c(2,10), lwid=c(2,4), keysize=0.75, key.title = "log2",Colv=FALSE)




#################################################################### HEATMAP BY GROUPS

###################### MEDIAN
# TWO GROUPS

sample <- c("166_86.txt","219_86.txt","350_88.txt","372_89.txt","374_90.txt","442_87.txt","443_88.txt","457_89.txt","585_88.txt","616_90.txt","631_88.txt","639_90.txt","118_86.txt","269_87.txt","278_90.txt","292_87.txt","300_88.txt","309_90.txt","340_87.txt","389_88.txt","419_87.txt","520_89.txt","603_87.txt","631_90.txt","635_89.txt")


d1<-read.table(file="108_88.txt", header = TRUE)
d2<-read.table(file="119_90.txt", header = TRUE)

m1<-d1 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)

m2<-d2 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)

df <- merge(m1,m2,by="gene",all.x = TRUE)
df

for(i in 1:length(sample)) {
  
  data<-read.table(file=sample[i], header = TRUE)
  
  medians<-data %>% 
    group_by(gene) %>%  
    summarise_all(.funs = median)
  
  df <- merge(df,medians,by="gene",all.x =TRUE)
}

df
colnames(df) <- c("gene","108_88","119_90","166_86","219_86","350_88","372_89","374_90","442_87","443_88","457_89","585_88","616_90","631_88","639_90","118_86","269_87","278_90","292_87","300_88","309_90","340_87","389_88","419_87","520_89","603_87","631_90","635_89")
df

rownames(df) <- df[,1]
df <- df[,-1]

data<-as.matrix(df)
heatmap.2(data, scale = "none", col = bluered(100), trace = "none", density.info = "none", main = "Median log2 CNA",sub="Dose: high vs low",  lhei=c(2,10), lwid=c(2,4), keysize=0.75, key.title = "log2",Colv=FALSE)



############################### MEAN

# TWO GROUPS
sample <- c("166_86.txt","219_86.txt","350_88.txt","372_89.txt","374_90.txt","442_87.txt","443_88.txt","457_89.txt","585_88.txt","616_90.txt","631_88.txt","639_90.txt","118_86.txt","269_87.txt","278_90.txt","292_87.txt","300_88.txt","309_90.txt","340_87.txt","389_88.txt","419_87.txt","520_89.txt","603_87.txt","631_90.txt","635_89.txt")

d1<-read.table(file="108_88.txt", header = TRUE)
d2<-read.table(file="119_90.txt", header = TRUE)

m1<-d1 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = mean)

m2<-d2 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = mean)

df <- merge(m1,m2,by="gene",all.x =TRUE)
df

for(i in 1:length(sample)) {
  
  data<-read.table(file=sample[i], header = TRUE)
  
  medians<-data %>% 
    group_by(gene) %>%  
    summarise_all(.funs = mean)
  
  df <- merge(df,medians,by="gene",all.x =TRUE)
}

df
colnames(df) <- c("gene","108_88","119_90","166_86","219_86","350_88","372_89","374_90","442_87","443_88","457_89","585_88","616_90","631_88","639_90","118_86","269_87","278_90","292_87","300_88","309_90","340_87","389_88","419_87","520_89","603_87","631_90","635_89")
df

rownames(df) <- df[,1]
df <- df[,-1]


data<-as.matrix(df)
heatmap.2(data, scale = "none", col = bluered(100), trace = "none", density.info = "none", main = "Mean log2 CNA",sub="Dose: high vs low",  lhei=c(2,10), lwid=c(2,4), keysize=0.75, key.title = "log2",Colv=FALSE)




















three_groups <- c("108_88.txt","119_90.txt","166_86.txt","219_86.txt","350_88.txt","372_89.txt","374_90.txt","442_87.txt","443_88.txt","457_89.txt","585_88.txt","616_90.txt","631_88.txt","639_90.txt","269_87.txt","278_90.txt","292_87.txt","300_88.txt","309_90.txt","340_87.txt","603_87.txt","631_90.txt","635_89.txt","118_86.txt","389_88.txt","419_87.txt","520_89.txt")

sample <- c("108_88.txt","119_90.txt","166_86.txt","219_86.txt","350_88.txt","372_89.txt","374_90.txt","442_87.txt","443_88.txt","457_89.txt","585_88.txt","616_90.txt","631_88.txt","639_90.txt","269_87.txt","278_90.txt","292_87.txt","300_88.txt","309_90.txt","340_87.txt","603_87.txt","631_90.txt","635_89.txt","118_86.txt","389_88.txt","419_87.txt","520_89.txt")






###################### MEDIAN
# THREE GROUPS

sample <- c("166_86.txt","219_86.txt","350_88.txt","372_89.txt","374_90.txt","442_87.txt","443_88.txt","457_89.txt","585_88.txt","616_90.txt","631_88.txt","639_90.txt","269_87.txt","278_90.txt","292_87.txt","300_88.txt","309_90.txt","340_87.txt","603_87.txt","631_90.txt","635_89.txt","118_86.txt","389_88.txt","419_87.txt","520_89.txt")


d1<-read.table(file="108_88.txt", header = TRUE)
d2<-read.table(file="119_90.txt", header = TRUE)

m1<-d1 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)

m2<-d2 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)

df <- merge(m1,m2,by="gene",all.x = TRUE)
df

for(i in 1:length(sample)) {
  
  data<-read.table(file=sample[i], header = TRUE)
  
  medians<-data %>% 
    group_by(gene) %>%  
    summarise_all(.funs = median)
  
  df <- merge(df,medians,by="gene",all.x =TRUE)
}

df
colnames(df) <- c("gene","108_88","119_90","166_86","219_86","350_88","372_89","374_90","442_87","443_88","457_89","585_88","616_90","631_88","639_90","269_87","278_90","292_87","300_88","309_90","340_87","603_87","631_90","635_89","118_86","389_88","419_87","520_89")
df

rownames(df) <- df[,1]
df <- df[,-1]

data<-as.matrix(df)
heatmap.2(data, scale = "none", col = bluered(100), trace = "none", density.info = "none", main = "Median log2 CNA",sub="Dose: high vs intermediate vs low",  lhei=c(2,10), lwid=c(2,4), keysize=0.75, key.title = "log2",Colv=FALSE)



############################### MEAN

# THREE GROUPS
sample <- c("166_86.txt","219_86.txt","350_88.txt","372_89.txt","374_90.txt","442_87.txt","443_88.txt","457_89.txt","585_88.txt","616_90.txt","631_88.txt","639_90.txt","118_86.txt","269_87.txt","278_90.txt","292_87.txt","300_88.txt","309_90.txt","340_87.txt","389_88.txt","419_87.txt","520_89.txt","603_87.txt","631_90.txt","635_89.txt")

d1<-read.table(file="108_88.txt", header = TRUE)
d2<-read.table(file="119_90.txt", header = TRUE)

m1<-d1 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = mean)

m2<-d2 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = mean)

df <- merge(m1,m2,by="gene",all.x =TRUE)
df

for(i in 1:length(sample)) {
  
  data<-read.table(file=sample[i], header = TRUE)
  
  medians<-data %>% 
    group_by(gene) %>%  
    summarise_all(.funs = mean)
  
  df <- merge(df,medians,by="gene",all.x =TRUE)
}

df
colnames(df) <- c("gene","108_88","119_90","166_86","219_86","350_88","372_89","374_90","442_87","443_88","457_89","585_88","616_90","631_88","639_90","269_87","278_90","292_87","300_88","309_90","340_87","603_87","631_90","635_89","118_86","389_88","419_87","520_89")
df

rownames(df) <- df[,1]
df <- df[,-1]


data<-as.matrix(df)
heatmap.2(data, scale = "none", col = bluered(100), trace = "none", density.info = "none", main = "Mean log2 CNA",sub="Dose: high vs intermediate vs low",  lhei=c(2,10), lwid=c(2,4), keysize=0.75, key.title = "log2",Colv=FALSE)
















sample <- c("166_86.txt","219_86.txt","350_88.txt","372_89.txt","374_90.txt","442_87.txt","443_88.txt","457_89.txt","585_88.txt","616_90.txt","631_88.txt","639_90.txt","118_86.txt","269_87.txt","278_90.txt","292_87.txt","300_88.txt","309_90.txt","340_87.txt","389_88.txt","419_87.txt","520_89.txt","603_87.txt","631_90.txt","635_89.txt")

d1<-read.table(file="108_88.txt", header = TRUE)
d2<-read.table(file="119_90.txt", header = TRUE)


df <- merge(d1,d2,by="gene",all.x =TRUE)
df



for(i in 1:length(sample)) {
  
  data<-read.table(file=sample[i], header = TRUE)
  
  df <- merge(df,data,by="gene",all.x =TRUE)
}

head(df)
rownames(df) <- df[,1]
df <- df[,-1]
data<-as.matrix(df)
min(data, na.rm = TRUE)
max(data, na.rm = TRUE)









































################################################# HEATMAP SD
sample <- c("119_90.txt","166_86.txt","219_86.txt","269_87.txt","278_90.txt","292_87.txt","300_88.txt","309_90.txt","340_87.txt","350_88.txt","372_89.txt","374_90.txt","389_88.txt","419_87.txt","442_87.txt","443_88.txt","457_89.txt","520_89.txt","585_88.txt","603_87.txt","616_90.txt","631_88.txt","631_90.txt","635_89.txt","639_90.txt")

d1<-read.table(file="108_88.txt", header = TRUE)
d2<-read.table(file="118_86.txt", header = TRUE)

m1<-d1 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = sd, na.rm = TRUE)

m2<-d2 %>% 
  group_by(gene) %>%  
  summarise_all(.funs = sd)

df <- merge(m1,m2,by="gene",all.x = TRUE)
df

for(i in 1:length(sample)) {
  
  data<-read.table(file=sample[i], header = TRUE)
  
  medians<-data %>% 
    group_by(gene) %>%  
    summarise_all(.funs = sd)
  
  df <- merge(df,medians,by="gene",all.x =TRUE)
}

df
colnames(df) <- c("gene","108_88","118_86","119_90","166_86","219_86","269_87","278_90","292_87","300_88","309_90","340_87","350_88","372_89","374_90","389_88","419_87","442_87","443_88","457_89","520_89","585_88","603_87","616_90","631_88","631_90","635_89","639_90")
df

rownames(df) <- df[,1]
df <- df[,-1]



df2<-filter(df, rowSums(is.na(df)) != ncol(df))

data<-as.matrix(df)
data[is.na(data)] <- 0

heatmap.2(data, scale = "none", col = bluered(100), trace = "none", density.info = "none", main = "SD log2 CNA",  lhei=c(2,10), lwid=c(2,4), keysize=0.75, key.title = "log2",Colv=FALSE)








############################################
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cnv/CNVkit/")


data<-read.table(file="all_dirvers.CNR", header = TRUE)

data<-data[,-1]
head(data)

max(data$log2)
min(data$log2)

medians<-data %>% 
  group_by(gene) %>%  
  summarise_all(.funs = median)


ggplot()+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  geom_jitter(data=data, aes(x=gene, y=log2, color=gene), show.legend = FALSE)+
  geom_line(data=medians, aes(x=gene, y=log2, group=1), show.legend = FALSE)+
  geom_hline(yintercept = median(data$log2),color = "black", linetype="dashed")+
  geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
  scale_shape_identity()+
  scale_size_identity()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  xlab(element_blank())+
  ggtitle("Radon rats all samples")+
  ylim(-24.8694, 21.0274)


####################################################
dir <- "C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cnv/CNVkit/"

all_files <- list.files(dir)

# Filter the files to only include those ending in ".txt"
txt_files <- all_files[grep("\\.txt$", all_files)]


for(i in 1:length(txt_files)) {
  data<-read.table(file=txt_files[i], header = TRUE)
  name<-tools::file_path_sans_ext(txt_files[i])
  
  head(data)
  
  medians<-data %>% 
    group_by(gene) %>%  
    summarise_all(.funs = median)
  
  filename<-paste0(name,".svg")
  
  p<-ggplot()+
    geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
    scale_shape_identity()+
    scale_size_identity()+
    geom_jitter(data=data, aes(x=gene, y=log2, color=gene), show.legend = FALSE)+
    geom_line(data=medians, aes(x=gene, y=log2, group=1), show.legend = FALSE)+
    geom_hline(yintercept = median(data$log2),color = "black", linetype="dashed")+
    geom_point(data=medians, aes(x=reorder(gene, -log2), y=log2, shape=18, size=3), show.legend = FALSE)+
    scale_shape_identity()+
    scale_size_identity()+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
    xlab(element_blank())+
    ggtitle(name)+
    ylim(-24.8694, 21.0274)
  
  
  ggsave(file=filename, plot=p)
  
}
