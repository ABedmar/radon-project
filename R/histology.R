setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea6")
getwd()

library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(patchwork)

################################################################################
data<-read.csv(file="adeno_histology.CSV",sep=",",header = T)

df<-data[,-c(7:11)]
df<-melt(df, id.vars = unique(c("group")))

a<-ggplot(data=df,aes(x=factor(group), y=value, fill=variable, label=value))+
  geom_bar(stat="identity")+
  geom_text(data = df %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 20, vjust = 1.3, hjust=1),legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("Adenocarcinoma")+
  xlab(element_blank())
a

df2<-data[,-c(2,3,4,5,6,9,10,11)]
df2<-melt(df2, id.vars = unique(c("group")))

b<-ggplot(data=df2,aes(x=factor(group), y=value, fill=variable, label=value))+
  geom_bar(stat="identity")+
  geom_text(data = df2 %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 20, vjust = 1.3, hjust=1),legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab(element_blank())+
  xlab(element_blank())
b

df3<-data[,-c(2,3,4,5,6,7,8)]
df3<-melt(df3, id.vars = unique(c("group")))

c<-ggplot(data=df3,aes(x=factor(group), y=value, fill=variable, label=value))+
  geom_bar(stat="identity")+
  geom_text(data = df3 %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 20, vjust = 1.3, hjust=1),legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab(element_blank())+
  xlab(element_blank())+
  scale_fill_manual(
    values = c("#f8766d", "#619cff", "#b3eac3")
  )
c

################################################################################
data<-read.csv(file="squamous_histology.CSV",sep=",")

df<-data[,-c(7:11)]
df<-melt(df, id.vars = unique(c("group")))

d<-ggplot(data=df,aes(x=factor(group), y=value, fill=variable, label=value))+
  geom_bar(stat="identity")+
  geom_text(data = df %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 20, vjust = 1.3, hjust=1),legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("Squamous")+
  xlab(element_blank())


df2<-data[,-c(2,3,4,5,6,9,10,11)]
df2<-melt(df2, id.vars = unique(c("group")))
e<-ggplot(data=df2,aes(x=factor(group), y=value, fill=variable, label=value))+
  geom_bar(stat="identity")+
  geom_text(data = df2 %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 20, vjust = 1.3, hjust=1),legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab(element_blank())+
  xlab(element_blank())

df3<-data[,-c(2,3,4,5,6,7,8)]
df3<-melt(df3, id.vars = unique(c("group")))

f<-ggplot(data=df3,aes(x=factor(group), y=value, fill=variable, label=value))+
  geom_bar(stat="identity")+
  geom_text(data = df3 %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 20, vjust = 1.3, hjust=1),legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab(element_blank())+
  xlab(element_blank())+
  scale_fill_manual(
    values = c("#f8766d", "#619cff", "#b3eac3")
  )
f

grid.arrange(a,b,c,d,e,f,nrow = 2)

