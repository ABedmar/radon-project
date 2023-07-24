setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/TNM")
getwd()

library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(patchwork)


data<-read.csv(file="T.CSV",sep=",")

df<-data[c(1:5),]
df<-melt(df, id.vars = unique(c("T")))
a<-ggplot(data=df,aes(x=factor(T), y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = df %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("T")

df2<-data[c(6,7),]
df2<-melt(df2, id.vars = unique(c("T")))
b<-ggplot(data=df2,aes(x=factor(T), y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = df2 %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("T")

df3<-data[c(8,9,10),]
df3<-melt(df3, id.vars = unique(c("T")))

c<-ggplot(data=df3,aes(x=factor(T), y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = df3 %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("T")

data<-read.csv(file="N.CSV",sep=",")

df<-data[c(1:5),]
df<-melt(df, id.vars = unique(c("N")))
d<-ggplot(data=df,aes(x=factor(N), y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = df %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("N")+
  scale_fill_manual(values = c("#b3b3b3","#e5d944","black","green","yellow"))

df2<-data[c(6,7),]
df2<-melt(df2, id.vars = unique(c("N")))
e<-ggplot(data=df2,aes(x=factor(N), y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = df2 %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("N")+
  scale_fill_manual(values = c("#b3b3b3","#e5d944","black","green","yellow"))

df3<-data[c(8,9,10),]
df3<-melt(df3, id.vars = unique(c("N")))

f<-ggplot(data=df3,aes(x=factor(N), y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = df3 %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("N")+
  scale_fill_manual(values = c("#b3b3b3","#e5d944","black","green","yellow"))

data<-read.csv(file="M.CSV",sep=",")

df<-data[c(1:5),]
df<-melt(df, id.vars = unique(c("M")))
g<-ggplot(data=df,aes(x=factor(M), y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = df %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("M")+
  scale_fill_manual(values = c("#b3b3b3","#e5d944","black","green","yellow"))
g
df2<-data[c(6,7),]
df2<-melt(df2, id.vars = unique(c("M")))
h<-ggplot(data=df2,aes(x=factor(M), y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = df2 %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("M")+
  scale_fill_manual(values = c("#b3b3b3","#e5d944","black","green","yellow"))

df3<-data[c(8,9,10),]
df3<-melt(df3, id.vars = unique(c("M")))

i<-ggplot(data=df3,aes(x=factor(M), y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = df3 %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("M")+
  scale_fill_manual(values = c("#b3b3b3","#e5d944","black","green","yellow"))


grid.arrange(a,b,c,d,e,f,g,h,i,ncol=3)



######################################### T GROUPED


data<-read.csv(file="T_grouped.CSV",sep=",")

df<-data[c(1:5),]
df<-melt(df, id.vars = unique(c("T")))
j<-ggplot(data=df,aes(x=factor(T), y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = df %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("T")+
  scale_fill_manual(values = c("black","#3bba5c","#3b80ba"))
j  

df2<-data[c(6,7),]
df2<-melt(df2, id.vars = unique(c("T")))
k<-ggplot(data=df2,aes(x=factor(T), y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = df2 %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("T")+
  scale_fill_manual(values = c("black","#3bba5c","#3b80ba"))

df3<-data[c(8,9,10),]
df3<-melt(df3, id.vars = unique(c("T")))

l<-ggplot(data=df3,aes(x=factor(T), y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = df3 %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("T")+
  scale_fill_manual(values = c("black","#3bba5c","#3b80ba"))




grid.arrange(j,d,g, ncol = 3)
grid.arrange(k,e,h, ncol = 3)
grid.arrange(l,f,i, ncol = 3)



setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea7")




tnm<-read.table(file="TNM.txt",header=T)

tnm<-melt(tnm, id.vars = unique(c("tnm")))

tnm1<-ggplot(data=tnm,aes(x=tnm, y=value, fill=variable, label = value))+
  geom_bar(position = "stack",stat="identity", position = "stack")+
  geom_text(data = tnm %>% dplyr::filter(value > 0), position = position_stack(vjust = 0.5))+
  theme_minimal()+
  ylab(element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("TNM")+
  guides(fill = guide_legend(title = "Duration"))


tnm2<-ggplot(data=tnm,aes(x=tnm, y=value, fill=variable, label = value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = tnm %>% dplyr::filter(value > 0), position = position_stack(vjust = 0.5))+
  theme_minimal()+
  ylab(element_blank())+
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())+
  xlab("TNM")+
  guides(fill = guide_legend(title = "Duration"))+
  scale_y_continuous(breaks = seq(0, 14, by = 1))


grid.arrange(tnm1,tnm2,ncol=2)



