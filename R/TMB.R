getwd()
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/TMB")
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/TMB")

library(ggplot2)
library(plotly)
library(gridExtra)
library(forcats)
library(ggpubr)


data<-read.csv(file="TMB-radon-rats.CSV",sep=",",header=T)
rownames(data)<-data$sample

data<-data[!(data$group=="Control"),]

head(data)

#Idea 1
idea1<-ggplot(data=data)+
  geom_col(aes(x=reorder(Tumor_Sample_Barcode, -total_perMB), y=total_perMB))+
  theme_minimal()+
  geom_hline(yintercept=median(data$total_perMB), color="red", linetype="dashed")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),)+
  theme(legend.position = "bottom")+
  xlab(element_blank())+
  ylab("Total (per MB)")
idea1


#Idea 2.1
#data$group<-factor(data$group, levels = unique(data$group))
#data$total_perMB<-factor(data$total_perMB, levels = unique(data$total_perMB))
data$group <- factor(data$group, levels = c("D1", "D3", "D6", "D12", "Pr"))
data$Tumor_Sample_Barcode <- factor(data$Tumor_Sample_Barcode, levels = c("D1-104","D1-132","D1-134","D1-162","D1-171","D1-189","D1-195","D1-207","D1-214","D1-36","D1-81","D1-98","D3-168","D3-203","D3-212","D3-232","D3-233","D3-71","D3-92","D6-132","D6-152","D6-179","D6-204","D12-145","D12-178","D12-184","D12-185-1","D12-185-2","D12-229","D12-231","D12-49","D12-51","Pr109","Pr140","Pr146","Pr214"))

idea2.1.a<-ggplot(data = data) +
  geom_col(aes(x = Tumor_Sample_Barcode, y = total_perMB, fill = group)) +
  theme_minimal() +
  geom_hline(yintercept = median(as.numeric(data$total_perMB)), color = "red", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.title = element_blank()) +
  xlab(element_blank()) +
  ylab("Total (per MB)") +
  labs(title = "A",
       subtitle = "Tumor Mutational Burden in Rn Rats by Sample")+
  scale_fill_viridis_d()+
  theme(plot.title = element_text(size = 12, face = "bold"))


means<-aggregate(data$total_perMB, list(data$group), FUN=mean)
means[2]<-round(means$x, digits = 2)
N <- aggregate(data$Tumor_Sample_Barcode, by=list(data$group), FUN=length)


idea2.1.b<-ggplot(data=data)+
  geom_boxplot(aes(x=fct_inorder(group), y=total_perMB, fill=group, alpha=0.3))+
  geom_text(data = means, aes(x=Group.1, label =  paste("mean:", x), y = 70))+
  geom_text(data = N, aes(x=Group.1, label =  paste("n =", x), y = 73))+
  geom_jitter(aes(x=group, y=total_perMB))+
  theme_minimal()+   
  geom_hline(yintercept=median(as.numeric(data$total_perMB)), color="red", linetype="dashed")+
  theme(legend.position = "none")+
  labs(title = "B",
       subtitle = "Tumor Mutational Burden in Rn Rats by Group")+
  xlab(element_blank())+
  ylab("Total (per MB)")+
  scale_fill_viridis_d()+
  theme(plot.title = element_text(size = 12, face = "bold"))

grid.arrange(idea2.1.a, idea2.1.b, ncol=2)

#Idea 2.2
idea2.2.a<-ggplot(data=data)+
  geom_col(aes(x=reorder(Tumor_Sample_Barcode, -total_perMB), y=total_perMB, fill=dose))+
  theme_minimal()+   
  geom_hline(yintercept=median(data$total_perMB), color="red", linetype="dashed")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(element_blank())+
  ylab("Total (per MB)")+
  scale_fill_viridis_d()


means<-aggregate(data$total_perMB, list(data$dose), FUN=mean)
means[2]<-round(means$x, digits = 2)
N <- aggregate(data$Tumor_Sample_Barcode, by=list(data$dose), FUN=length)
compare_means(total_perMB ~ dose, data = data)

idea2.2.b<-ggplot(data=data)+
  geom_boxplot(aes(x=fct_inorder(dose), y=total_perMB, fill= dose, alpha =0.3))+
  geom_text(data = means, aes(x=Group.1, label =  paste("mean:", x), y = 70))+
  geom_text(data = N, aes(x=Group.1, label =  paste("n =", x), y = 73))+
  geom_jitter(aes(x=dose, y=total_perMB))+
  theme_minimal()+   
  geom_hline(yintercept=median(data$total_perMB), color="red", linetype="dashed")+ 
  theme(legend.position = "none")+
  ggtitle("TMB radon-rats (by dose)")+
  xlab(element_blank())+
  ylab("Total (per MB)")+
  scale_fill_viridis_d()
  

idea2.2.b + stat_compare_means(method = "t.test")
compare_means(total_perMB ~ dose, data = data)

grid.arrange(idea2.2.a, idea2.2.b, ncol=2)


#Idea 3
idea3.a<-ggplot(data=data)+
  geom_col(aes(x=reorder(Tumor_Sample_Barcode, -total_perMB), y=total_perMB, fill=duration))+
  theme_minimal()+   
  geom_hline(yintercept=median(data$total_perMB), color="red", linetype="dashed")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(element_blank())+
  ylab("Total (per MB)")+
  scale_fill_viridis_d()


means<-aggregate(data$total_perMB, list(data$duration), FUN=mean)
means[2]<-round(means$x, digits = 2)
N <- aggregate(data$Tumor_Sample_Barcode, by=list(data$duration), FUN=length)

idea3.b<-ggplot(data=data)+
  geom_boxplot(aes(x=fct_inorder(duration), y=total_perMB))+
  geom_text(data = means, aes(x=Group.1, label =  paste("mean:", x), y = 70))+
  geom_text(data = N, aes(x=Group.1, label =  paste("n =", x), y = 73))+
  geom_jitter(aes(x=duration, y=total_perMB))+
  theme_minimal()+   
  geom_hline(yintercept=median(data$total_perMB), color="red", linetype="dashed")+ 
  theme(legend.position = "none")+
  ggtitle("TMB radon-rats (by duration)")+
  xlab(element_blank())+
  ylab("Total (per MB)")+
  scale_fill_viridis_d()


grid.arrange(idea3.a, idea3.b, ncol=2)


#Idea 5
idea5.a<-ggplot(data=data)+
  geom_col(aes(x=reorder(Tumor_Sample_Barcode, -total_perMB), y=total_perMB, fill=histology))+
  theme_minimal()+   
  geom_hline(yintercept=median(data$total_perMB), color="red", linetype="dashed")+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(element_blank())+
  ylab("Total (per MB)")+
  scale_fill_viridis_d()


means<-aggregate(data$total_perMB, list(data$histology), FUN=mean)
means[2]<-round(means$x, digits = 2)
N <- aggregate(data$Tumor_Sample_Barcode, by=list(data$histology), FUN=length)

idea5.b<-ggplot(data=data)+
  geom_boxplot(aes(x=reorder(histology, -total_perMB), y=total_perMB))+
  geom_jitter(aes(x=histology, y=total_perMB))+
  theme_minimal()+   
  geom_hline(yintercept=median(data$total_perMB), color="red", linetype="dashed")+ 
  theme(legend.position = "none")+
  ggtitle("TMB radon-rats (by histology)")+
  xlab(element_blank())+
  ylab("Total (per MB)")+
  scale_fill_viridis_d()


grid.arrange(idea5.a, idea5.b, ncol=2)




#Idea 7
idea7.a<-ggplot(data=data)+
  geom_col(aes(x=reorder(Tumor_Sample_Barcode, -total_perMB), y=total_perMB, fill=estadioIV))+
  theme_minimal()+   
  geom_hline(yintercept=median(data$total_perMB), color="red", linetype="dashed")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(element_blank())+
  ylab("Total (per MB)")+
  scale_fill_viridis_d()


means<-aggregate(data$total_perMB, list(data$estadioIV), FUN=mean)
means[2]<-round(means$x, digits = 2)
N <- aggregate(data$Tumor_Sample_Barcode, by=list(data$estadioIV), FUN=length)

idea7.b<-ggplot(data=data)+
  geom_boxplot(aes(x=fct_inorder(estadioIV), y=total_perMB))+
  geom_jitter(aes(x=estadioIV, y=total_perMB))+
  theme_minimal()+   
  geom_hline(yintercept=median(data$total_perMB), color="red", linetype="dashed")+ 
  theme(legend.position = "none")+
  ggtitle("TMB radon-rats (by histology)")+
  xlab(element_blank())+
  ylab("Total (per MB)")+
  scale_fill_viridis_d()


grid.arrange(idea7.a, idea7.b, ncol=2)




#REMOVE CONTROLS

data<-data[-c(1,2,3,43),]
means<-aggregate(data$total_perMB, list(data$group), FUN=mean)
means[2]<-round(means$x, digits = 2)

box<-ggplot(data)+
  geom_boxplot(aes(x=reorder(group, -total_perMB), y=total_perMB, fill=group))+
  geom_text(data = means, aes(x=Group.1, label =  paste("Mean:", x), y = 70))+
  geom_jitter(aes(x=group, y=total_perMB))+
  theme_minimal()+ 
  theme(legend.position = "none")+
  ggtitle("TMB radon-rats")+
  xlab(element_blank())+
  ylab("Total (per MB)")



col<-ggplot(data=data)+
  geom_col(aes(x=reorder(Tumor_Sample_Barcode, -total_perMB), y=total_perMB, fill=group))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "bottom")+
  xlab(element_blank())+
  ylab("Total (per MB)")

grid.arrange(box,col,nrow=2)


ggplot(data)+
  geom_violin(aes(x=cohort,y=total_perMB))+
  geom_jitter(aes(x=cohort,y=total_perMB, color=group))+
  theme_minimal()




#BOXPLOT OF ALL GROUPS OF TMB
head(data)
a<-ggplot(data=data, aes(x=reorder(interaction(group, dose, duration, histology, papillary, estadioIV), -total_perMB), y=total_perMB, color=group))+
  geom_boxplot(lwd=0.8)+
  geom_jitter()+
  theme_minimal()+
  theme(legend.position = c(0.8,0.7), axis.text.x = element_blank(),plot.margin=unit(c(0,3,0,3), 'cm'), legend.title = element_blank())+
  ylab("TMB (per MB)")+
  xlab(element_blank())

b<-ggplot(data=data, aes(x=reorder(interaction(group, dose, duration, histology, papillary, estadioIV), -total_perMB), y=total_perMB, color=dose))+
  geom_boxplot(lwd=0.8)+
  geom_jitter()+
  theme_minimal()+
  theme(legend.position = c(0.8,0.7), axis.text.x = element_blank(),plot.margin=unit(c(0,3,0,3), 'cm'), legend.title = element_blank())+
  ylab("TMB (per MB)")+
  xlab(element_blank())

c<-ggplot(data=data, aes(x=reorder(interaction(group, dose, duration, histology, papillary, estadioIV), -total_perMB), y=total_perMB, color=duration))+
  geom_boxplot(lwd=0.8)+
  geom_jitter()+
  theme_minimal()+
  theme(legend.position = c(0.8,0.7), axis.text.x = element_blank(),plot.margin=unit(c(0,3,0,3), 'cm'), legend.title = element_blank())+
  ylab("TMB (per MB)")+
  xlab(element_blank())

d<-ggplot(data=data, aes(x=reorder(interaction(group, dose, duration, histology, papillary, estadioIV), -total_perMB), y=total_perMB, color=histology))+
  geom_boxplot(lwd=0.8)+
  geom_jitter()+
  theme_minimal()+
  theme(legend.position = c(0.8,0.7), axis.text.x = element_blank(),plot.margin=unit(c(0,3,0,3), 'cm'), legend.title = element_blank())+
  ylab("TMB (per MB)")+
  xlab(element_blank())

e<-ggplot(data=data, aes(x=reorder(interaction(group, dose, duration, histology, papillary, estadioIV), -total_perMB), y=total_perMB, color=papillary))+
  geom_boxplot(lwd=0.8)+
  geom_jitter()+
  theme_minimal()+
  theme(legend.position = c(0.8,0.7), axis.text.x = element_blank(),plot.margin=unit(c(0,3,0,3), 'cm'), legend.title = element_blank())+
  ylab("TMB (per MB)")+
  xlab(element_blank())

f<-ggplot(data=data, aes(x=reorder(interaction(group, dose, duration, histology, papillary, estadioIV), -total_perMB), y=total_perMB, color=estadioIV))+
  geom_boxplot(lwd=0.8)+
  geom_jitter()+
  theme_minimal()+
  theme(legend.position = c(0.8,0.7), axis.text.x = element_blank(),plot.margin=unit(c(0,3,0,3), 'cm'), legend.title = element_blank())+
  ylab("TMB (per MB)")+
  xlab(element_blank())


grid.arrange(a,b,c,d,e,f,nrow=3)


#interactive
library(plotly)

fig <- plot_ly(data, x = ~reorder(interaction(group, dose, duration, histology, papillary, estadioIV,drop = TRUE)), y = ~total_perMB, color = ~interaction(group,histology,dose,drop = TRUE), type = "box")
fig



  


data[data == ''] <- NA
head(data)

ggplot(data=data)+
  geom_boxplot(aes(x=group, y=total_perMB, color="group"),lwd=0.8)+
  geom_boxplot(aes(x=dose, y=total_perMB, color="dose"),lwd=0.8)+
  geom_boxplot(aes(x=duration, y=total_perMB, color="duration"),lwd=0.8)+
  geom_boxplot(data=subset(data, !is.na(papillary)),aes(x=papillary, y=total_perMB, color="adeno.papillary"),lwd=0.8)+
  geom_boxplot(data=subset(data, !is.na(estadioIV)),aes(x=estadioIV, y=total_perMB, color="estadioIV"),lwd=0.8)+
  geom_boxplot(aes(x=histology, y=total_perMB, color="histology"),lwd=0.8)+
  theme_minimal()+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 25, vjust = 1, hjust=0.9),plot.margin=unit(c(0,0,0,2), 'cm'))+
  ylab("TMB (per MB)")+
  xlab(element_blank())




#STATISTICAL TEST IDEA 6
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea6")

data<-read.csv(file="id6_group.CSV",sep=",",header=T)
subtipo<-data[,c(1,2,3,4)]
pattern<-data[,c(1,5,6,7,8,9,10,11,12,13)]

subtipo<-melt(subtipo, id.vars = unique(c("group")))
pattern<-melt(pattern, id.vars = unique(c("group")))


a<-ggplot(data=subtipo,aes(x=group, y=value, fill=variable))+
  geom_bar(position="stack", stat="identity")+
  theme_minimal()+
  ylab(element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
a
b<-ggplot(data=pattern,aes(x=group, y=value, fill=variable))+
  geom_bar(position="stack", stat="identity")+
  theme_minimal()+
  ylab(element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
b

#DOSE
data<-read.csv(file="id6_dose.CSV",sep=",",header=T)
subtipo<-data[,c(1,2,3,4)]
pattern<-data[,c(1,5,6,7,8,9,10,11,12,13)]

subtipo<-melt(subtipo, id.vars = unique(c("group")))
pattern<-melt(pattern, id.vars = unique(c("group")))

c<-ggplot(data=subtipo,aes(x=group, y=value, fill=variable))+
  geom_bar(position="stack", stat="identity")+
  theme_minimal()+
  xlab("dose")+
  ylab(element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

d<-ggplot(data=pattern,aes(x=group, y=value, fill=variable))+
  geom_bar(position="stack", stat="identity")+
  theme_minimal()+
  xlab("dose")+
  ylab(element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#DURATION
data<-read.csv(file="id6_duration.CSV",sep=",",header=T)
subtipo<-data[,c(1,2,3,4)]
pattern<-data[,c(1,5,6,7,8,9,10,11,12,13)]

subtipo<-melt(subtipo, id.vars = unique(c("group")))
pattern<-melt(pattern, id.vars = unique(c("group")))

e<-ggplot(data=subtipo,aes(x=group, y=value, fill=variable))+
  geom_bar(position="stack", stat="identity")+
  theme_minimal()+
  xlab("duration")+
  ylab(element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

f<-ggplot(data=pattern,aes(x=group, y=value, fill=variable))+
  geom_bar(position="stack", stat="identity")+
  theme_minimal()+
  xlab("duration")+
  ylab(element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


a
c
e

b
d
f

#####################################################################

adeno<-data[,c(1,2,3,4)]
adeno<-melt(adeno, id.vars = unique(c("group")))

f<-ggplot(data=pattern,aes(x=group, y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "fill")+
  geom_text(data = pattern %>% dplyr::filter(value > 0),size = 3, position = position_fill(vjust = 0.5))+
  theme_minimal()+
  theme(legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("count")+
  xlab("N")
f
#####################################################################
#IDEA 7

setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea7")

data<-read.csv(file="radon-rats-metadata.CSV",sep=",",header=T)
med1<-median(data$macroscopic_lung_lesions.mm.)

data$group <- as.character(data$group)
data$group <- factor(data$group, levels=c("D1", "D3", "D6", "D12", "Pr"))

data$duration <- as.character(data$duration)
data$duration <- factor(data$duration, levels=c("short", "D6", "long"))

a<-ggplot(data=data)+
  geom_boxplot(aes(x=group,y=macroscopic_lung_lesions.mm.,fill=group), alpha=0.4)+
  geom_hline(yintercept=med1, color="red", linetype="dashed")+
  geom_jitter(aes(x=group,y=macroscopic_lung_lesions.mm.,color=group))+
  annotate('text', x = 4, y = 29, label = 'cohort median = 12.5', color="red")+
  theme_minimal()+
  theme(panel.grid.major = element_blank())

b<-ggplot(data=data)+
  geom_boxplot(aes(x=dose,y=macroscopic_lung_lesions.mm.,fill=dose), alpha=0.4)+
  geom_hline(yintercept=med1, color="red", linetype="dashed")+
  geom_jitter(aes(x=dose,y=macroscopic_lung_lesions.mm.,color=dose))+
  annotate('text', x = 2, y = 29, label = 'cohort median = 12.5', color="red")+
  theme_minimal()+
  theme(panel.grid.major = element_blank())

c<-ggplot(data=data)+
  geom_boxplot(aes(x=duration,y=macroscopic_lung_lesions.mm.,fill=duration), alpha=0.4)+
  geom_hline(yintercept=med1, color="red", linetype="dashed")+
  geom_jitter(aes(x=duration,y=macroscopic_lung_lesions.mm.,color=duration))+
  annotate('text', x = 2.6, y = 29, label = 'cohort median = 12.5', color="red")+
  theme_minimal()+
  theme(panel.grid.major = element_blank())

grid.arrange(a,b,c,ncol=3)








med<-median(data$Rn_surv)

d<-ggplot(data=data)+
  geom_boxplot(aes(x=group,y=Rn_surv,fill=group), alpha=0.4)+
  geom_hline(yintercept=med, color="red", linetype="dashed")+
  geom_jitter(aes(x=group,y=Rn_surv,color=group))+
  annotate('text', x = 4, y = 1000, label = 'cohort median = 800.5', color="red")+
  theme_minimal()+
  theme(panel.grid.major = element_blank())+
  ylab("Rn surv (days)")

e<-ggplot(data=data)+
  geom_boxplot(aes(x=dose,y=Rn_surv,fill=dose), alpha=0.4)+
  geom_hline(yintercept=med, color="red", linetype="dashed")+
  geom_jitter(aes(x=dose,y=Rn_surv,color=dose))+
  annotate('text', x = 1.5, y = 1000, label = 'cohort median = 800.5', color="red")+
  theme_minimal()+
  theme(panel.grid.major = element_blank())+
  ylab("Rn surv (days)")

f<-ggplot(data=data)+
  geom_boxplot(aes(x=duration,y=Rn_surv,fill=duration), alpha=0.4)+
  geom_hline(yintercept=med, color="red", linetype="dashed")+
  geom_jitter(aes(x=duration,y=Rn_surv,color=duration))+
  annotate('text', x = 2.6, y = 1000, label = 'cohort median = 800.5', color="red")+
  theme_minimal()+
  theme(panel.grid.major = element_blank())+
  ylab("Rn surv (days)")

grid.arrange(d,e,f,ncol=3)

grid.arrange(a,b,c,d,e,f,ncol=3,nrow=2)

library(ggrepel)
ggplot(data=data,aes(x=Rn_surv,y=macroscopic_lung_lesions.mm., color=dose))+
  geom_point(aes(shape=factor(group)),size=3)+
  scale_shape_manual(values=c(15,18,17,16,4))+
  geom_text_repel(aes(x=Rn_surv, y=macroscopic_lung_lesions.mm., label = Sample.code))+
  geom_smooth(method='lm')+
  geom_hline(yintercept=med1, color="red", linetype="dashed")+
  geom_vline(xintercept=med, color="#00bfc4", linetype="dashed")+
  theme_minimal()+
  xlab("Rn surv (days)")+
  ylab("Macroscopic lung lesions (mm)")

head(data)


ggplot(data=data)+
  geom_col(aes(x=group,y=interaction(factor(X_T),factor(N),factor(M)),fill=group),size=2.5)+
  theme_minimal()+
  xlab("T, N ,M")


interaction(data$X_T,data$N,data$M, drop=T)
lvl<-levels(interaction(data$X_T,data$N,data$M, drop=T))

D1<-data[data$group %in% c('D1'),]
D3<-data[data$group %in% c('D3'),]
D6<-data[data$group %in% c('D6'),]
D12<-data[data$group %in% c('D12'),]
Pr<-data[data$group %in% c('Pr'),]

d1<-table(interaction(D1$X_T,D1$N,D1$M, drop=T))
d3<-table(interaction(D3$X_T,D3$N,D3$M, drop=T))
d6<-table(interaction(D6$X_T,D6$N,D6$M, drop=T))
d12<-table(interaction(D12$X_T,D12$N,D12$M, drop=T))
pr<-table(interaction(Pr$X_T,Pr$N,Pr$M, drop=T))

d1
d3
d6
d12
pr

TNM<-table(interaction(data$X_T,data$N,data$M, drop=T))


tnm<-read.table(file="TNM.txt",header=T)

tnm<-melt(tnm, id.vars = unique(c("tnm")))

tnm1<-ggplot(data=tnm,aes(x=tnm, y=value, fill=variable, label = value))+
  geom_bar(position = "fill",stat="identity")+
  geom_text(data = tnm %>% dplyr::filter(value > 0), position = position_fill(vjust = 0.5))+
  theme_minimal()+
  ylab(element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("TNM")+
  guides(fill = guide_legend(title = "Duration"))


tnm2<-ggplot(data=tnm,aes(x=tnm, y=value, fill=variable, label = value))+
  geom_bar(stat="identity")+
  geom_text(data = tnm %>% dplyr::filter(value > 0), position = position_stack(vjust = 0.5))+
  theme_minimal()+
  ylab(element_blank())+
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())+
  xlab("TNM")+
  guides(fill = guide_legend(title = "Duration"))+
  scale_y_continuous(breaks = seq(0, 14, by = 1))


grid.arrange(tnm1,tnm2,ncol=2)



data<-read.csv(file="radon-rats-metadata.CSV",sep=",",header=T)

data<-data[,c(2,7)]
data<-melt(data, id.vars = unique(c("group")))
head(data)

a<-ggplot(data=data, aes(x=value, fill=group))+
  geom_bar()+
  geom_text(stat="count", aes(label=..count..), vjust=-1)+
  theme_minimal()+
  xlab("T")


data<-read.csv(file="radon-rats-metadata.CSV",sep=",",header=T)

data<-data[,c(2,8)]
data<-melt(data, id.vars = unique(c("group")))
head(data)

b<-ggplot(data=data, aes(x=value, fill=group))+
  geom_bar()+
  geom_text(stat="count", aes(label=..count..), vjust=-1)+
  theme_minimal()+
  xlab("N")+
  scale_x_continuous(breaks = seq(0, 2, by = 1))




data<-read.csv(file="radon-rats-metadata.CSV",sep=",",header=T)

data<-data[,c(2,9)]
data<-melt(data, id.vars = unique(c("group")))
head(data)

c<-ggplot(data=data, aes(x=value, fill=group))+
  geom_bar()+
  geom_text(stat="count", aes(label=..count..), vjust=-1)+
  theme_minimal()+
  xlab("M")+
  scale_x_continuous(breaks = seq(0, 2, by = 1))


grid.arrange(a,b,c,ncol=3)



setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea7")




################################################################################
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea6")
data<-read.csv(file="adeno_histology.CSV",sep=",",header=T)


adeno<-data[c(1,2,3,4,5),]
adeno<-melt(adeno, id.vars = unique(c("group")))

a<-ggplot(data=adeno,aes(x=group, y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = adeno %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("")

adeno<-data[c(6,7),]
adeno<-melt(adeno, id.vars = unique(c("group")))

b<-ggplot(data=adeno,aes(x=group, y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = adeno %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("")

adeno<-data[c(8,9,10),]
adeno<-melt(adeno, id.vars = unique(c("group")))

c<-ggplot(data=adeno,aes(x=group, y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = adeno %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("")

data<-read.csv(file="squamous_histology.CSV",sep=",",header=T)


squamous<-data[c(1,2,3,4,5),]
squamous<-melt(squamous, id.vars = unique(c("group")))

d<-ggplot(data=squamous,aes(x=group, y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = squamous %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("")

squamous<-data[c(6,7),]
squamous<-melt(squamous, id.vars = unique(c("group")))

e<-ggplot(data=squamous,aes(x=group, y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = squamous %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("")

squamous<-data[c(8,9,10),]
squamous<-melt(squamous, id.vars = unique(c("group")))

f<-ggplot(data=squamous,aes(x=group, y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position = "stack")+
  geom_text(data = squamous %>% dplyr::filter(value > 0),size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = "none",legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("")+
  xlab("")

grid.arrange(a,b,c,d,e,f,ncol=3)



