setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea6")
getwd()
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/TMB")

library(maftools)
library(data.table)
library(rstatix)

data<-read.csv(file="TMB-radon-rats.CSV",sep=",",header=T)
rownames(data)<-data$sample

  pairwise_t_test(data=data,total_perMB ~ group, p.adjust.method = "bonferroni")
pairwise_t_test(data=data,total_perMB ~ dose, p.adjust.method = "bonferroni")
pairwise_t_test(data=data,total_perMB ~ duration, p.adjust.method = "bonferroni")
pairwise_t_test(data=data,total_perMB ~ subtipo, p.adjust.method = "bonferroni")
pairwise_t_test(data=data,total_perMB ~ papillary, p.adjust.method = "bonferroni")


compare_means(total_perMB ~ dose, data = data)
compare_means(total_perMB ~ duration, data = data)
compare_means(total_perMB ~ subtipo, data = data)
compare_means(total_perMB ~ papillary, data = data)

setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cancer_genes_filtered/plots/groups")
getwd()
load("maf.Rdata")

df<-data.frame(table(laml@data$Source_MAF))
df[3]<-c(12,9,7,4,4)
colnames(df)<-c("group", "variants", "n")
head(df)


pairwise_t_test(data=df,variants ~ group, p.adjust.method = "bonferroni")



################################################################################

#excess relative risk


library(linERR)

setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea7")

data<-read.csv(file="radon-rats-metadata.CSV",sep=",",header=T)

data<-data[,-c(1,3,4,5,7,8,9,10,11)]
data[,4]<-c(105,105,105,105,105,105,105,105,105,105,105,105,105,105,105,105,42,42,42,42,105,105,105,105,105,105,105,105,105,105,105,105,42,42,42,42)
data[,5]<-c(8.3,8.3,8.3,8.3,8.3,8.3,8.3,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1,1.7,1.7,1.7,1.7,26,26,26,26,26,26,26,26,26,26,26,26,8.9,8.9,8.9,8.9)
data[,6]<-c(13,13,13,13,13,13,13,52.1,52.1,52.1,52.1,52.1,52.1,52.1,52.1,52.1,26,26,26,26,4.3,4.3,4.3,4.3,4.3,4.3,4.3,4.3,4.3,4.3,4.3,4.3,21,21,21,21)

colnames(data)<- c("Sample","Rn_surv","group","WLM","Weekly_dose(Bq/m3)","exposure(weeks)")
head(data)

pairwise.t.test(data=data,Rn_surv ~ group, p.adjust.method = "bonferroni")



fit.1 <- fit.linERR(data=data,beta=NULL,
                    ages=cohort1[, 7:38], lag=2)
fit.linERR()


data(cohort1)
fit.1 <- fit.linERR(Surv(entryage, exitage, leu)~sex|dose1+dose2+dose3+dose4+dose5+dose6+
                      dose7+dose8+dose9+dose10+dose11+dose12+dose13+dose14+dose15+dose16+
                      dose17+dose18+dose19+dose20+dose21+dose22+dose23+dose24+dose25+dose26+
                      dose27+dose28+dose29+dose30+dose31+dose32, data=cohort1, beta=NULL,
                    ages=cohort1[, 7:38], lag=2)
ERRci(fit.1, prob=0.9)
summary(fit.1)

library(swimplot)
library(ggplot2)
library(tidyverse)
library(scales)
library(grid)

df<-read.csv(file="swimming_plot_surv.CSV",sep=",",header=T)
head(df)
#colnames(df)<-c("sample","date_of_autopsy","DS","cause_of_DS","Birth_surv","Rn_surv")

ggplot(df, aes(x = 0, xend = Rn.surv., y = Sample.code, yend = Sample.code, color = D.S)) + 
  geom_segment(size = 4) +
  theme_classic() + 
  labs(x = "Rn survival", y = "Sample code", color = "D.S") + 
  ggtitle("Swimmer plot of survival by sample")

grid.force()
# change shape of arrows
grid.gedit("segments", gp=gpar(linejoin ='mitre'))
grid.gedit("layout", gp=gpar(linejoin ='mitre'))




df$date.of.autopsy <- as.Date(df$date.of.autopsy, format = "%m/%d/%Y")

ggplot(df, aes(x = date.of.autopsy, xend = date.of.autopsy, y = Sample.code, yend = Sample.code, color = D.S)) + 
  geom_segment(aes(x = date.of.autopsy, xend = date.of.autopsy, y = Birth.surv., yend = Rn.surv.), size = 4) + 
  theme_classic() + 
  labs(x = "Date of Autopsy", y = "Sample code", color = "D.S") + 
  ggtitle("Swimmer plot of survival by sample") +
  scale_x_date(date_labels = "%b %Y")







