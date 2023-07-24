#MIRROR
library("factoextra")
library("missMDA")
library("FactoMineR")
library("PCAmixdata")
library("NanoTube")
library(org.Hs.eg.db)
library(RColorBrewer)
library(gplots)
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/")

################################################################################
################################################################################
################################################################################

#TRANSCRIPTOMICS ANALYSIS
#prepare data
files<-list.files("C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/RCC/")
mirror_number<-c("2","9","15","1","14-1","16","14-2","12","3","4","6","37","21","20","26","36","18","32","22","35","28","30","34","77","79","90","99","101","92","HCUV02","HCUV03","HCUV04","HCUV07","HCUV09","HCUV12","HCUV14","HCUV15","HCUV16","HCUV17","HCUV18","HCUV20","HCUV21","HCUV22","HCUV23","HCUV24","HCUV27","HCUV28","126","HCUV31","HCUV33","HCUV34","HCUV35","HCUV36","HCUV37","HCUV38","HCUV40")

diagnosis <- c("Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenosquamous","Squamous","Adenocarcinoma","Squamous","Mixte","Adenocarcinoma","Adenocarcinoma","Squamous","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Squamous","Adenocarcinoma","Squamous","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Squamous","Squamous","Squamous","Squamous","Adenocarcinoma","Squamous","Squamous","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Squamous","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma")

df<-data.frame(
  Sample_Title = paste("Mirror patient", rep(1:56)),
  GEO_Accession = sub(".RCC", "", sub(".*?_.*?_(.*)", "\\1", files)),
  Sample_Status = "Private on Jul 5 2023",
  Sample_Type = "RNA",
  Sample_Source = "tissue RNA",
  Sample_Diagnosis = paste("Lung Cancer", diagnosis),
  RCC_Name = files)
write.csv(df, file = "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/MIRROR_sample_data.csv")

example_data <- "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/RCC/"
sample_info <- "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/MIRROR_sample_data.csv"

dat <- processNanostringData(example_data,
                             sampleTab = sample_info,
                             idCol = "RCC_Name",
                             groupCol = "Sample_Diagnosis",
                             normalization = "nSolver", 
                             bgType = "t.test", 
                             bgPVal = 0.01)

ruv::ruv_rle(t(log2(exprs(dat))), ylim = c(-2, 2))+ggtitle("Non normalized data")

table(dat$groups)
limmaResults <- runLimmaAnalysis(dat, base.group = "Lung Cancer Adenocarcinoma")

limmaStats <- makeDiffExprFile(limmaResults, filename = NULL, returns = "stats")
limmaStats <- as.data.frame(limmaStats)
write.table(limmaStats, file = "limma_results_histology.tsv")

# Rounding for clarity
limmaTab <- head(limmaStats[order(limmaStats$`p-val (Lung.Cancer.Squamous)`, 
                                  decreasing = FALSE), 1:4])
limmaTab[,1] <- format(limmaTab[,1], digits = 2, nsmall = 1)
limmaTab[,3] <- format(limmaTab[,3], digits = 1, scientific = TRUE)
limmaTab[,4] <- format(limmaTab[,4], digits = 1, nsmall = 1)

# Order by lowest to highest p-value for 'Autoimmune Retinopathy' vs. 'None'
knitr::kable(head(limmaTab), 
             row.names = TRUE, format = "html", align = "c")


limmaStats$diffexpressed <- "NO"
limmaStats$diffexpressed[limmaStats$`Log2FC (Lung.Cancer.Squamous)` > 0.5 & limmaStats$`p-val (Lung.Cancer.Squamous)` < 0.05] <- "UP"
limmaStats$diffexpressed[limmaStats$`Log2FC (Lung.Cancer.Squamous)`  < -0.5 &  limmaStats$`p-val (Lung.Cancer.Squamous)` < 0.05] <- "DOWN"

limmaStats$delabel <- NA
non_missing_genes <- !is.na(rownames(limmaStats))
limmaStats$delabel[limmaStats$diffexpressed != "NO" & non_missing_genes] <- rownames(limmaStats)[limmaStats$diffexpressed != "NO" & non_missing_genes]

svg(paste0("volcano_",'Lung.Cancer.Squamous.svg',sep="."), height = 5, width = 9)
ggplot(data=limmaStats, aes(x=`Log2FC (Lung.Cancer.Squamous)`, y=-log10(`p-val (Lung.Cancer.Squamous)`), col = limmaStats$diffexpressed, label=delabel)) +
  geom_point()+
  geom_text(size = 2.5, vjust=-0.8)+
  geom_vline(xintercept=c(-0.5, 0.5), col="red", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype= "dashed") +
  theme_minimal()+
  xlab("log2(Fold Change)")+
  ylab("-log10(pval)")+
  labs(color='Expression regulation')+
  ggtitle("Histology: Adenocarcinoma vs Squamous")
dev.off()




#normalize the data
exprs(dat)=log(exprs(dat)+1, 2)

std<-apply(dat@assayData$exprs, 1, sd)
summary(std)

plot <- nanostringPCA(dat[std > quantile(std, 0.25),], interactive.plot = TRUE)$pl
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[2], 2), ", N genes = ", nrow(dat[std > quantile(std, 0.25),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Histology_1qt.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

plot <- nanostringPCA(dat[std>mean(std),], interactive.plot = TRUE)$plt
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[4], 2), ", N genes = ", nrow(dat[std > mean(std),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Histology_mean.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

plot <- nanostringPCA(dat[std>quantile(std, 0.75),], interactive.plot = TRUE)$plt
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[5], 2), ", N genes = ", nrow(dat[std > quantile(std, 0.75),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Histology_3qt.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

ruv::ruv_rle(t(log2(exprs(dat))), ylim = c(-2, 2))+ggtitle("Normalized data")

#QUALITY CONTROL
dat <- processNanostringData(example_data, 
                             idCol = "RCC_Name",
                             output.format = "list",
                             includeQC = TRUE)

head(dat$qc)[,1:5]

posQC <- positiveQC(dat)
knitr::kable(head(posQC$tab), 
             row.names = FALSE, format = "html", align = "c", digits = 2)

posQC2 <- positiveQC(dat, samples = 1:6)

posQC2$plt

negQC <- negativeQC(dat, interactive.plot = FALSE)

knitr::kable(head(negQC$tab), 
             row.names = TRUE, format = "html", align = "c")

negQC$plt

signif(head(dat$hk.scalefactors), digits = 2)


################################################################################
################################################################################
################################################################################
diagnosis <- c("low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","high","low","low","low","low","low","low","low","low","low","low","low","low","high","low","low","high","low","low","low","low","low","low","low","low")

df<-data.frame(
  Sample_Title = paste("Mirror patient", rep(1:56)),
  GEO_Accession = sub(".RCC", "", sub(".*?_.*?_(.*)", "\\1", files)),
  Sample_Status = "Private on Jul 5 2023",
  Sample_Type = "RNA",
  Sample_Source = "tissue RNA",
  Sample_Diagnosis = paste("Radon Dose", diagnosis),#*****
  RCC_Name = files)

write.csv(df, file = "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/MIRROR_sample_data.csv")
sample_info <- "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/MIRROR_sample_data.csv"

dat <- processNanostringData(example_data,
                             sampleTab = sample_info,
                             idCol = "RCC_Name",
                             groupCol = "Sample_Diagnosis",
                             normalization = "nSolver", 
                             bgType = "t.test", 
                             bgPVal = 0.01)



table(dat$groups)


limmaResults <- runLimmaAnalysis(dat, base.group = "Radon Dose low")#*****

limmaStats <- makeDiffExprFile(limmaResults, filename = NULL, returns = "stats")
limmaStats <- as.data.frame(limmaStats)
write.table(limmaStats, file = "limma_results_radon_dose.tsv", sep = "\t")

limmaTab <- head(limmaStats[order(limmaStats$`p-val (Radon.Dose.high)`, #*****
                                  decreasing = FALSE), 1:4])
limmaTab[,1] <- format(limmaTab[,1], digits = 2, nsmall = 1)
limmaTab[,3] <- format(limmaTab[,3], digits = 1, scientific = TRUE)
limmaTab[,4] <- format(limmaTab[,4], digits = 1, nsmall = 1)

knitr::kable(head(limmaTab), 
             row.names = TRUE, format = "html", align = "c")



limmaStats$diffexpressed <- "NO"
limmaStats$diffexpressed[limmaStats$`Log2FC (Radon.Dose.high)` > 0.5 & limmaStats$`p-val (Radon.Dose.high)` < 0.05] <- "UP"
limmaStats$diffexpressed[limmaStats$`Log2FC (Radon.Dose.high)`  < -0.5 &  limmaStats$`p-val (Radon.Dose.high)` < 0.05] <- "DOWN"

limmaStats$delabel <- NA
non_missing_genes <- !is.na(rownames(limmaStats))
limmaStats$delabel[limmaStats$diffexpressed != "NO" & non_missing_genes] <- rownames(limmaStats)[limmaStats$diffexpressed != "NO" & non_missing_genes]

svg(paste0("volcano_",'Radon.Dose.high.svg',sep="."), height = 5, width = 9)
ggplot(data=limmaStats, aes(x=`Log2FC (Radon.Dose.high)`, y=-log10(`p-val (Radon.Dose.high)`), col = limmaStats$diffexpressed, label=delabel)) +
  geom_point()+
  geom_text(size = 2.5, vjust=-0.8)+
  geom_vline(xintercept=c(-0.5, 0.5), col="red", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype= "dashed") +
  theme_minimal()+
  xlab("log2(Fold Change)")+
  ylab("-log10(pval)")+
  labs(color='Expression regulation')+
  ggtitle("Radon Dose: low vs high")
dev.off()


#normalize the data
exprs(dat)=log(exprs(dat)+1, 2)

std<-apply(dat@assayData$exprs, 1, sd)
summary(std)

plot <- nanostringPCA(dat[std > quantile(std, 0.25),], interactive.plot = TRUE)$pl
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[2], 2), ", N genes = ", nrow(dat[std > quantile(std, 0.25),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Radon_Dose_1qt.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

plot <- nanostringPCA(dat[std>mean(std),], interactive.plot = TRUE)$plt
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[4], 2), ", N genes = ", nrow(dat[std > mean(std),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Radon_Dose_mean.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

plot <- nanostringPCA(dat[std>quantile(std, 0.75),], interactive.plot = TRUE)$plt
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[5], 2), ", N genes = ", nrow(dat[std > quantile(std, 0.75),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Radon_Dose_3qt.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)


###############################################################################
################################################################################
################################################################################
diagnosis <- c("smoking","non smoking","smoking","non smoking","non smoking","smoking","non smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","non smoking","smoking","non smoking","non smoking","smoking","smoking","smoking","smoking","smoking","smoking","non smoking","non smoking","radon + smoking","smoking","smoking","smoking","smoking","smoking","smoking","non smoking","smoking","smoking","smoking","smoking","smoking","radon + smoking","smoking","smoking","radon + smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking")

df<-data.frame(
  Sample_Title = paste("Mirror patient", rep(1:56)),
  GEO_Accession = sub(".RCC", "", sub(".*?_.*?_(.*)", "\\1", files)),
  Sample_Status = "Private on Jul 5 2023",
  Sample_Type = "RNA",
  Sample_Source = "tissue RNA",
  Sample_Diagnosis = paste("Carcinogen Group", diagnosis),#*****
  RCC_Name = files)

write.csv(df, file = "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/MIRROR_sample_data.csv")
sample_info <- "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/MIRROR_sample_data.csv"

dat <- processNanostringData(example_data,
                             sampleTab = sample_info,
                             idCol = "RCC_Name",
                             groupCol = "Sample_Diagnosis",
                             normalization = "nSolver", 
                             bgType = "t.test", 
                             bgPVal = 0.01)


table(dat$groups)

limmaResults <- runLimmaAnalysis(dat, base.group = "Carcinogen Group smoking")#*****Carcinogen Group non smoking

limmaStats <- makeDiffExprFile(limmaResults, filename = NULL, returns = "stats")
limmaStats <- as.data.frame(limmaStats)
write.table(limmaStats, file = "limma_results_carcinogen_group.tsv", sep = "\t")

limmaTab <- head(limmaStats[order(limmaStats$`p-val (Carcinogen.Group.radon.+.smoking)`, #*****
                                  decreasing = FALSE), 1:4])
limmaTab[,1] <- format(limmaTab[,1], digits = 2, nsmall = 1)
limmaTab[,3] <- format(limmaTab[,3], digits = 1, scientific = TRUE)
limmaTab[,4] <- format(limmaTab[,4], digits = 1, nsmall = 1)

knitr::kable(head(limmaTab), 
             row.names = TRUE, format = "html", align = "c")



limmaStats$diffexpressed <- "NO"
limmaStats$diffexpressed[limmaStats$`Log2FC (Carcinogen.Group.radon.+.smoking)` > 0.5 & limmaStats$`p-val (Carcinogen.Group.radon.+.smoking)` < 0.05] <- "UP"
limmaStats$diffexpressed[limmaStats$`Log2FC (Carcinogen.Group.radon.+.smoking)`  < -0.5 &  limmaStats$`p-val (Carcinogen.Group.radon.+.smoking)` < 0.05] <- "DOWN"

limmaStats$delabel <- NA
non_missing_genes <- !is.na(rownames(limmaStats))
limmaStats$delabel[limmaStats$diffexpressed != "NO" & non_missing_genes] <- rownames(limmaStats)[limmaStats$diffexpressed != "NO" & non_missing_genes]

svg(paste0("volcano_",'Carcinogen.Group.radon.+.smoking.svg',sep="."), height = 5, width = 9)
ggplot(data=limmaStats, aes(x=`Log2FC (Carcinogen.Group.radon.+.smoking)`, y=-log10(`p-val (Carcinogen.Group.radon.+.smoking)`), col = limmaStats$diffexpressed, label=delabel)) +
  geom_point()+
  geom_text(size = 2.5, vjust=-0.8)+
  geom_vline(xintercept=c(-0.5, 0.5), col="red", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype= "dashed") +
  theme_minimal()+
  xlab("log2(Fold Change)")+
  ylab("-log10(pval)")+
  labs(color='Expression regulation')+
  ggtitle("Radon Carcinogen: Smoker vs Radon + Smoking")
dev.off()


#normalize the data
exprs(dat)=log(exprs(dat)+1, 2)

std<-apply(dat@assayData$exprs, 1, sd)
summary(std)

plot <- nanostringPCA(dat[std > quantile(std, 0.25),], interactive.plot = TRUE)$pl
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[2], 2), ", N genes = ", nrow(dat[std > quantile(std, 0.25),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Radon_Carcinogen_1qt.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

plot <- nanostringPCA(dat[std>mean(std),], interactive.plot = TRUE)$plt
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[4], 2), ", N genes = ", nrow(dat[std > mean(std),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Radon_Carcinogen_mean.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

plot <- nanostringPCA(dat[std>quantile(std, 0.75),], interactive.plot = TRUE)$plt
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[5], 2), ", N genes = ", nrow(dat[std > quantile(std, 0.75),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Radon_Carcinogen_3qt.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)


###############################################################################
################################################################################
################################################################################
diagnosis <- c("yes","no","yes","no","no","yes","no","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","no","yes","no","no","yes","yes","yes","yes","yes","yes","no","no","yes","yes","yes","yes","yes","yes","yes","no","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes")

df<-data.frame(
  Sample_Title = paste("Mirror patient", rep(1:56)),
  GEO_Accession = sub(".RCC", "", sub(".*?_.*?_(.*)", "\\1", files)),
  Sample_Status = "Private on Jul 5 2023",
  Sample_Type = "RNA",
  Sample_Source = "tissue RNA",
  Sample_Diagnosis = paste("Smoking Habit", diagnosis),#*****
  RCC_Name = files)

write.csv(df, file = "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/MIRROR_sample_data.csv")
sample_info <- "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/MIRROR_sample_data.csv"

dat <- processNanostringData(example_data,
                             sampleTab = sample_info,
                             idCol = "RCC_Name",
                             groupCol = "Sample_Diagnosis",
                             normalization = "nSolver", 
                             bgType = "t.test", 
                             bgPVal = 0.01)


table(dat$groups)

limmaResults <- runLimmaAnalysis(dat, base.group = "Smoking Habit no")#*****

limmaStats <- makeDiffExprFile(limmaResults, filename = NULL, returns = "stats")
limmaStats <- as.data.frame(limmaStats)
write.table(limmaStats, file = "limma_results_smoking_habit.tsv", sep = "\t")

limmaTab <- head(limmaStats[order(limmaStats$`p-val (Smoking.Habit.yes)`, #*****
                                  decreasing = FALSE), 1:4])
limmaTab[,1] <- format(limmaTab[,1], digits = 2, nsmall = 1)
limmaTab[,3] <- format(limmaTab[,3], digits = 1, scientific = TRUE)
limmaTab[,4] <- format(limmaTab[,4], digits = 1, nsmall = 1)

knitr::kable(head(limmaTab), 
             row.names = TRUE, format = "html", align = "c")



limmaStats$diffexpressed <- "NO"
limmaStats$diffexpressed[limmaStats$`Log2FC (Smoking.Habit.yes)` > 0.5 & limmaStats$`p-val (Smoking.Habit.yes)` < 0.05] <- "UP"
limmaStats$diffexpressed[limmaStats$`Log2FC (Smoking.Habit.yes)`  < -0.5 &  limmaStats$`p-val (Smoking.Habit.yes)` < 0.05] <- "DOWN"

limmaStats$delabel <- NA
non_missing_genes <- !is.na(rownames(limmaStats))
limmaStats$delabel[limmaStats$diffexpressed != "NO" & non_missing_genes] <- rownames(limmaStats)[limmaStats$diffexpressed != "NO" & non_missing_genes]

svg(paste0("volcano_",'Smoking.Habit.yes.svg',sep="."), height = 5, width = 9)
ggplot(data=limmaStats, aes(x=`Log2FC (Smoking.Habit.yes)`, y=-log10(`p-val (Smoking.Habit.yes)`), col = limmaStats$diffexpressed, label=delabel)) +
  geom_point()+
  geom_text(size = 2.5, vjust=-0.8)+
  geom_vline(xintercept=c(-0.5, 0.5), col="red", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype= "dashed") +
  theme_minimal()+
  xlab("log2(Fold Change)")+
  ylab("-log10(pval)")+
  labs(color='Expression regulation')+
  ggtitle("Smoking Habit: no vs yes")
dev.off()


#normalize the data
exprs(dat)=log(exprs(dat)+1, 2)

std<-apply(dat@assayData$exprs, 1, sd)
summary(std)

plot <- nanostringPCA(dat[std > quantile(std, 0.25),], interactive.plot = TRUE)$pl
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[2], 2), ", N genes = ", nrow(dat[std > quantile(std, 0.25),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Smoking_Habit_1qt.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

plot <- nanostringPCA(dat[std>mean(std),], interactive.plot = TRUE)$plt
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[4], 2), ", N genes = ", nrow(dat[std > mean(std),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Smoking_Habit_mean.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

plot <- nanostringPCA(dat[std>quantile(std, 0.75),], interactive.plot = TRUE)$plt
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[5], 2), ", N genes = ", nrow(dat[std > quantile(std, 0.75),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Smoking_Habit_3qt.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

###############################################################################
################################################################################
################################################################################
diagnosis <- c("heavy","non_smoking","heavy","non_smoking","non_smoking","very heavy","non_smoking","heavy","heavy","medium","heavy","heavy","heavy","medium","light","heavy","very heavy","medium","heavy","non_smoking","heavy","non_smoking","non_smoking","medium","medium","heavy","medium","medium","heavy","non_smoking","non_smoking","light","medium","heavy","heavy","heavy","heavy","heavy","non_smoking","heavy","heavy","heavy","heavy","medium","heavy","medium","heavy","heavy","light","light","medium","heavy","heavy","heavy","heavy","heavy")

df<-data.frame(
  Sample_Title = paste("Mirror patient", rep(1:56)),
  GEO_Accession = sub(".RCC", "", sub(".*?_.*?_(.*)", "\\1", files)),
  Sample_Status = "Private on Jul 5 2023",
  Sample_Type = "RNA",
  Sample_Source = "tissue RNA",
  Sample_Diagnosis = paste("Smoking Habit", diagnosis),#*****
  RCC_Name = files)

write.csv(df, file = "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/MIRROR_sample_data.csv")
sample_info <- "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/MIRROR_sample_data.csv"

dat <- processNanostringData(example_data,
                             sampleTab = sample_info,
                             idCol = "RCC_Name",
                             groupCol = "Sample_Diagnosis",
                             normalization = "nSolver", 
                             bgType = "t.test", 
                             bgPVal = 0.01)


table(dat$groups)

limmaResults <- runLimmaAnalysis(dat, base.group = "Smoking Habit non_smoking")#*****

limmaStats <- makeDiffExprFile(limmaResults, filename = NULL, returns = "stats")
limmaStats <- as.data.frame(limmaStats)
write.table(limmaStats, file = "limma_results_smoking_habit2.tsv", sep = "\t")

limmaTab <- head(limmaStats[order(limmaStats$`p-val (Smoking.Habit.heavy)`, #*****
                                  decreasing = FALSE), 1:4])
limmaTab[,1] <- format(limmaTab[,1], digits = 2, nsmall = 1)
limmaTab[,3] <- format(limmaTab[,3], digits = 1, scientific = TRUE)
limmaTab[,4] <- format(limmaTab[,4], digits = 1, nsmall = 1)

knitr::kable(head(limmaTab), 
             row.names = TRUE, format = "html", align = "c")



limmaStats$diffexpressed <- "NO"
limmaStats$diffexpressed[limmaStats$`Log2FC (Smoking.Habit.heavy)` > 0.5 & limmaStats$`p-val (Smoking.Habit.heavy)` < 0.05] <- "UP"
limmaStats$diffexpressed[limmaStats$`Log2FC (Smoking.Habit.heavy)`  < -0.5 &  limmaStats$`p-val (Smoking.Habit.heavy)` < 0.05] <- "DOWN"

limmaStats$delabel <- NA
non_missing_genes <- !is.na(rownames(limmaStats))
limmaStats$delabel[limmaStats$diffexpressed != "NO" & non_missing_genes] <- rownames(limmaStats)[limmaStats$diffexpressed != "NO" & non_missing_genes]

svg(paste0("volcano_",'Smoking.Habit.heavy.svg',sep="."), height = 5, width = 9)
ggplot(data=limmaStats, aes(x=`Log2FC (Smoking.Habit.heavy)`, y=-log10(`p-val (Smoking.Habit.heavy)`), col = limmaStats$diffexpressed, label=delabel)) +
  geom_point()+
  geom_text(size = 2.5, vjust=-0.8)+
  geom_vline(xintercept=c(-0.5, 0.5), col="red", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype= "dashed") +
  theme_minimal()+
  xlab("log2(Fold Change)")+
  ylab("-log10(pval)")+
  labs(color='Expression regulation')+
  ggtitle("Smoking Habit: no vs heavy")
dev.off()


#normalize the data
exprs(dat)=log(exprs(dat)+1, 2)

std<-apply(dat@assayData$exprs, 1, sd)
summary(std)

plot <- nanostringPCA(dat[std > quantile(std, 0.25),], interactive.plot = TRUE)$pl
plot <- plotly::layout(plot, title = "sd > 1st quartile")
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Smoking_Habit2_1qt.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

plot <- nanostringPCA(dat[std>mean(std),], interactive.plot = TRUE)$plt
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[4], 2), ", N genes = ", nrow(dat[std > mean(std),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Smoking_Habit2_mean.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

plot <- nanostringPCA(dat[std>quantile(std, 0.75),], interactive.plot = TRUE)$plt
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[5], 2), ", N genes = ", nrow(dat[std > quantile(std, 0.75),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'Smoking_Habit2_3qt.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)



###############################################################################
################################################################################
################################################################################
diagnosis <- c("low","low","low","low","mod","low","low","mod","mod","mod","high","low","low","low","mod","mod","low","mod","low","mod","low","low","mod","mod","low","high","mod","high","low","low","low","low","mod","low","low","low","low","high","low","low","low","low","low","low","high","low","low","low","mod","low","high","low","high","low","high","low")

df<-data.frame(
  Sample_Title = paste("Mirror patient", rep(1:56)),
  GEO_Accession = sub(".RCC", "", sub(".*?_.*?_(.*)", "\\1", files)),
  Sample_Status = "Private on Jul 5 2023",
  Sample_Type = "RNA",
  Sample_Source = "tissue RNA",
  Sample_Diagnosis = paste("TILs", diagnosis),#*****
  RCC_Name = files)

write.csv(df, file = "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/MIRROR_sample_data.csv")
sample_info <- "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/MIRROR_sample_data.csv"

dat <- processNanostringData(example_data,
                             sampleTab = sample_info,
                             idCol = "RCC_Name",
                             groupCol = "Sample_Diagnosis",
                             normalization = "nSolver", 
                             bgType = "t.test", 
                             bgPVal = 0.01)


table(dat$groups)

limmaResults <- runLimmaAnalysis(dat, base.group = "TILs low")#*****

limmaStats <- makeDiffExprFile(limmaResults, filename = NULL, returns = "stats")
limmaStats <- as.data.frame(limmaStats)
write.table(limmaStats, file = "limma_results_TILs.tsv", sep = "\t")#*****

limmaTab <- head(limmaStats[order(limmaStats$`p-val (TILs.high)`, #*****
                                  decreasing = FALSE), 1:4])
limmaTab[,1] <- format(limmaTab[,1], digits = 2, nsmall = 1)
limmaTab[,3] <- format(limmaTab[,3], digits = 1, scientific = TRUE)
limmaTab[,4] <- format(limmaTab[,4], digits = 1, nsmall = 1)

knitr::kable(head(limmaTab), 
             row.names = TRUE, format = "html", align = "c")



limmaStats$diffexpressed <- "NO"
limmaStats$diffexpressed[limmaStats$`Log2FC (TILs.high)` > 0.5 & limmaStats$`p-val (TILs.high)` < 0.05] <- "UP"
limmaStats$diffexpressed[limmaStats$`Log2FC (TILs.high)`  < -0.5 &  limmaStats$`p-val (TILs.high)` < 0.05] <- "DOWN"

limmaStats$delabel <- NA
non_missing_genes <- !is.na(rownames(limmaStats))
limmaStats$delabel[limmaStats$diffexpressed != "NO" & non_missing_genes] <- rownames(limmaStats)[limmaStats$diffexpressed != "NO" & non_missing_genes]

svg(paste0("volcano_",'TILs.yes.svg',sep="."), height = 5, width = 9)
ggplot(data=limmaStats, aes(x=`Log2FC (TILs.high)`, y=-log10(`p-val (TILs.high)`), col = limmaStats$diffexpressed, label=delabel)) +
  geom_point()+
  geom_text(size = 2.5, vjust=-0.8)+
  geom_vline(xintercept=c(-0.5, 0.5), col="red", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype= "dashed") +
  theme_minimal()+
  xlab("log2(Fold Change)")+
  ylab("-log10(pval)")+
  labs(color='Expression regulation')+
  ggtitle("TILs: low vs high")
dev.off()


#normalize the data
exprs(dat)=log(exprs(dat)+1, 2)

std<-apply(dat@assayData$exprs, 1, sd)
summary(std)

plot <- nanostringPCA(dat[std > quantile(std, 0.25),], interactive.plot = TRUE)$pl
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[2], 2), ", N genes = ", nrow(dat[std > quantile(std, 0.25),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'TILs_1qt.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

plot <- nanostringPCA(dat[std>mean(std),], interactive.plot = TRUE)$plt
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[4], 2), ", N genes = ", nrow(dat[std > mean(std),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'TILs_mean.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)

plot <- nanostringPCA(dat[std>quantile(std, 0.75),], interactive.plot = TRUE)$plt
plot <- plotly::layout(plot, title = paste0("sd > ", round(summary(std)[5], 2), ", N genes = ", nrow(dat[std > quantile(std, 0.75),])))
htmlwidgets::saveWidget(
  widget = plot, #the plotly object
  file = paste0("PCA_",'TILs_3qt.html'), #the path & file name
  selfcontained = TRUE #creates a single html file
)


#############################################################

#############################################################

#############################################################













#HEATMAP
limmaStats <- limmaStats[limmaStats$diffexpressed != "NO", ]
limmaStats <- limmaStats[limmaStats$`p-val (Radon.Dose.high)` < 0.01, ]

genes<-limmaStats$delabel

ref<-read.table(file = "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/ref.RCC", sep = ",", header = T)
ref<-subset(ref, Name %in% genes)[,2:3]

rownames_dat <- rownames(dat@assayData$exprs)
subset_exprs <- dat@assayData$exprs[rownames_dat %in% ref$Accession, ]
subset_exprs<-log(subset_exprs+1, 2)

# Find the indices of matching Accession values between subset_exprs and ref
indices <- match(rownames(subset_exprs), ref$Accession)
rownames(subset_exprs) <- ref$Name[indices]

#RADON DOSE
diagnosis <- c("low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","high","low","low","low","low","low","low","low","low","low","low","low","low","high","low","low","high","low","low","low","low","low","low","low","low")

cols<-paste(mirror_number, diagnosis, sep = ".")

colnames(subset_exprs) <- cols


diagnosis <- factor(diagnosis, levels = c("high", "low"))
subset_exprs <- subset_exprs[, order(diagnosis)]

pal <- colorRampPalette(brewer.pal(11, "RdYlGn"))(100)
heatmap.2(as.matrix(subset_exprs), col = pal, trace = "none", key.xlab = "log2(count)", key.title = "pval<0.01", dendrogram = "none", Colv = FALSE)


#GENERATE HEATMAP WITYH THE DATA
heatmap.2(as.matrix(subset_exprs), col = pal, trace = "none", key.xlab="log2(count)", key.title = "pval<0.01")





#HEATMAP
df <- data.frame(gene = rownames(limmaStats),
                 pval1 = limmaStats$`p-val (Carcinogen.Group.non.smoking)`,
                 pval2 = limmaStats$`p-val (Carcinogen.Group.radon.+.smoking)`
)

filtered_df <- subset(df, apply(df[, -1], 1, function(x) any(x < 0.05)))


genes<-filtered_df$gene

ref<-read.table(file = "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/ref.RCC", sep = ",", header = T)
ref<-subset(ref, Name %in% genes)[,2:3]

rownames_dat <- rownames(dat@assayData$exprs)
subset_exprs <- dat@assayData$exprs[rownames_dat %in% ref$Accession, ]
subset_exprs<-log(subset_exprs+1, 2)

# Find the indices of matching Accession values between subset_exprs and ref
indices <- match(rownames(subset_exprs), ref$Accession)
rownames(subset_exprs) <- ref$Name[indices]

#SMOKING 
diagnosis <- c("smoking","non smoking","smoking","non smoking","non smoking","smoking","non smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","non smoking","smoking","non smoking","non smoking","smoking","smoking","smoking","smoking","smoking","smoking","non smoking","non smoking","radon + smoking","smoking","smoking","smoking","smoking","smoking","smoking","non smoking","smoking","smoking","smoking","smoking","smoking","radon + smoking","smoking","smoking","radon + smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking","smoking")
cols<-paste(mirror_number, diagnosis, sep = ".")

#TILS
diagnosis <- c("low","low","low","low","mod","low","low","mod","mod","mod","high","low","low","low","mod","mod","low","mod","low","mod","low","low","mod","mod","low","high","mod","high","low","low","low","low","mod","low","low","low","low","high","low","low","low","low","low","low","high","low","low","low","mod","low","high","low","high","low","high","low")
cols<-paste(mirror_number, diagnosis, sep = ".")

colnames(subset_exprs) <- cols


diagnosis <- factor(diagnosis, levels = c("radon + smoking", "smoking", "non smoking"))
subset_exprs <- subset_exprs[, order(diagnosis)]

pal <- colorRampPalette(brewer.pal(11, "RdYlGn"))(100)
#GENERATE HEATMAP WITYH THE DATA
heatmap.2(as.matrix(subset_exprs), col = pal, trace = "none", key.xlab="log2(count)", key.title = "pval<0.05", margins=c(7,5))

heatmap.2(as.matrix(subset_exprs), col = pal, trace = "none", key.xlab = "log2(count)", key.title = "pval<0.05", dendrogram = "none", Colv = FALSE, margins=c(9,5))






















#HEATMAP FILTER STANDARD DEVIATION AND STUFF
ref<-read.table(file = "C:/Users/Bedmar/Desktop/IDIBAPS/MIRROR/ref.RCC", sep = ",", header = T)
ref<-subset(ref, Name %in% genes)[,2:3]

rownames_dat <- rownames(dat@assayData$exprs)
subset_exprs <- dat@assayData$exprs[rownames_dat %in% ref$Accession, ]
subset_exprs<-log(subset_exprs+1, 2)

# Find the indices of matching Accession values between subset_exprs and ref
indices <- match(rownames(subset_exprs), ref$Accession)
rownames(subset_exprs) <- ref$Name[indices]

#RADON DOSE
diagnosis <- c("low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","low","high","low","low","low","low","low","low","low","low","low","low","low","low","high","low","low","high","low","low","low","low","low","low","low","low")

#smoking habit
diagnosis <- c("yes","no","yes","no","no","yes","no","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","no","yes","no","no","yes","yes","yes","yes","yes","yes","no","no","yes","yes","yes","yes","yes","yes","yes","no","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes","yes")

cols<-paste(mirror_number, diagnosis, sep = ".")

colnames(subset_exprs) <- cols


diagnosis <- factor(diagnosis, levels = c("yes", "no"))
subset_exprs <- subset_exprs[, order(diagnosis)]

pal <- colorRampPalette(brewer.pal(11, "RdYlGn"))(100)
heatmap.2(as.matrix(subset_exprs), col = pal, trace = "none", key.xlab = "log2(count)", key.title = "", dendrogram = "none", Colv = FALSE)


#GENERATE HEATMAP WITYH THE DATA
heatmap.2(as.matrix(subset_exprs), col = pal, trace = "none", key.xlab="log2(count)", key.title = "")















































go_annotations <- select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", columns = "GO")
printGraph(GOdata, resultWeight, firstSigNodes = 10, resultFis, fn.prefix = "tGO", useInfo = "def")

library(topGO)
library(ALL)

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)



GO_cancer <- c('GO:0030307', # Positive regulation of cell growth
               'GO:0008284', # Positive regulation of cell proliferation
               'GO:0045787', # Positivie regulation of cell cycle
               'GO:0045786', # negative regulation of cell cycle
               'GO:0008283', # Cell Proliferation
               'GO:0007049', # Cell Cycle
               'GO:0051301', # Cell Division
               'GO:0051781', # Positive Regulation of Cell Division
               'GO:0009968', # Negative regulation of signal transduction
               'GO:0030308', # Negative regulation of cell growth
               'GO:0008285', # Negative regulation of cell proliferation
               'GO:0043069', # Negative regualtion of programmed cell death
               'GO:0012501', # programmed cell death
               'GO:0043067', # regulation of programmed cell death
               'GO:0001302', # Replicative cell aging           
               'GO:0090398', # Cellular senescence          
               'GO:0000723', # telomere maintenance         
               'GO:0032204', # regulation of telomere maintenance           
               'GO:2000772', # Regulation of cellular senescence            
               'GO:0032200', # telomere Organization
               'GO:0045766', # Positive regulation of angiogenesis
               'GO:0001570', # Vasculogenesis           
               'GO:0001525', # angiogenesis 
               'GO:0007162', # Negative regulation of cell adhesion         
               'GO:0001837', # Epithelial to mesenchymal transition         
               'GO:0016477', # Cell migration           
               'GO:0007155', # Cell adhesion            
               'GO:0034330', # Cell junction Organization           
               'GO:0030030', # Cell projection Organization
               'GO:0030155', # Regulation of Cell Adhesion
               'GO:0045005', # Maintenance of fidelity involved in DNA-dependent DNA replication            
               'GO:0006281', # DNA repair           
               'GO:0006282', # regulation of DNA repair         
               'GO:0031570', # DNA integrity Checkpoint
               'GO:0002367', # Cytokine production involved in immune response          
               'GO:0050727', # regulation of inflammatory response          
               'GO:0006096', # Glycolysis
               'GO:0071456', # Cellular response to hypoxia
               'GO:0002837', # Regulation of immune response to tumor cell          
               'GO:0002418', # Immune response to tumor cells           
               'GO:0050776', # Regulation of immune response            
               'GO:0006955') # Immune Response          

GO:0000165
go_annotations[go_annotations$GO %in% GO_cancer,]


#GSEA ANALYSIS
data("ExamplePathways")

fgseaResults <- limmaToFGSEA(limmaResults, gene.sets = ExamplePathways, 
                             min.set = 5, rank.by = "t",
                             skip.first = TRUE)




names(fgseaResults)



fgseaTab <- head(fgseaResults$Smoking.Habit.yes[
  order(fgseaResults$Smoking.Habit.yes$pval, 
        decreasing = FALSE),])
fgseaTab[,2:6] <- lapply(fgseaTab[,2:6], format, digits = 1, nsmall = 1)

knitr::kable(fgseaTab, 
             row.names = FALSE, format = "html", align = "c")



# Leading edge for pathways with adjusted p < 0.05
leading.edge <- fgseaToLEdge(fgsea.res = fgseaResults,
                             cutoff.type = "pval", cutoff = 0.05) 

pheatmap::pheatmap(t(leading.edge$Smoking.Habit.yes),
                   legend = FALSE,
                   color = c("white", "black"),
                   main = "p-value cutoff = 0.05")


# Leading edge for pathways with abs(NES) > 1
leading.edge.nes <- fgseaToLEdge(fgsea.res = fgseaResults,
                                 cutoff.type = "NES", cutoff = 1,
                                 nes.abs.cutoff = TRUE) 
pheatmap::pheatmap(t(leading.edge.nes$Smoking.Habit.yes),
                   legend = FALSE,
                   color = c("white", "black"),
                   main = "NES cutoff = 1")





