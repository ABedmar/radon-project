setwd('/Users/alexb/OneDrive/Escritorio/IDIBAPS')

TMB <- read.csv(file = "Radon_TMB.csv", sep = ";", header = T)

rownames(TMB)<-TMB$Tumor_Sample_Barcode
head(TMB)

library(ggplot2)
library(plotly)

fig <-
  ggplot(TMB, aes(x=total_perMB , y=total, col=Group))+
  geom_point()+
  geom_text(
    label=rownames(TMB), 
    nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = T
  )+
  ggtitle("Radon Tumor Mutational Burden by sample")+
  xlab("")+
  ylab("TMB")+
  theme_minimal()

ggplotly(fig)


TMB_barplot_by_group <- ggplot(TMB, aes(x=Tumor_Sample_Barcode , y=total, fill=Group))+
  geom_col()+
  ggtitle("Radon Tumor Mutational Burden by sample")+
  ylab("TMB")+
  xlab("")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
ggsave("TMB_barplot_by_group.png")
  


  
TMB_ordered_barplot <- ggplot(TMB, aes(x=reorder(Tumor_Sample_Barcode, -total) , y=total, fill=Group))+
    geom_col()+
    ggtitle("Radon Tumor Mutational Burden by sample")+
    ylab("TMB")+
    xlab("")+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
ggsave("TMB_ordered_barplot.png")

TMB_grouped_boxplot <- ggplot(TMB, aes(x=Group , y=total, fill=Group))+
  geom_boxplot(show.legend = FALSE)+
  ggtitle("Radon Tumor Mutational Burden by sample")+
  ylab("TMB")+
  xlab("")+
  theme_minimal()

ggsave("TMB_grouped_boxplot.png")

