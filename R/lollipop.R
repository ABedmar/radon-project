#module load pcre2/10.35 gcc/9.2.0 R/4.2.2 bedtools
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/lolliplots/")
library(trackViewer)
library(randomcoloR)

genes<-c("Nf1","Nras","Met","Ctnnb1","Erbb3","Erbb4","Fgfr2","Keap1","Akt1","Fgfr3","Fgfr4","Mapk1","Pten","Alk","Braf","Egfr","Erbb2","Fgfr1","Ntrk3","Ret","Ros1","Stk11")
genes<-c("Met")


gannot=read.table("uniq_Rattus_norvegicus.mRatBN7.2.109.csv")
annot=read.table("rn6.exons.bed",sep='\t')
samples=read.table("AA_changes_EXTENDED_drivers.txt",sep='\t',h=T)

#samples[samples$Hugo_Symbol %in% "Erbb2",]


#pdf("drivers_plots.pdf")
for (gene in genes) 
  {
  for (transcript in gannot$V2[gannot$V3 == gene])
    {
    exons=annot[grep(transcript,annot$V4),]
    gs=samples[samples$Hugo_Symbol == gene,]
    
    if (nrow(exons != 0))
      {

      SNP <- gs$Start_Position - exons$V2
      
      sample.gr <- GRanges(exons$V1, IRanges(SNP, width=1, names=paste0(gs$Protein_Change, " ", gsub("\\.maf", "",gs$Source_MAF))))

      histology_colors <- c("#ef40d0", "#ef8240", "#40efc8") # Define colors for each histology
      sample.gr$border <- histology_colors[as.factor(gs$Histology)]
      
      group_colors <- c(D1 = "#440154", D3 = "#3b528b", D6 = "#21908c", D12 = "#5dc863", Pr = "#fde725")# Define colors for each group
      sample.gr$color <- group_colors[as.character(gs$Group)]
      
      
      
      sample.gr.rot <- sample.gr
      sample.gr.rot$label.parameter.rot <- 80
      
      
      features <- GRanges(exons$V1, 
                          IRanges(as.numeric(unlist(strsplit(exons$V12,","))),
                                            width=as.numeric(unlist(strsplit(exons$V11,","))),
                                            names=paste0("block", 1:length(as.numeric(unlist(strsplit(exons$V12,",")))))))
      
      palette <- randomColor(exons$V10)
      
      features$color <- palette
      features$fill <- palette
      
      
      svg(paste(gene,transcript,'lolliplot.svg',sep="."), height = 7, width = 20)
      lolliplot(sample.gr.rot, features, main=gene, ylab="")
      grid.text(gene, x=.5, y=.6, just="top", gp=gpar(cex=1.5, fontface="bold"))
      dev.off()
      
    }
  }
}
#dev.off()



help(lolliplot)
