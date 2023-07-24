setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cancer_genes_filtered")

library(maftools)
library(dplyr)
library(stringr)
library(ggplot2)
library(plotly)
library(forcats)

load("maf.Rdata")

data<-read.csv(file="radon-rats-cancer-genes.CSV",sep=",",header=T)

ggplot(data=data)+
  geom_col(aes(x=reorder(Tumor_Sample_Barcode, -total_perMB), y=total_perMB, fill=group))+
  theme_minimal()+
  geom_hline(yintercept=median(data$total_perMB), color="red", linetype="dashed")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(element_blank())+
  ylab("Total (per MB)")

means<-aggregate(data$total_perMB, list(data$group), FUN=mean)
means[2]<-round(means$x, digits = 2)
N <- aggregate(data$Tumor_Sample_Barcode, by=list(data$group), FUN=length)

ggplot(data=data)+
  geom_boxplot(aes(x=fct_inorder(group), y=total_perMB,color=group))+
  geom_text(data = means, aes(x=Group.1, label =  paste("mean:", x), y = 7.5))+
  geom_text(data = N, aes(x=Group.1, label =  paste("n =", x), y = 8))+
  geom_jitter(aes(x=group, y=total_perMB))+
  theme_minimal()+   
  geom_hline(yintercept=median(data$total_perMB), color="red", linetype="dashed")+ 
  theme(legend.position = "none")+
  ggtitle("TMB radon-rats (by group)")+
  xlab(element_blank())+
  ylab("Total (per MB)")


cancer_genes<-c("Ada", "Ahr", "Akt1", "Akt2", "Alb", "Aldoa", "Alox5", "Amd1", "Anxa1", "Anxa5", "Apc", "Apex1", "Apoa1", "Apoc3", "Apoe", "Aqp4", "Araf", "Areg", "Arg1", "Atf3", "Avp", "Azgp1", "B2m", "Bax", "Bcl2", "Bmp2", "Brca1", "Bsg", "Ddr1", "Casp3", "Cat", "Runx2", "Runx1", "Ccnb1", "Ccne1", "Ccng1", "Cd44", "Cd74", "Cd9", "Cdkn2a", "Cdkn2b", "Cftr", "Chga", "Chrna3", "Chrna5", "Chrna7", "Chrnb4", "Ckb", "Abcc2", "Col10a1", "Crebbp", "Csf1r", "Csf3", "Cycs", "Cyp1a1", "Cyp1b1", "Cyp24a1", "Cyp2a1", "Cyp2e1", 
         "Ace", "Nqo1", "Eef2k", "Egf", "Egfr", "Egr1", "Eno1", "Eno2", "Ephx1", "Erbb2", "Esr1", "Fgf18", "Fgf9", "Fgfr2", "Fgfr4", "Fh", "Fosl2", "Xrcc6", "Gapdh", "Gata1", "Gipr", "Gnas", "Gpi", "Gpx1", "Grik2", "Nr3c1", "Gstm1", "Gstm2", "Gstp1", "Gstt1", "H19", "Hk1", "Hk2", "Foxa1", "Foxa2", "Hp", "Hras", "Icam1", "Id2", "Idh1", "Igf2r", "Il10", "Il18", "Il1a", "Il1b", "Il6r", "Il7", "Irf1", "Itgb1", "Itgb4", "Jak2", "Junb", "Kdr", "Kras", "Ldha", "Lep", "Lta", "Smad3", "Smad4", "Slc3a2", "Met", "Mmp11", 
         "Mmp7", "Abcc1", "Muc1", "Myc", "Ncl", "Nf1", "Nf2", "Nfkbia", "Notch1", "Npm1", "Nras", "Ntrk3", "Oxtr", "Pcbd1", "Pdgfra", "Pgam1", "Abcb1b", "Serpina1", "Pik3r1", "Plau", "Plcg1", "Plod2", "Ppia", "Ppp3ca", "Prkaa1", "Prkcb", "Ptgis", "Pthlh", "Ptk2", "Ptpn11", "Ptprj", "Raf1", "Rarb", "Rasa1", "Rb1", "Ros1", "Sdc1", "Sdc4", "Sftpc", "Sftpd", "Shox2", "Slc19a1", "Sod2", "Serpina3n", "Spp1", "Sult1a1", "Stat3", "Stat5b", "Syp", "Tf", "Tgfa", "Timp3", "Tnf", "Faslg", "Tp53", "Tpi1", "Ttr", "Tyms", 
         "Uchl1", "Ugt1a1", "Vav1", "Wt1", "Xrcc5", "Yy1", "Hspb1", "Prkn", "Foxm1", "Casp9", "Map2k2", "Pecam1", "Hif1a", "Robo1", "Slc15a2", "Eef2", "Phgdh", "Pten", "Top2a", "Ptges", "Hes1", "Akt3", "Mbl2", "Akr1a1", "Mtor", "Ccnd1", "Epas1", "Slc4a4", "Jak1", "Cdkn1b", "Gpx3", "Hsd17b10", "Hnrnpa1", "Htra1", "Ager", "Axin2", "Ahcy", "Mre11", "Socs1", "Cdh1", "Hmgb2", "Slco1b2", "Slit2", "Erbb3", "Cdkn1a", "Lgals1", "Rrad", "Acsl4", "Ifi27", "Ccnh", "Ctnnb1", "Nat1", "Nat2", "Tert", "Map2k1", "Mapk14", 
         "Nfkb1", "Birc5", "Acsl3", "Bhlhe41", "Cd164", "Epha7", "Casp6", "Smad9", "Ascl1", "Rpsa", "Grik3", "Hnrnpl", "Bag6", "Csrp3", "Cdh13", "Rac1", "L1cam", "Erp29", "Nek2", "Atrx", "Serpina5", "Xrcc1", "Fas", "Ddr2", "Grm8", "Gclc", "Braf", "Rhoa", "Ccl22", "Abcb1a", "Hdac3", "Mat2a", "Vegfa", "Rims2", "Axl", "Ep300", "Prdx1", "Rxrg", "Chd4", "Runx3", "Atp5pd", "Bad", "Ncoa2", "Ncoa3", "Ncoa6", "Nnat", "Eif2b4", "Eif4ebp1", "Calr", "Pdpk1", "Ptgs2", "Nfe2l2", "Rock1", "Fhit", "Lmna", "Cxcr4", "Sox9", 
         "Nop58", "Erbb4", "Rpl10a", "Xpo1", "Lmnb1", "Mcl1", "Cblb", "Kit", "Apaf1", "Cux1", "Slc7a5", "Fasn", "Strbp", "Birc2", "Fgfr1", "Fgfr3", "Uck2", "Notch3", "Gnaq", "Msh2", "Src", "Srm", "Eif2b1", "Epha5", "Axin1", "Tp63", "Chi3l1", "Slc34a2", "Ppp2r1a", "Pik3ca", "Mlh1", "Casp7", "Casp8", "Map3k1", "Map3k8", "Dnmt1", "Dab2", "Tk1", "Lamc2", "Tnc", "Bard1", "Ccnd2", "Max", "Cdk12", "Cdk4", "Cdk6", "Il1r2", "Il1rn", "Mif", "Ogg1", "Slc4a7", "Serpina10", "Daxx", "Ubr5", "Tnfrsf1a", "Adamts1", "Mt3", 
         "Fap", "Fat1", "Trpc4", "Naip6", "Birc3", "Thbd", "Snrpb", "Mmp2", "Mmp3", "Mmp9", "Muc4", "E2f5", "Kpna2", "Ptch1", "Kcnj4", "Klf5", "Mef2d", "Bzw2", "Gng11", "Chrna9", "Sdha", "Kif5b", "Nrp1", "Keap1", "Stk39", "Dab2ip", "Rassf5", "Sftpb", "Smarca4", "Timm10", "Timm8a1", "Khsrp", "Bambi", "Pold1", "Rheb", "Elavl2", "Alk", "Pycard", "Rpl13a", "Tnfsf10", "Lama3", "Cst6", "Trpc7", "Itgb3", "Selenos", "Tlr9", "Rbm10", "Serpinf1", "E2f6", "Gsc", "Bbc3", "Eif2b5", "Ppp6c", "Serpina12", "Arid1b", "Slc1a5", 
         "Serpina4", "Hs3st2", "Trmt11", "Adam28", "Cox8c", "Gpr135", "Rela", "Cdkn1c", "Ddx24", "Fbxo11", "Gjb2", "E2f1", "Cbx7", "Myd88", "Tmem132d", "Impdh2", "Srr", "Cct5", "Psat1", "Snap47", "Acvr1b", "Ehmt2", "Osmr", "Msh5", "Tomm40", "Kif5a", "Dnmt3b", "Slc22a18", "Dnmt3a", "Slc22a22", "Eml4", "Phlda2", "Palb2", "Ppp2r1b", "Serpina9", "Ctdspl", "Smarcc1", "Pfas", "Rad52", "Suclg2", "Pabpc4", "Lrrc59", "Kif21a", "Rnf43", "Ipo4", "Abce1", "Rcl1", "Map4k4", "Alox12b", "Mdm2", "Lect2", "Kdm5a", "Ltf", 
         "Mki67", "Runx1t1", "Fbl", "Anapc2", "Satb1", "Setd2", "Pagr1", "Tymp", "Pxn", "Xpc", "Nudt21", "Atr", "Ptprb", "Mbd1", "Lrp1b", "Fkbp11", "Btbd7", "Ppan", "Chfr", "Nod2", "Kdm4a", "Asb2", "Car9", "Car12", "Cbx5", "Deup1", "Rpl36a", "Dot1l", "Otub2", "Traf4", "Rtel1", "Stim1", "Akr1c1", "Ercc1", "Ikzf3", "Map2k4", "Cntnap2", "Tp73", "Prr13", "Ercc3", "Epha2", "Rictor", "Kat2a", "Birc6", "Rnaseh2a", "Manf", "Ndrg1", "Prex2", "Wee1", "Clptm1l", "Serpina1f", "Mmp1", "Nuf2", "Nsd2", "Rbm7", "Tp53bp1", 
         "Col4a2", "E2f8", "Ergic3", "Mrps5", "Cyld", "Setdb1", "Sfxn1", "Arhgap5", "Crkl", "Shmt2", "Stk11", "Foxp1", "Ikzf2", "Ipo5", "Gart", "Plbd1", "Arhgap35", "Rassf8", "Dvl3", "Reck", "Irak3", "Nop56", "Meox1", "Col22a1", "Ercc2", "Cert1", "Foxo3", "Adcy1", "Dicer1", "Cmtr2", "Kdm2a", "Sox4", "Dok1", "Tab2", "Dll4", "Rpl27a", "Il7r", "Hells", "Brd7", "Rnaset2", "Bcl11a", "Aven", "Tmem135", "Mthfr", "Foxl2", "Ccdc116", "Cmpk1", "Rfc4", "Ndc80", "Alx4", "Xrn2", "E2f7", "Hnrnpa2b1", "Hsp90b1", "Arid1a", 
         "Itgav", "Grwd1", "Diablo", "Sos1", "Dok2", "Mpzl3", "Cadm1", "Recql", "Cda", "Arid2", "Strap", "Unc79", "Itga9", "Dnai7", "Frzb", "Pign", "Tent5a", "Rrs1", "Col7a1", "Aldh18a1", "Ercc6", "Eif3e", "Xpo5", "Spop", "Tet2", "Dapk1", "Lrrc56", "Plscr4", "Capg", "Prdm1", "Cdc73", "Rptor", "Arpc5", "Dok3", "Bap1", "Ms4a1", "Shmt1", "Ubr7", "Serpina11", "Srsf6", "Rassf1", "Srsf2", "Sdhc", "Dagla", "Ctnna1", "Wrap53", "Btk", "Cdyl", "Mdc1", "Ifi27l2b", "Hoxb7", "Zdbf2", "Map2k7", "Dusp3", "Dip2c", "Zfhx3", 
         "Ppp4r4", "Pacrg", "Paox", "Cbl", "E2f3", "Tnik", "Lyrm9", "Zfp462", "Macir", "Satb2", "Irs4", "Skp2", "Pml", "Ncapg", "Lamb3", "Pot1", "Kdm1a", "Ikzf1", "Gnpnat1", "Rangrf", "Ptpn13", "Cenpa", "Msh3", "Gar1", "Polr1d", "Fes", "Coro1c", "Scgb3a2", "Emx2", "Naa10", "Npm3", "Pbrm1", "H1f10", "Sox2", "Ttf1", "Ercc6l", "Bcorl1", "Kdm3b", "Prima1", "U2af1", "Snrpg", "Bend4", "Abl1", "Rnf168", "Mbd4", "Med12", "Kmt2a", "Ercc5", "Sox30", "Kmt2c", "Srsf1", "Mrpl12", "Kdm4b", "Mycl", "Golm1", "Fubp1", "Mpo", 
         "Rrm1", "Atm", "E2f2", "Rad21", "Pole", "Ezh2", "Serpina6", "Iqsec1", "Ncaph", "Prame", "Idh2", "Rrm2", "Mir210", "Prdm14", "Hilpda", "Fbxw7", "Snrpe", "Msh6", "Kmt2d", "Mir484", "Mir377", "Mir423", "Mir497", "Mir205", "Mir146a", "Mir136", "Mir130a", "Mir126b", "Mir31", "Mir29a", "Mir20a", "Mir346", "Mir144", "Mir127", "Mir665", "Mir431", "Mir19a", "Mir9-3", "Mir345", "Mir671", "Mir410", "Mir381", "Mir223", "Mir222", "Mir154", "Mir152", "Mir145", "Mir122", "Mir107", "Mir98", "Mir30c1", "Mir29b1", 
         "Mir21", "Mir337", "Mir193a", "Insm1", "Fam181a", "Kdm6a", "Snord22", "Mapk3", "Plcg2", "Prkca", "Rxra", "Rxrb", "Pik3r2", "Mapk1", "Grb2", "Pik3cb", "Pik3r3", "Pik3cd", "Pik3r5")

unique(laml@data$Hugo_Symbol)

df<-data[data$Hugo_Symbol %in% cancer_genes, ]

dim(df)
nrow(df[duplicated(df$Hugo_Symbol), ])


library(tibble)
library(dplyr)

top50<-tibble(df$Hugo_Symbol) %>% 
  group_by(df$Hugo_Symbol) %>% 
  count(sort = TRUE) %>%
  ungroup() %>% 
  slice_max(n, n=50)

write.table(top50,file = "TOP50.txt")

top50<-tibble(df$Variant_Classification) %>% 
  group_by(df$Hugo_Symbol) %>% 
  count(sort = TRUE) %>%
  ungroup() %>% 
  slice_max(n, n=50)

head(top50)
allmaf<-df[,c(2,10)]
top<-melt(allmaf, id.vars = unique(c("Hugo_Symbol")))


factor(top30$`df$Hugo_Symbol`)

#MATRIX OF ALL NUMBER OF MUTATIONS
all<-tibble(df$Hugo_Symbol) %>% 
  group_by(df$Hugo_Symbol) %>% 
  count(sort = TRUE) %>%
  ungroup() %>% 
  slice_max(n, n=nrow(df))


ggplot(data=top50,aes(x=reorder(factor(`df$Hugo_Symbol`), n),y=n, color=df$Source_MAF))+
  geom_bar(position = "stack",stat='identity')+
  coord_flip()+
  theme_minimal()+
  theme(legend.position = "none",panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 75, by = 5))+
  geom_hline(yintercept=median(all$n), color="red", linetype="dashed")+
  geom_hline(yintercept=mean<-mean(all$n), color="yellow", linetype="dotted")+
  ggtitle("Top 50 mutated genes in rat lung cancer related genes")+
  ylab("# mutations")+
  xlab("genes")





#maftools


vcs = getSampleSummary(laml)
vcs = vcs[,colnames(vcs)[!colnames(x = vcs) %in% c('total', 'Amp', 'Del', 'CNV_total')], with = FALSE]

vcs = vcs[,c(1,order(colSums(x = vcs[,2:(ncol(vcs)), with =FALSE]), decreasing = TRUE)+1), with =FALSE] #order based on most event
vcs.m = data.table::melt(data = vcs, id = 'Tumor_Sample_Barcode')
colnames(vcs.m) = c('Tumor_Sample_Barcode', 'Variant_Classification', 'N')

data.table::setDF(vcs)
rownames(x = vcs) = vcs$Tumor_Sample_Barcode
vcs = vcs[,-1]
vcs = t(vcs)


boxH = vcs.m[,boxplot.stats(N)$stat[5], by = .(Variant_Classification)]
colnames(boxH)[ncol(boxH)] = 'boxStat'
b = boxplot(N ~ Variant_Classification, data = vcs.m, col = col[levels(vcs.m$Variant_Classification)],
            axes = FALSE, outline = FALSE, lwd = 1, border = grDevices::adjustcolor(col = "black", alpha.f = 0.6))
axis(side = 2, at = as.integer(seq(0, max(boxH[,boxStat], na.rm = TRUE), length.out = 4)),
     lwd = 2, font = 2, cex.axis = fs, las = 2)
title(main = "Variant Classification summary", adj = 0, cex.main = fs, font = 2, line = 1)
b


























##########################################################################
df <- laml@data


# group the data by gene and count the number of occurrences
gene_count <- df %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>%   # sort by descending order of counts
  head(30)                   # select the top 30 genes

sort(gene_count$Hugo_Symbol)


genes <- df %>%
  group_by(Hugo_Symbol, Variant_Classification) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  arrange(Hugo_Symbol, desc(n)) %>%
  group_by(Hugo_Symbol)

filtered_genes <- genes %>%
  filter(Hugo_Symbol %in% gene_count$Hugo_Symbol)

ggplot(data=filtered_genes, aes(x=n, y=reorder(Hugo_Symbol, n, FUN=sum), fill=Variant_Classification))+
  geom_bar(position="stack", stat="identity")+
  theme_minimal()+
  xlab("Count")+
  ylab("Hugo Symbol")+
  ggtitle("Top Recurrent Cancer Genes in LAML")


# create a bar plot ################
all<-ggplot(data = gene_count, aes(x = count, y = reorder(Hugo_Symbol, count))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Count") +
  ylab("Gene") +
  theme_minimal()+
  ggtitle("Top Recurrent Cancer Genes in LAML") +
  theme(plot.title = element_text(hjust = 0.5))

################################################################
df_RnD6 <- df[grepl("^RnD6", df$Source_MAF), ]

genes_D6 <- df_RnD6 %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>%   # sort by descending order of counts
  head(30)

gene_count_D6 <- df_RnD6 %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>%   # sort by descending order of counts
  head(30)                   # select the top 30 genes

filtered_genes_D6 <- genes_D6 %>%
  filter(Hugo_Symbol %in% gene_count_D6$Hugo_Symbol)

ggplot(data=filtered_genes_D6, aes(x=n, y=reorder(Hugo_Symbol, n, FUN=sum), fill=Variant_Classification))+
  geom_bar(position="stack", stat="identity")+
  theme_minimal()+
  xlab("Count")+
  ylab("Hugo Symbol")+
  ggtitle("Top Recurrent Cancer Genes in LAML RnD6")




# create a bar plot
D6<-ggplot(data = gene_count_D6, aes(x = count, y = reorder(Hugo_Symbol, count))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Count") +
  ylab("Gene") +
  theme_minimal()+
  ggtitle("Top Recurrent Cancer Genes in low comulative") +
  theme(plot.title = element_text(hjust = 0.5))





################################################################
df_RnD3 <- df[grepl("^RnD3", df$Source_MAF), ]

gene_count_D3 <- df_RnD3 %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>%   # sort by descending order of counts
  head(30)  

# create a bar plot
D3<-ggplot(data = gene_count_D3, aes(x = count, y = reorder(Hugo_Symbol, count))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Count") +
  ylab("Gene") +
  theme_minimal()+
  ggtitle("Top Recurrent Cancer Genes in intermediate high") +
  theme(plot.title = element_text(hjust = 0.5))

################################################################
df_RnD12 <- df[grepl("^RnD12", df$Source_MAF), ]

gene_count_D12 <- df_RnD12 %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>%   # sort by descending order of counts
  head(30)  

# create a bar plot
D12<-ggplot(data = gene_count_D12, aes(x = count, y = reorder(Hugo_Symbol, count))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Count") +
  ylab("Gene") +
  theme_minimal()+
  ggtitle("Top Recurrent Cancer Genes in intermediate low") +
  theme(plot.title = element_text(hjust = 0.5))

################################################################
df_high <- df[grepl("^RnD1-|^RnPr", df$Source_MAF), ]

gene_count_high <- df_high %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>%   # sort by descending order of counts
  head(30)  

# create a bar plot
D1_Pr<-ggplot(data = gene_count_high, aes(x = count, y = reorder(Hugo_Symbol, count))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Count") +
  ylab("Gene") +
  theme_minimal()+
  ggtitle("Top Recurrent Cancer Genes in high PAEC") +
  theme(plot.title = element_text(hjust = 0.5))

################################################################
df_high_comulative <- df[grepl("^RnD1-|^RnPr|^RnD3|^RnD12-", df$Source_MAF), ]

gene_count_high_comulative <- df_high_comulative %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>%   # sort by descending order of counts
  head(30)  

# create a bar plot
D1_Pr_D12_D3<-ggplot(data = gene_count_high_comulative, aes(x = count, y = reorder(Hugo_Symbol, count))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Count") +
  ylab("Gene") +
  theme_minimal()+
  ggtitle("Top Recurrent Cancer Genes in high comulative") +
  theme(plot.title = element_text(hjust = 0.5))


################################################################
df_low_PAEC <- df[grepl("^RnD3-|^Rn6|^RnD12-", df$Source_MAF), ]

gene_count_low_PAEC <- df_low_PAEC %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>%   # sort by descending order of counts
  head(30)  

# create a bar plot
D6_D3_D12<-ggplot(data = gene_count_low_PAEC, aes(x = count, y = reorder(Hugo_Symbol, count))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Count") +
  ylab("Gene") +
  theme_minimal()+
  ggtitle("Top Recurrent Cancer Genes in low PAEC") +
  theme(plot.title = element_text(hjust = 0.5))


################################################################
##########################################################################
#BY COMULATIVE DOSE:
grid.arrange(D1_Pr_D12_D3, D6, ncol=2)

#BY DOSE RATE:
grid.arrange(D1_Pr, D6_D3_D12, ncol=2)


#BY BOTH:
grid.arrange(D1_Pr, D3, D12, D6, ncol=2)
