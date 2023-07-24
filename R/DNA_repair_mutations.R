#DNA REPAIR PATHWAYS
library(gridExtra)
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1")

load("maf_v2.Rdata")

DBS_response<-c("Chek2","TP53","Brcc3","Rad50","Abraxas1","Eya4","Mre11","Kdm4b","Babam2","Eya3","Uimc1","Ppp5c","Apbb1","Rnf168","Tp53bp1","Rnf8","Pias4","Eya2","Nbn","Atm","Kdm4a","Herc2","Eya1","Ube2i","Babam1","Kpna2","Nsd2","Baz1b","Smarca5","Sumo1","Ube2n","Ube2v2","Hsaz","Rps27a","Ubc")
NHEJ<-c("Rnf168","Poll","Tp53bp1","Pias4","Lig4","Rnf8","Xrcc4","Dclre1c","Tdp1","Polm","Nhej1","Mre11","Rad50","Xrcc6","Nbn","Atm","Prkdc","Rif1","Herc2","Xrcc5","Tdp2","Nsd1","Sumo1","Pxip1","Ube2n","Ube2v2","H2ax")
ATM_signaling<-c("Nbn","Atm","Baz1b","Brca1","Brcc3l2","Eya1","Eya3","H2ax","Kat5","Kdm4a","L3mbtl1","Mdc1","Mre11","Otub1","Pias1","Pias4","Ppm1d","Ppp4c","Ppp6c","Psmd14","Pwwp3a","Rad50","Rbbp8","Rif1","Rnf168","Rnf20","Rnf4","Rnf40","Tp53bp1","Trip12","Ube2n","Ubr5","Usp16","Usp44","Vcp")
HRR<-c("Ubc","Topbp1","Rfc3","Rad50","Mre11","Wrm","Rfc1","Rpa1","Pole","Pold1","Rad9a","Rad51c","Xrcc3","Brip1","Brca2","Eme1","Hus1","Palb2","Exo1","Rad9b","Rpa2","Rpa3","Xrcc2","Slx1b","Rad51","Slc25a16","Pole4","Stt3a","Rps27a")
NER<-c("ERCC1 ","Pole2","Actb","Ddb1","Rad23b","Pcna","Rps27a","Cops2","Cul4a","Gps1","Cops4","Cops3","Polr2e","Porl2b","Ubc","cops5","Parp1","Rpa1","Polr2k","Prpf10","Cdk7","Cops8","Porl2g","Cetn2","Rfc1","Usp7","Rpa3","Rad23a","Hmgn1","Mcrs1","Isy1","Porl2c","Rfc4","Lig3","Pold4","Aqr","Cul4b","Polr2a","Xpc","Ube2v2","Ino80c","Polr2d","Tcea1","Rnf111","Ccnh","Ercc4","Ercc3","Ercc5","Xrcc1","Pcna","Lig1")
BER<-c("H2afb3","Neil2","Neil3","Pole2","Pcn1","H2az1","H2ax","Parp1","Rpa1","Rpa3","Apex1","Rfc1","Parg","Pole4","Tdg","Pold4","Rfc4","Lig3","Terf2ip","Adprhl2","Xrcc1","Rpa2","Rfc3","Terf1","Pnkp","Pold3","Pot1","Tinf2","Rfc5","Polb","Ung","Acd","Rfc2","Neil1","Parp2","Pold1","Mutyh","Pold2","Fen1","Nthl1","Pole","H2az2","Smug1","Terf2ip","Pole3","Mpg")
MMR<-c("Nbn ","Exo1","Pcna","Rpa1","Rpa3","Msh6","Pold4","Rpa2","Msh2","Pold3","Mlh1","Pms2","Pold1","Pold2")


DBS_response<-toupper(DBS_response)
NHEJ<-toupper(NHEJ)
ATM_signaling<-toupper(ATM_signaling)
HRR<-toupper(HRR)
NER<-toupper(NER)
BER<-toupper(BER)
MMR<-toupper(MMR)



dbs <- laml@data[laml@data$Hugo_Symbol %in% DBS_response]
nhej <- laml@data[laml@data$Hugo_Symbol %in% NHEJ]
atm <- laml@data[laml@data$Hugo_Symbol %in% ATM_signaling]
hrr <- laml@data[laml@data$Hugo_Symbol %in% HRR]
ner <- laml@data[laml@data$Hugo_Symbol %in% NER]
ber <- laml@data[laml@data$Hugo_Symbol %in% BER]
mmr <- laml@data[laml@data$Hugo_Symbol %in% MMR]


dbs<-table(dbs$Hugo_Symbol, dbs$Tumor_Sample_Barcode)
nhej<-table(nhej$Hugo_Symbol, nhej$Tumor_Sample_Barcode)
atm<-table(atm$Hugo_Symbol, atm$Tumor_Sample_Barcode)
hrr<-table(hrr$Hugo_Symbol, hrr$Tumor_Sample_Barcode)
ner<-table(ner$Hugo_Symbol, ner$Tumor_Sample_Barcode)
ber<-table(ber$Hugo_Symbol, ber$Tumor_Sample_Barcode)
mmr<-table(mmr$Hugo_Symbol, mmr$Tumor_Sample_Barcode)



# List of dataframe names
dataframe_names <- c("dbs", "nhej", "atm", "hrr", "ner", "ber", "mmr")

# Loop through the dataframe names
for (name in dataframe_names) {
  # Get the current dataframe
  current_df <- get(name)
  
  # Perform the operations on the current dataframe
  current_df_melt <- melt(current_df)
  current_df_melt$group <- substr(current_df_melt$Var2, 1, 3)
  
  
  current_df_melt <- current_df_melt %>%
    mutate(group = case_when(
      as.character(group) %in% substitution_guidelines$group_value ~
        substitution_guidelines$group_substitute[match(as.character(group), substitution_guidelines$group_value)],
      TRUE ~ as.character(group)
    ))
  
  
  #current_df_melt$group <- gsub("-", "", current_df_melt$group)
  #current_df_melt$group <- gsub("Pr1", "Pr", current_df_melt$group)
  #current_df_melt$group <- gsub("Pr2", "Pr", current_df_melt$group)
  #current_df_melt$group <- factor(current_df_melt$group, levels=c("D1", "D3", "D6", "D12", "Pr"))
  
  # Save the current dataframe as a new dataframe with a specific name
  assign(paste0(name, "_long"), current_df_melt)
}


long_dataframe_names <- c("dbs_long", "nhej_long", "atm_long", "hrr_long", "ner_long", "ber_long", "mmr_long")

# Loop through the dataframe names
for (name in long_dataframe_names) {
  # Get the current dataframe
  current_df <- get(name)
  assign(paste0(name, "_heatmap"), heatmap)
  
  heatmap<-ggplot(data = current_df, mapping = aes(x = Var2,y = reorder(Var1, value),fill = value)) +
    geom_tile()+
    xlab("")+
    ylab("")+
    scale_fill_gradient(name = "# mut.",
                        low = "#FFFFFF",
                        high = "#a52a2a")+
    facet_grid(~ group, switch = "x", scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          strip.placement = "outside")+
    ggtitle(label = "")
}



dbs_long_heatmap <- dbs_long_heatmap + ggtitle("DBS_response") + theme(legend.key.size = unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),  legend.title = element_text(size = 9))
nhej_long_heatmap <- nhej_long_heatmap + ggtitle("NHEJ") + theme(legend.key.size = unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),  legend.title = element_text(size = 9))
atm_long_heatmap <- atm_long_heatmap + ggtitle("ATM_signaling") + theme(legend.key.size = unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),  legend.title = element_text(size = 9))
hrr_long_heatmap <- hrr_long_heatmap + ggtitle("HRR") + theme(legend.key.size = unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),  legend.title = element_text(size = 9))
ner_long_heatmap <- ner_long_heatmap + ggtitle("NER") + theme(legend.key.size = unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),  legend.title = element_text(size = 9))
ber_long_heatmap <- ber_long_heatmap + ggtitle("BER") + theme(legend.key.size = unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),  legend.title = element_text(size = 9))
mmr_long_heatmap <- mmr_long_heatmap + ggtitle("MMR") + theme(legend.key.size = unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),  legend.title = element_text(size = 9))



grid.arrange(dbs_long_heatmap, nhej_long_heatmap, atm_long_heatmap, hrr_long_heatmap, ner_long_heatmap, ber_long_heatmap, mmr_long_heatmap, ncol = 4)



#OVERALL
OVERALL<-c("Chek2","TP53","Brcc3","Rad50","Abraxas1","Eya4","Mre11","Kdm4b","Babam2","Eya3","Uimc1","Ppp5c","Apbb1","Rnf168","Tp53bp1","Rnf8","Pias4","Eya2","Nbn","Atm","Kdm4a","Herc2","Eya1","Ube2i","Babam1","Kpna2","Nsd2","Baz1b","Smarca5","Sumo1","Ube2n","Ube2v2","Hsaz","Rps27a","Ubc ","Poll","Lig4","Xrcc4","Dclre1c","Tdp1","Polm","Nhej1","Xrcc6","Prkdc","Rif1","Xrcc5","Tdp2","Nsd1","Pxip1","H2ax","Brca1","Brcc3l2","Kat5","L3mbtl1","Mdc1","Otub1","Pias1","Ppm1d","Ppp4c","Ppp6c","Psmd14","Pwwp3a","Rbbp8","Rnf20","Rnf4","Rnf40","Trip12","Ubr5","Usp16","Usp44","Vcp","Ubc","Topbp1","Rfc3","Wrm","Rfc1","Rpa1","Pole","Pold1","Rad9a","Rad51c","Xrcc3","Brip1","Brca2","Eme1","Hus1","Palb2","Exo1","Rad9b","Rpa2","Rpa3","Xrcc2","Slx1b","Rad51","Slc25a16","Pole4","Stt3a","ERCC1 ","Pole2","Actb","Ddb1","Rad23b","Pcna","Cops2","Cul4a","Gps1","Cops4","Cops3","Polr2e","Porl2b","cops5","Parp1","Polr2k","Prpf10","Cdk7","Cops8","Porl2g","Cetn2","Usp7","Rad23a","Hmgn1","Mcrs1","Isy1","Porl2c","Rfc4","Lig3","Pold4","Aqr","Cul4b","Polr2a","Xpc","Ino80c","Polr2d","Tcea1","Rnf111","Ccnh","Ercc4","Ercc3","Ercc5","Xrcc1","Lig1","H2afb3","Neil2","Neil3","Pcn1","H2az1","Apex1","Parg","Tdg","Terf2ip","Adprhl2","Terf1","Pnkp","Pold3","Pot1","Tinf2","Rfc5","Polb","Ung","Acd","Rfc2","Neil1","Parp2","Mutyh","Pold2","Fen1","Nthl1","H2az2","Smug1","Pole3","Mpg","Nbn ","Msh6","Msh2","Mlh1","Pms2")
OVERALL<-toupper(OVERALL)

overall <- laml@data[laml@data$Hugo_Symbol %in% OVERALL]
overall<-table(overall$Hugo_Symbol, overall$Tumor_Sample_Barcode)

overall_melt<-melt(overall)
overall_melt$group <- substr(overall_melt$Var2, 1, 3)


#overall_melt$group<-gsub("-", "", overall_melt$group)
#overall_melt$group<-gsub("Pr1", "Pr", overall_melt$group)
#overall_melt$group<-gsub("Pr2", "Pr", overall_melt$group)
#overall_melt$group <- factor(overall_melt$group, levels=c("D1", "D3", "D6", "D12", "Pr"))



overall_long_heatmap<-ggplot(data = overall_melt, mapping = aes(x = Var2,y = reorder(Var1, value),fill = value)) +
  geom_tile()+
  xlab("")+
  ylab("")+
  scale_fill_gradient(name = "# mut.",
                      low = "#FFFFFF",
                      high = "#a52a2a")+
  facet_grid(~ group, switch = "x", scales = "free_x", space = "free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        strip.placement = "outside",
        legend.key.size = unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),  legend.title = element_text(size = 9))+
  ggtitle("DNA repair")

overall_long_heatmap










# MINERSS DATA
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS")

load("maf_pathogenic.Rdata")

library(dplyr)

# Define the guidelines for substitution
substitution_guidelines <- tribble(
  ~group_value, ~group_substitute,
  "108", "HIGH",
  "118", "LOW",
  "119", "HIGH",
  "166", "HIGH",
  "219", "HIGH",
  "269", "INTERMEDIATE",
  "278", "INTERMEDIATE",
  "292", "INTERMEDIATE",
  "300", "INTERMEDIATE",
  "309", "INTERMEDIATE",
  "340", "INTERMEDIATE",
  "350", "HIGH",
  "372", "HIGH",
  "374", "HIGH",
  "389", "LOW",
  "419", "LOW",
  "442", "HIGH",
  "443", "HIGH",
  "457", "HIGH",
  "520", "LOW",
  "585", "HIGH",
  "603", "INTERMEDIATE",
  "616", "HIGH",
  "631", "HIGH",
  "635", "INTERMEDIATE",
  "639", "HIGH"
)



substitution_guidelines <- tribble(
  ~group_value, ~group_substitute,
  "108", "SMOKER",
  "118", "SMOKER",
  "119", "NON-SMOKER",
  "166", "SMOKER",
  "219", "SMOKER",
  "269", "NON-SMOKER",
  "278", "SMOKER",
  "292", "SMOKER",
  "300", "SMOKER",
  "309", "NON-SMOKER",
  "340", "SMOKER",
  "350", "NON-SMOKER",
  "372", "SMOKER",
  "374", "SMOKER",
  "389", "NON-SMOKER",
  "419", "SMOKER",
  "442", "SMOKER",
  "443", "SMOKER",
  "457", "NON-SMOKER",
  "520", "SMOKER",
  "585", "SMOKER",
  "603", "SMOKER",
  "616", "SMOKER",
  "631", "SMOKER",
  "635", "NON-SMOKER",
  "639", "SMOKER"
)




# Perform the substitution
overall_melt <- overall_melt %>%
  mutate(group = case_when(
    as.character(group) %in% substitution_guidelines$group_value ~
      substitution_guidelines$group_substitute[match(as.character(group), substitution_guidelines$group_value)],
    TRUE ~ as.character(group)
  ))


laml@data$Tumor_Sample_Barcode



