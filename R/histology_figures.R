library(extrafont)
font_import()
loadfonts(device="win")       #Register fonts for Windows bitmap output
fonts() 


#HISTOLOGY

# Dataframe 1: Location
Location <- data.frame(
  Location = c("Right", "Left"),
  `Overall_(N=36)` = c(67, 33),
  `RnD1_(N=12)` = c(67, 33),
  `RnD3_(N=7)` = c(86, 14),
  `RnD6_(N=4)` = c(75, 25),
  `RnD12_(N=9)` = c(56, 44),
  `RnPr_(N=4)` = c(50, 50)
)

# Dataframe 2: Histology
Histology <- data.frame(
  Histology = c("Squamous", "Adenocarcinoma", "Mixed"),
  `Overall_(N=36)` = c(19, 75, 6),
  `RnD1_(N=12)` = c(42, 58, 0),
  `RnD3_(N=7)` = c(0, 100, 0),
  `RnD6_(N=4)` = c(0, 100, 0),
  `RnD12_(N=9)` = c(22, 67, 11),
  `RnPr_(N=4)` = c(0, 75, 25)
)

# Dataframe 3: Adenocarcinoma
Adenocarcinoma <- data.frame(
  Adenocarcinoma = c("Acinar", "Papillary", "Solid", "Not specified"),
  `Overall_(N=28)` = c(28, 54, 7, 11),
  `RnD1_(N=7)` = c(43, 43, 14, 0),
  `RnD3_(N=7)` = c(43, 43, 14, 0),
  `RnD6_(N=4)` = c(25, 75, 0, 0),
  `RnD12_(N=6)` = c(0, 83, 0, 17),
  `RnPr_(N=4)` = c(25, 25, 0, 50)
)

# Dataframe 4: Stage
Stage <- data.frame(
  Stage = c("T1", "T2", "T3", "T4"),
  `Overall_(N=36)` = c(3, 3, 39, 55),
  `RnD1_(N=12)` = c(0, 0, 25, 75),
  `RnD3_(N=7)` = c(0, 14, 57, 29),
  `RnD6_(N=4)` = c(25, 0, 25, 50),
  `RnD12_(N=9)` = c(0, 0, 44, 56),
  `RnPr_(N=4)` = c(0, 0, 50, 50)
)


# Dataframe 5: Stage - N
Stage_N <- data.frame(
  `Stage_N` = c("N0", "N1"),
  `Overall_(N=36)` = c(89, 11),
  `RnD1_(N=12)` = c(75, 25),
  `RnD3_(N=7)` = c(100, 0),
  `RnD6_(N=4)` = c(100, 0),
  `RnD12_(N=9)` = c(100, 0),
  `RnPr_(N=4)` = c(75, 25)
)


# Dataframe 6: Stage - Pleural inv.
Stage_Pleural <- data.frame(
  `Stage_Pleural` = c("Pleural invasion", "No"),
  `Overall_(N=36)` = c(83, 17),
  `RnD1_(N=12)` = c(100, 0),
  `RnD3_(N=7)` = c(57, 43),
  `RnD6_(N=4)` = c(100, 0),
  `RnD12_(N=9)` = c(78, 22),
  `RnPr_(N=4)` = c(75, 25)
)



melted_df1 <- melt(Location, id.vars = "Location")
melted_df2 <- melt(Histology, id.vars = "Histology")
melted_df3 <- melt(Adenocarcinoma, id.vars = "Adenocarcinoma")
melted_df4 <- melt(Stage, id.vars = "Stage")
melted_df5 <- melt(Stage_N, id.vars = "Stage_N")
melted_df6 <- melt(Stage_Pleural, id.vars = "Stage_Pleural")


A<-ggplot(melted_df1, aes(x =variable , y =value , fill =Location, label = value))+
  geom_bar(stat="identity", position = "stack", width = 0.4)+
  geom_text(data=subset(melted_df1, value != 0), aes(y = value, label = value, family="Calibri"), size = 5, position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c("#9E8AE3", "#E0D358"))+
  scale_x_discrete(labels=c("Overall (N=36)", "RnD1 (N=12)", "RnD3 (N=7)", "RnD6 (N=4)", "RnD12 (N=9)", "RnPr (N=4)"))+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank(), legend.position = "top", legend.title = element_blank())+
  xlab("")+
  ylab("")+
  labs(title = "A")+
  theme(text=element_text(family="Calibri", size = 17), plot.title = element_text(size = 18, face = "bold"))


B<-ggplot(melted_df2, aes(x =variable , y =value , fill =Histology, label = value ))+
  geom_bar(stat="identity", position = "stack", width = 0.4)+
  geom_text(data=subset(melted_df2, value != 0), aes(y = value, label = value, family="Calibri"), size = 5, position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c("#E0D358", "#9E8AE3", "#D675E0"))+
  scale_x_discrete(labels=c("Overall (N=36)", "RnD1 (N=12)", "RnD3 (N=7)", "RnD6 (N=4)", "RnD12 (N=9)", "RnPr (N=4)"))+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank(), legend.position = "top", legend.title = element_blank())+
  xlab("")+
  ylab("")+
  labs(title = "B")+
  theme(text=element_text(family="Calibri", size = 17), plot.title = element_text(size = 18, face = "bold"))


C<-ggplot(melted_df3, aes(x =variable , y =value , fill =Adenocarcinoma, label = value ))+
  geom_bar(stat="identity", position = "stack", width = 0.4)+
  geom_text(data=subset(melted_df3, value != 0), aes(y = value, label = value, family="Calibri"), size = 5, position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c("#9E8AE3", "#B3B3B3", "#E0D358", "#D675E0"))+
  scale_x_discrete(labels=c("Overall (N=28)", "RnD1 (N=7)", "RnD3 (N=7)", "RnD6 (N=4)", "RnD12 (N=6)", "RnPr (N=4)"))+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank(), legend.position = "top", legend.title = element_blank())+
  xlab("")+
  ylab("")+
  labs(title = "C")+
  theme(text=element_text(family="Calibri", size = 17), plot.title = element_text(size = 18, face = "bold"))






D<-ggplot(melted_df4, aes(x =variable , y =value , fill =Stage, label = value ))+
  geom_bar(stat="identity", position = "stack", width = 0.4)+
  geom_text(data=subset(melted_df4, value > 3 ), aes(y = value, label = value, family="Calibri"), size = 5, position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c("#d5ffff", "#9ce0db", "#5dc1b9", "#338b85"))+
  scale_x_discrete(labels=c("Overall (N=36)", "RnD1 (N=12)", "RnD3 (N=7)", "RnD6 (N=4)", "RnD12 (N=9)", "RnPr (N=4)"))+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank(), legend.position = "top", legend.title = element_blank())+
  xlab("")+
  ylab("")+
  labs(title = "D")+
  theme(text=element_text(family="Calibri", size = 17), plot.title = element_text(size = 18, face = "bold"))


E<-ggplot(melted_df5, aes(x =variable , y =value , fill =Stage_N, label = value ))+
  geom_bar(stat="identity", position = "stack", width = 0.4)+
  geom_text(data=subset(melted_df5, value != 0), aes(y = value, label = value, family="Calibri"), size = 5, position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c("#9ce0db", "#338b85"))+
  scale_x_discrete(labels=c("Overall (N=36)", "RnD1 (N=12)", "RnD3 (N=7)", "RnD6 (N=4)", "RnD12 (N=9)", "RnPr (N=4)"))+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank(), legend.position = "top", legend.title = element_blank())+
  xlab("")+
  ylab("")+
  labs(title = "E")+
  theme(text=element_text(family="Calibri", size = 17), plot.title = element_text(size = 18, face = "bold"))


F_<-ggplot(melted_df6, aes(x =variable , y =value , fill =Stage_Pleural, label = value ))+
  geom_bar(stat="identity", position = "stack", width = 0.4)+
  geom_text(data=subset(melted_df6, value != 0), aes(y = value, label = value, family="Calibri"), size = 5, position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c("#9ce0db", "#338b85"))+
  scale_x_discrete(labels=c("Overall (N=36)", "RnD1 (N=12)", "RnD3 (N=7)", "RnD6 (N=4)", "RnD12 (N=9)", "RnPr (N=4)"))+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank(), legend.position = "top", legend.title = element_blank())+
  xlab("")+
  ylab("")+
  labs(title = "F")+
  theme(text=element_text(family="Calibri", size = 17), plot.title = element_text(size = 18, face = "bold"))



grid.arrange(A, B, C, nrow=3)
grid.arrange(D, E, F_, nrow=3)









# legend changes 
ua <- data.frame(
  `x` = c("C>T", "T>C", "C>A", "T>G", "C>G", "T>A"),
  `Overall_(N=36)` = c(60, 20, 10, 5, 3, 2)
)
ua$x <- factor(ua$x, levels = c("C>T", "T>C", "C>A", "T>G", "C>G", "T>A"))

melted_df7 <- melt(ua, id.vars = "x")


ggplot(melted_df7, aes(x =variable , y =value , fill =x))+
  geom_bar(stat="identity", position = "stack", width = 0.4)+
  scale_fill_manual(values = c("#f6695e","#ffcd39","#4dabf5","#ffad33","#6574c4","#70bf73"))

