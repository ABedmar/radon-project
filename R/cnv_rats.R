
library(circlize)


################################################################################
###### CNV
setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-RATS/cnv/cns/D1")
RnD1_36<-read.table(file = "RnD1-36.call.cns", sep="\t", head=T)
RnD1_36<-RnD1_36[,-c(4,6,7,8,9,10)]
colnames(RnD1_36)<-c("chr", "start", "end", "value1")
RnD1_36$chr <- paste0("chr", RnD1_36$chr)

RnD1_81<-read.table(file = "RnD1-81.call.cns", sep="\t", head=T)
RnD1_81<-RnD1_81[,-c(4,6,7,8,9,10)]
colnames(RnD1_81)<-c("chr", "start", "end", "value1")
RnD1_81$chr <- paste0("chr", RnD1_81$chr)

RnD1_98<-read.table(file = "RnD1-98.call.cns", sep="\t", head=T)
RnD1_98<-RnD1_98[,-c(4,6,7,8,9,10)]
colnames(RnD1_98)<-c("chr", "start", "end", "value1")
RnD1_98$chr <- paste0("chr", RnD1_98$chr)

RnD1_104<-read.table(file = "RnD1-104.call.cns", sep="\t", head=T)
RnD1_104<-RnD1_104[,-c(4,6,7,8,9,10)]
colnames(RnD1_104)<-c("chr", "start", "end", "value1")
RnD1_104$chr <- paste0("chr", RnD1_104$chr)


##############################################

RnD1_132<-read.table(file = "RnD1-132.call.cns", sep="\t", head=T)
RnD1_132<-RnD1_132[,-c(4,6,7,8,9,10)]
colnames(RnD1_132)<-c("chr", "start", "end", "value1")
RnD1_132$chr <- paste0("chr", RnD1_132$chr)

RnD1_134<-read.table(file = "RnD1-134.call.cns", sep="\t", head=T)
RnD1_134<-RnD1_134[,-c(4,6,7,8,9,10)]
colnames(RnD1_134)<-c("chr", "start", "end", "value1")
RnD1_134$chr <- paste0("chr", RnD1_134$chr)

RnD1_162<-read.table(file = "RnD1-162.call.cns", sep="\t", head=T)
RnD1_162<-RnD1_162[,-c(4,6,7,8,9,10)]
colnames(RnD1_162)<-c("chr", "start", "end", "value1")
RnD1_162$chr <- paste0("chr", RnD1_162$chr)

RnD1_171<-read.table(file = "RnD1-171.call.cns", sep="\t", head=T)
RnD1_171<-RnD1_171[,-c(4,6,7,8,9,10)]
colnames(RnD1_171)<-c("chr", "start", "end", "value1")
RnD1_171$chr <- paste0("chr", RnD1_171$chr)

##############################################

RnD1_189<-read.table(file = "RnD1-189.call.cns", sep="\t", head=T)
RnD1_189<-RnD1_189[,-c(4,6,7,8,9,10)]
colnames(RnD1_189)<-c("chr", "start", "end", "value1")
RnD1_189$chr <- paste0("chr", RnD1_189$chr)

RnD1_195<-read.table(file = "RnD1-195.call.cns", sep="\t", head=T)
RnD1_195<-RnD1_195[,-c(4,6,7,8,9,10)]
colnames(RnD1_195)<-c("chr", "start", "end", "value1")
RnD1_195$chr <- paste0("chr", RnD1_195$chr)

RnD1_207<-read.table(file = "RnD1-207.call.cns", sep="\t", head=T)
RnD1_207<-RnD1_207[,-c(4,6,7,8,9,10)]
colnames(RnD1_207)<-c("chr", "start", "end", "value1")
RnD1_207$chr <- paste0("chr", RnD1_207$chr)

RnD1_214<-read.table(file = "RnD1-214.call.cns", sep="\t", head=T)
RnD1_214<-RnD1_214[,-c(4,6,7,8,9,10)]
colnames(RnD1_214)<-c("chr", "start", "end", "value1")
RnD1_214$chr <- paste0("chr", RnD1_214$chr)

##############################################

# Create a circos plotGROUP D1
layout(matrix(1:3, 1, 3))
for(i in 1:3) {
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  circos.par(cell.padding = c(0, 0, 0, 0))
  circos.clear()
  circos.initializeWithIdeogram(species = "rn6")
  circos.genomicTrack(RnD1_36,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD1_81,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD1_98,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD1_104,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  ##############################################
  
  # Create a circos plot GROUP D1
  circos.clear()
  circos.initializeWithIdeogram(species = "rn6")
  circos.genomicTrack(RnD1_132,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD1_134,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD1_162,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD1_171,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  
  ##############################################
  
  circos.clear()
  circos.initializeWithIdeogram(species = "rn6")
  circos.genomicTrack(RnD1_189,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD1_195,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD1_207,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD1_214,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
}

##############################################
##############################################
##############################################
# Create a circos plot GROUP D3
setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-RATS/cnv/cns/D3")

RnD3_71<-read.table(file = "RnD3-71.call.cns", sep="\t", head=T)
RnD3_71<-RnD3_71[,-c(4,6,7,8,9,10)]
colnames(RnD3_71)<-c("chr", "start", "end", "value1")
RnD3_71$chr <- paste0("chr", RnD3_71$chr)

RnD3_92<-read.table(file = "RnD3-92.call.cns", sep="\t", head=T)
RnD3_92<-RnD3_92[,-c(4,6,7,8,9,10)]
colnames(RnD3_92)<-c("chr", "start", "end", "value1")
RnD3_92$chr <- paste0("chr", RnD3_92$chr)

RnD3_168<-read.table(file = "RnD3-168.call.cns", sep="\t", head=T)
RnD3_168<-RnD3_168[,-c(4,6,7,8,9,10)]
colnames(RnD3_168)<-c("chr", "start", "end", "value1")
RnD3_168$chr <- paste0("chr", RnD3_168$chr)

RnD3_203<-read.table(file = "RnD3-203.call.cns", sep="\t", head=T)
RnD3_203<-RnD3_203[,-c(4,6,7,8,9,10)]
colnames(RnD3_203)<-c("chr", "start", "end", "value1")
RnD3_203$chr <- paste0("chr", RnD3_203$chr)

##############################################

RnD3_212<-read.table(file = "RnD3-212.call.cns", sep="\t", head=T)
RnD3_212<-RnD3_212[,-c(4,6,7,8,9,10)]
colnames(RnD3_212)<-c("chr", "start", "end", "value1")
RnD3_212$chr <- paste0("chr", RnD3_212$chr)

RnD3_232<-read.table(file = "RnD3-232.call.cns", sep="\t", head=T)
RnD3_232<-RnD3_232[,-c(4,6,7,8,9,10)]
colnames(RnD3_232)<-c("chr", "start", "end", "value1")
RnD3_232$chr <- paste0("chr", RnD3_232$chr)

RnD3_233<-read.table(file = "RnD3-233.call.cns", sep="\t", head=T)
RnD3_233<-RnD3_233[,-c(4,6,7,8,9,10)]
colnames(RnD3_233)<-c("chr", "start", "end", "value1")
RnD3_233$chr <- paste0("chr", RnD3_233$chr)

##############################################
# Create a circos plot GROUP D3
layout(matrix(1:2, 1, 2))
for(i in 1:2) {
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  circos.par(cell.padding = c(0, 0, 0, 0))
  circos.clear()
  circos.initializeWithIdeogram(species = "rn6")
  circos.genomicTrack(RnD3_71,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD3_92,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD3_168,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD3_203,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  
  ##############################################
  circos.clear()
  circos.initializeWithIdeogram(species = "rn6")
  circos.genomicTrack(RnD3_212,  track.height = 0.2,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD3_232,  track.height = 0.2,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD3_233,  track.height = 0.2,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
}




##############################################
##############################################
##############################################
# Create a circos plot GROUP D6
setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-RATS/cnv/cns/D6")
RnD6_132<-read.table(file = "RnD6-132.call.cns", sep="\t", head=T)
RnD6_132<-RnD6_132[,-c(4,6,7,8,9,10)]
colnames(RnD6_132)<-c("chr", "start", "end", "value1")
RnD6_132$chr <- paste0("chr", RnD6_132$chr)

RnD6_152<-read.table(file = "RnD6-152.call.cns", sep="\t", head=T)
RnD6_152<-RnD6_152[,-c(4,6,7,8,9,10)]
colnames(RnD6_152)<-c("chr", "start", "end", "value1")
RnD6_152$chr <- paste0("chr", RnD6_152$chr)

RnD6_179<-read.table(file = "RnD6-179.call.cns", sep="\t", head=T)
RnD6_179<-RnD6_179[,-c(4,6,7,8,9,10)]
colnames(RnD6_179)<-c("chr", "start", "end", "value1")
RnD6_179$chr <- paste0("chr", RnD6_179$chr)

RnD6_204<-read.table(file = "RnD6-204.call.cns", sep="\t", head=T)
RnD6_204<-RnD6_204[,-c(4,6,7,8,9,10)]
colnames(RnD6_204)<-c("chr", "start", "end", "value1") 
RnD6_204$chr <- paste0("chr", RnD6_204$chr)


##############################################
# Create a circos plot GROUP D6
layout(matrix(1:1, 1, 1))
circos.clear()
circos.initializeWithIdeogram(species = "rn6")
circos.genomicTrack(RnD6_132,  track.height = 0.175,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "red", "green"), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })
circos.genomicTrack(RnD6_152,  track.height = 0.175,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "red", "green"), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })
circos.genomicTrack(RnD6_179,  track.height = 0.175,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "red", "green"), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })
circos.genomicTrack(RnD6_204,  track.height = 0.175,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "red", "green"), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })



##############################################
##############################################
##############################################
# Create a circos plot GROUP D12
setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-RATS/cnv/cns/D12")
RnD12_49<-read.table(file = "RnD12-49.call.cns", sep="\t", head=T)
RnD12_49<-RnD12_49[,-c(4,6,7,8,9,10)]
colnames(RnD12_49)<-c("chr", "start", "end", "value1")
RnD12_49$chr <- paste0("chr", RnD12_49$chr)

RnD12_51<-read.table(file = "RnD12-51.call.cns", sep="\t", head=T)
RnD12_51<-RnD12_51[,-c(4,6,7,8,9,10)]
colnames(RnD12_51)<-c("chr", "start", "end", "value1")
RnD12_51$chr <- paste0("chr", RnD12_51$chr)

RnD12_145<-read.table(file = "RnD12-145.call.cns", sep="\t", head=T)
RnD12_145<-RnD12_145[,-c(4,6,7,8,9,10)]
colnames(RnD12_145)<-c("chr", "start", "end", "value1")
RnD12_145$chr <- paste0("chr", RnD12_145$chr)

##############################################

RnD12_178<-read.table(file = "RnD12-178.call.cns", sep="\t", head=T)
RnD12_178<-RnD12_178[,-c(4,6,7,8,9,10)]
colnames(RnD12_178)<-c("chr", "start", "end", "value1")
RnD12_178$chr <- paste0("chr", RnD12_178$chr)

RnD12_184<-read.table(file = "RnD12-184.call.cns", sep="\t", head=T)
RnD12_184<-RnD12_184[,-c(4,6,7,8,9,10)]
colnames(RnD12_184)<-c("chr", "start", "end", "value1")
RnD12_184$chr <- paste0("chr", RnD12_184$chr)

RnD12_229<-read.table(file = "RnD12-229.call.cns", sep="\t", head=T)
RnD12_229<-RnD12_229[,-c(4,6,7,8,9,10)]
colnames(RnD12_229)<-c("chr", "start", "end", "value1")
RnD12_229$chr <- paste0("chr", RnD12_229$chr)

RnD12_231<-read.table(file = "RnD12-231.call.cns", sep="\t", head=T)
RnD12_231<-RnD12_231[,-c(4,6,7,8,9,10)]
colnames(RnD12_231)<-c("chr", "start", "end", "value1")
RnD12_231$chr <- paste0("chr", RnD12_231$chr)

##############################################

RnD12_185_1<-read.table(file = "RnD12-185-1.call.cns", sep="\t", head=T)
RnD12_185_1<-RnD12_185_1[,-c(4,6,7,8,9,10)]
colnames(RnD12_185_1)<-c("chr", "start", "end", "value1")
RnD12_185_1$chr <- paste0("chr", RnD12_185_1$chr)

RnD12_185_2<-read.table(file = "RnD12-185-2.call.cns", sep="\t", head=T)
RnD12_185_2<-RnD12_185_2[,-c(4,6,7,8,9,10)]
colnames(RnD12_185_2)<-c("chr", "start", "end", "value1")
RnD12_185_2$chr <- paste0("chr", RnD12_185_2$chr)

##############################################
# Create a circos plot GROUP D12

layout(matrix(1:3, 1, 3))
for(i in 1:3) {
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  circos.par(cell.padding = c(0, 0, 0, 0))
  circos.clear()
  circos.initializeWithIdeogram(species = "rn6")
  circos.genomicTrack(RnD12_49,  track.height = 0.2,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD12_51,  track.height = 0.2,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD12_145,  track.height = 0.2,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  
  ##############################################
  
  circos.clear()
  circos.initializeWithIdeogram(species = "rn6")
  circos.genomicTrack(RnD12_178,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD12_184,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD12_229,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD12_231,  track.height = 0.175,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  
  ##############################################
  
  circos.clear()
  circos.initializeWithIdeogram(species = "rn6")
  circos.genomicTrack(RnD12_185_1,  track.height = 0.3,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(RnD12_185_2,  track.height = 0.3,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
}






##############################################
##############################################
##############################################
# Create a circos plot GROUP Pr
setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-RATS/cnv/cns/Pr")
RnPr109<-read.table(file = "RnPr109.call.cns", sep="\t", head=T)
RnPr109<-RnPr109[,-c(4,6,7,8,9,10)]
colnames(RnPr109)<-c("chr", "start", "end", "value1")
RnPr109$chr <- paste0("chr", RnPr109$chr)

RnPr140<-read.table(file = "RnPr140.call.cns", sep="\t", head=T)
RnPr140<-RnPr140[,-c(4,6,7,8,9,10)]
colnames(RnPr140)<-c("chr", "start", "end", "value1")
RnPr140$chr <- paste0("chr", RnPr140$chr)

RnPr146<-read.table(file = "RnPr146.call.cns", sep="\t", head=T)
RnPr146<-RnPr146[,-c(4,6,7,8,9,10)]
colnames(RnPr146)<-c("chr", "start", "end", "value1")
RnPr146$chr <- paste0("chr", RnPr146$chr)

RnPr214<-read.table(file = "RnPr214.call.cns", sep="\t", head=T)
RnPr214<-RnPr214[,-c(4,6,7,8,9,10)]
colnames(RnPr214)<-c("chr", "start", "end", "value1")
RnPr214$chr <- paste0("chr", RnPr214$chr)

##############################################
# Create a circos plot GROUP Pr

layout(matrix(1:1, 1, 1))
circos.clear()
circos.initializeWithIdeogram(species = "rn6")
circos.genomicTrack(RnPr109,  track.height = 0.15,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "red", "green"), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })
circos.genomicTrack(RnPr140,  track.height = 0.15,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "red", "green"), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })
circos.genomicTrack(RnPr214,  track.height = 0.15,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "red", "green"), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })
circos.genomicTrack(RnPr214,  track.height = 0.15,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "red", "green"), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })


########################################################################################################################################################################################
# PERCENTAGE OF GAINS OR LOSSES IN EACH CHROMOSOME OVERALL

RnD1_104
Rn_names <- ls(pattern = "^Rn")

# filter to only keep data frames
Rn_df_names <- Rn_names[sapply(Rn_names, function(x) is.data.frame(get(x)))]

# retrieve each data frame and store in a list
Rn_dfs <- lapply(Rn_df_names, get)

# concatenate all data frames into a single data frame
combined_df <- do.call(rbind, Rn_dfs)

combined_df$chr <- sub("chr", "", combined_df$chr)
combined_df$chr <- as.character(combined_df$chr)
combined_df$chr <- factor(combined_df$chr, levels=unique(combined_df$chr))

library(magrittr)
library(ggforestplot)
library(ggplot2)
library(tidyverse)

ggplot(data = combined_df, aes(x = chr, y = value1, color = value1 > 0)) +
  theme_minimal()+
  theme(legend.title = element_blank(),panel.grid.major.y=element_blank(), panel.grid.major.x = element_line(color = "lightgray", linewidth=1))+
  scale_x_discrete(breaks = c(2,4,6,8,10,12,14,16,18,20,"Y"))+
  geom_jitter(alpha=0.5) +
  labs(x = "Chromosome", y = "log2")+
  scale_color_manual(values = c("#4ae153", "#e1534a"), labels = c("Loss", "Gain"))




# calculate the percentage of positive and negative value1 for each unique chr
result <- aggregate(
  value1 ~ chr,
  data = combined_df,
  FUN = function(x) {
    n <- length(x)
    pos_pct <- mean(x > 0) * 100
    neg_pct <- mean(x < 0) * 100
    c(pos_pct, neg_pct)
  }
)

# print the result
head(result)

res<-as.data.frame(result)

typeof(result)

setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-RATS/cnv/cns")
data<-read.csv(file="percentage_gain_loss.CSV",sep=",",header=T)

data$chromosome <- sub("chr", "", data$chromosome)
data$chromosome <- as.character(data$chromosome)
data$chromosome <- factor(data$chromosome, levels=unique(data$chromosome))

data<-melt(data,id.vars = "chromosome")

ggplot(data, aes(x = chromosome, y = value, fill = variable)) +
  geom_bar(stat="identity", position="stack") +
  labs(x = "Chromosome", y = "Percentage", title = "Gains and losses by chromosome") +
  scale_fill_manual(values = c("#e1534a", "#4ae153"))+
  theme_minimal()

########################################################################################################################################################################################

##############################################
##############################################
##############################################
# Create a circos plot GROUP D12
setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-RATS/cnv/cns/D12")

bed<-read.table(file = "RnD12.call.cns", sep="\t", head=T)
colnames(bed)<-c("chr", "start", "end", "value1")
bed$chr <- paste0("chr", bed$chr)
head(bed)

layout(matrix(1:1, 1, 1))
circos.clear()
circos.initializeWithIdeogram(species = "rn6")
circos.genomicTrack(bed,  track.height =0.5,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "red", "green"), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })

layout(matrix(1:1, 1, 1))
circos.clear()
circos.initializeWithIdeogram(species = "rn6")
circos.genomicTrack(bed, track.height =0.5,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, area = T)
                    })
