#!/usr/bin/env Rscript
library(readr)
library(data.table)
library(ppcor)

options(readr.num_columns = 0)

if (interactive()) {
  optMod="K04/K04_F1NA1"
  optConn="P" # Correlation, Partial Correlation, Semi-partial Correlation
  optROI="rann_LAI,rann_LCH6p,rann_LFEF,rann_LIFJ,rann_LIP3,rann_LMFC,rann_MCV6r,rann_RAI,rann_RCH6,rann_RFEF,rann_RIFJ,rann_RIPS,rann_RMFC"
  subjs <- "14533,14555,14556,14581,14582,14583,14707,14708,14710,14727,14728,14730,14765,14783,14828,14837,14877,14890,14891,14892,14904,14905,14906,14907,14909,14911,14913,14914,15669,19202,19203,19204,19205,19207,19209"
} else {
  args<-commandArgs(TRUE)
  optROI <- args[1]
  optConn <- args[2]
  optMod <- args[3]
  subjs <- args[4]
}

subjs<-strsplit(subjs, ",")[[1]]
optROI<-sub("^","T_",strsplit(optROI, ",")[[1]])

cat("<<<<<")
cat(paste0(optROI,collapse = " "))
cat(">>>>>\n")

#ts <- read_table2(paste0(optMod,"/",subjs[1],"A/ROI_R1.1D")) # get the columns from the first subject
cat(paste0(optMod,"/",subjs[1],"A/ROI_R1.1D\n"))

ts <- read_table2(paste0(optMod,"/",subjs[1],"A/ROI_R1.1D"),n_max = 1) # get the columns from the first subject
#colnames(ts)
ts<-ts[,optROI]
roiNames <- colnames(ts)
pathNames=sprintf("%s_%s",rep(roiNames,length(roiNames)), rep(roiNames,each=length(roiNames)))
#pathNames<-matrix(pathNames,length(roiNames))
#pathNames<-pathNames[upper.tri(pathNames)]

conn<-matrix(nrow = length(subjs)*4, ncol = length(pathNames)) # 2 modalities x 2 conditions


if (optConn=="P") {
  connest<-function(ts) {
    pcor(ts)$estimate
  }
} else if (optConn=="S") {
  connest<-function(ts) {
    spcor(ts)$estimate
  }
} else {
  connest<-cor
}

n=0
for (s in subjs){
  for (m in c("A","V")) {
    for(cond in c("/ROI_R1.1D","/ROI_R2.1D")) {
      n=n+1
      ts <- read_table2(paste0(optMod,"/",s,m,cond))
      ts<-ts[,optROI]
      conn[n,] <- connest(ts)
    }
  }
}

#conn
stat.ttest <- matrixTests::col_t_onesample(conn, alternative = "two.sided", conf.level = 0.95)
grpP <- matrix(data=stat.ttest$pvalue,nrow = length(roiNames), ncol = length(roiNames))
connMat <- matrix(data=stat.ttest$statistic,nrow = length(roiNames), ncol = length(roiNames))
rownames(grpP)<-colnames(grpP)<-roiNames
if (optConn=="C") grpP[lower.tri(grpP)]<-NA
zMat <- conn
if (!dir.exists(optMod)) dir.create(optMod)
cat(optMod)
save(connMat,grpP,zMat,roiNames,file = paste0(optMod,"/zMat.RData"))
#load("L01/L01_A1AA1/zMat.RData")

