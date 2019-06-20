#!/usr/bin/env Rscript

#########################################################################################
# R script to get DMCs using BiSeq (tBS-seq) 
#
# Collaboration with Mirco Menigatti and Giancarlo Marra
#
# Stephany Orjuela, January 2019
#########################################################################################


library("BiSeq")

#From previously generated BSraw-objects (Parker, et al., 2018)
#NOTE: I did all these loops because it was taking too long an I was getting annoyed

for(i in 1:22){ 
  
  load(file=paste0("../../../../../sorjuela/serrated_pathway_paper/BiSulf/Data/CRCdata/", 
                   paste0("chr",i,"BSraw.CRC.RData")))

  #add patient effect
  colData(BSrchr)$patient <- rep(1:6, each = 2)
  
  ### Define clusters ###----------------------------------------------------------------------------- 
  
  BSr.clust.unlim <- clusterSites(BSrchr,
                                  groups = factor(colData(BSrchr)$group),
                                  perc.samples = 0.8,
                                  min.sites = 3,
                                  max.dist =  500,
                                  mc.cores = 10)
  
  ### smooth the methylation values of CpG sites within the clusters ###----------------------------- 
  predictedMeth <- predictMeth(object = BSr.clust.unlim, h=1000, mc.cores = 10)
  #save(predictedMeth, file=paste0("betaRegs_Te_nocimpCRCsVsNorm/predictedMethchr",i,".RData"))
  
  ### Model methylation within a beta regression for each comparison ###-----------------------------
  
  #subset for groups of interest
  s <- grep("^non$", colData(predictedMeth)$group)
  n <- grep("NORMAL.non", colData(predictedMeth)$group)
  
  SNpMeth <- predictedMeth[,c(s,n)]
  colData(SNpMeth)$group <- factor(colData(SNpMeth)$group, levels = c("non", "NORMAL.non"))
  #save(SNpMeth, file=paste0("betaRegs_Te_nocimpCRCsVsNorm/predictedMethchr",i,".RData"))
  
  betaResultsSN <- betaRegression(formula = ~group+patient, link = "probit", object = SNpMeth, type = "BR", mc.cores = 10) 
  save(betaResultsSN, file=paste0("betaRegs_Te_nocimpCRCsVsNorm/betaRes.",i,".RData"))
}

#Then rbind all betaResults into single table

load("betaRegs_Te_nocimpCRCsVsNorm/betaRes.1.RData")

betaregAll <- betaResultsSN
for(i in 2:22){ 
  load(paste0("betaRegs_Te_nocimpCRCsVsNorm/betaRes.",i,".RData"))
  betaregAll <- rbind(betaregAll,betaResultsSN)
}
