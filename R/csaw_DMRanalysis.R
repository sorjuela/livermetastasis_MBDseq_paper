#!/usr/bin/env Rscript

#########################################################################################
# R script to detect differential methylation using the csaw package.
# csaw is for differential binding, so the assumption is that this analysis 
# is DB for the methyl-binding protein. I think this makes sense.
#
# MBD-seq with 7 metastatis untreated, 5 metastasis treated, 6 CRC, 3 normal mucosa
# Collaboration with Mirco Menigatti and Giancarlo Marra
#
# Useful package explanation and guidelines here:
# https://master.bioconductor.org/help/course-materials/2015/BioC2015/csaw_lab.html
#
# Stephany Orjuela, August 2018
#########################################################################################

library("csaw")
library("edgeR")
library("annotatr")
library("matrixStats")
library("parallel")
library("ggplot2")


#Get metadata
samples <- read.table("../Data/metadata.txt", header =TRUE,sep = "\t")
samples$file <- paste0("../bam_dedup/dedups_local/", 
                       list.files("../bam_dedup/dedups_local/", "_s.bam$"))

#### Counting Number of Reads in Each Bin ####----------------------------------

#cross-correlation plot
max.delay <- 500
param <- readParam(restrict = paste0("chr",seq(1:22)))
dedup.on <- reform(param, dedup=TRUE)
x <- correlateReads(as.character(samples$file), max.delay, param=dedup.on)

pdf("myFigs/cross-correlation.pdf")
plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)")
dev.off()

maximizeCcf(x)
#180

#window count
data <- windowCounts(as.character(samples$file), ext = 180,
                     width = 200, spacing = 10, shift = 0,
                     param=param, filter=30) #743,44,976

colnames(data) <- samples$sample

#save(data, file = "../Data/csaw_data_nomapqfilt.RData")
#load("../Data/csaw_data_nomapqfilt.RData")

#### Filter ####----------------------------------------------------------------

#Independent filtering for count data
abundances <- aveLogCPM(asDGEList(data))
keep.simple <- abundances > -1
data2 <- data[keep.simple,]

rm(data)

#### Normalize ####-------------------------------------------------------------

#Deal with trended bias
data2 <- normOffsets(data2, type = "loess", se.out=T)
data.off <- assay(data2, "offset")

pdf("myFigs/pre_norm_verify.pdf")

# MA plot without normalization.
ac.y <- asDGEList(data2)
adjc <- cpm(ac.y, log=TRUE)
abval <- aveLogCPM(ac.y)
mval <- adjc[,1]-adjc[,2]
fit <- loessFit(x=abval, y=mval)
smoothScatter(abval, mval, ylab="M", xlab="Average logCPM",
              main="Raw", ylim=c(-2,2), xlim=c(0, 7))
o <- order(abval)
lines(abval[o], fit$fitted[o], col="red")

# Repeating after normalization.
re.adjc <- log2(assay(data2)+0.5) - data.off/log(2)
mval <- re.adjc[,1]-re.adjc[,2]
fit <- loessFit(x=abval, y=mval)
smoothScatter(abval, re.adjc[,1]-re.adjc[,2], ylab="M", xlab="Average logCPM",
              main="Normalized", ylim=c(-2,2), xlim=c(0, 7))
lines(abval[o], fit$fitted[o], col="red")
dev.off()
 
# save(data2, file = "../Data/csaw_data_nomapqfilt_abundfilt.RData")

#### Testing for differential binding ####--------------------------------------

#load("../Data/csaw_data_nomapqfilt_abundfilt.RData")

#Turn to DGElist
y <- asDGEList(data2)

#Make model
grp <- as.character(samples$condition)
grp2 <- factor(ifelse(grp == "metastasis_untreated" | grp == "metastasis_treated", 
                      "metastasis", grp))
grp2
 
mm <- model.matrix(~ 0 + grp2)

#Estimate dispersions
y <- estimateDisp(y, mm) 
#load("csaw_DGEwithDisp_nomapqfilt.RData")

fit <- glmQLFit(y, mm, robust=TRUE)
save(fit, file = "csaw_fit.RData")


png("myFigs/AveLogCPM.png")
o <- order(y$AveLogCPM)
plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
     ylim=c(0, 1), xlab=expression("Ave."~Log[2]~"CPM"),
     ylab=("Biological coefficient of variation"))
plotQLDisp(fit)
dev.off()

#Testing for each contrast
#treated and untreated together

mc <- makeContrasts(cn = "colorectal_cancer-normal_mucosa",
                    mn = "metastasis-normal_mucosa",
                    mc = "metastasis-colorectal_cancer",
                    levels=levels(grp2))


lrts1 <- mclapply(as.data.frame(mc), function(u) {
  lrt <- glmQLFTest(fit, contrast=u)
}, mc.cores = 3)

#### Correct multiple testing ####-----------------------------------------------

#merge windows and combine p-vals when logical
merged <- mergeWindows(rowRanges(data2), tol=50L)

tables <- lapply(lrts1,function(x){ return(x$table)})

mrg_id <- list(merged$id)[rep(1,3)]
tabcom <- Map( combineTests, mrg_id, tables)

#calculate mean logFC per merged window for cancer Vs normal
#change to predFC
grp2 <- relevel(grp2, "normal_mucosa")
mm <- model.matrix(~grp2)
predlfc <- predFC(y, mm, prior.count = 1)

# cntable <- lrts1$cn$table$logFC
cntable <- predlfc[,2]
uniqueid <- unique(sort(merged$id))
meanlogFC <- lapply(uniqueid, function(x){
  w <- which(merged$id == x)
  out <- mean(cntable[w])
  return(out)
})

#### MDS ####------------------------------------------------------------------

#Calculate mean cpm per merged window to calculate MDS

cpms <- cpm(y, log=TRUE)

uniqueid <- unique(sort(merged$id))
mean_cpms <- sapply(uniqueid, function(i){
  w <- merged$id == i
  if(sum(w) > 1) {
      matrixStats::colMeans2(cpms[merged$id == i,])
  } else{
    cpms[merged$id == i,]
  }
})

tmean_cpms <- t(mean_cpms)

#get MDS
mds <- plotMDS(tmean_cpms, top = 10000, plot = FALSE)$cmdscale.out
df <- data.frame(dim1 = mds[,1], dim2 = mds[,2],
                 names = samples$sample, Tissue = samples$condition)

ggplot()+
  geom_point(data = df, mapping = aes_(x=~dim1, y=~dim2, color=~Tissue),
             size = 5) +
  geom_text(data = df, mapping = aes_(x=~dim1, y=~dim2-0.1, label=~names),
            size = 3) +
  theme_bw() +
  
ggsave("myFigs/MDS_2D_samp.png")

#### Export tables ####---------------------------------------------------------
#Get locations
res <- as.data.frame(merged$region)

#export bed file of regions
write.table(res[,1:3], "MBD_regions.bed", row.names=FALSE, quote=FALSE, sep="\t")

# Select annotations for intersection with regions
#CpG annot
annotscpg = c("hg19_cpg_islands", "hg19_cpg_shores",
              "hg19_cpg_shelves", "hg19_cpg_inter")

annotations = build_annotations(genome = 'hg19', annotations = annotscpg)

dm_annotated = annotate_regions(
  regions = merged$region,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

df_dm_annotated <- data.frame(dm_annotated)

over <- findOverlaps(merged$region, dm_annotated)
uniQ <- unique(queryHits(over))
Q <- queryHits(over)
S <- subjectHits(over)

cpgAnnot <- sapply(uniQ, function(w){
  sub <- df_dm_annotated[S[Q == w],]
  paste(unique(sub$annot.id), collapse="/")
})

res$CpG <- "NA"
res$CpG[uniQ] <- cpgAnnot
# 
# #genes
annotsgene <- c("hg19_genes_promoters", "hg19_genes_3UTRs", "hg19_genes_introns",
                "hg19_genes_exons", "hg19_genes_5UTRs", "hg19_genes_cds", "hg19_genes_intergenic",
                "hg19_genes_1to5kb")

annotations = build_annotations(genome = 'hg19', annotations = annotsgene)

dm_annotated = annotate_regions(
  regions = merged$region,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

df_dm_annotated <- data.frame(dm_annotated)

over <- findOverlaps(merged$region, dm_annotated)
uniQ <- unique(queryHits(over))
Q <- queryHits(over)
S <- subjectHits(over)

geneAnnot <- sapply(uniQ, function(w){
  sub <- df_dm_annotated[S[Q == w],]
  splited <- limma::strsplit2(sub$annot.id, ":")
  j <- paste(unique(paste(splited[,1], sub$annot.symbol, sep = ".")), collapse = "/")
  return(j)
})

res$gene <- "NA"
res$gene[uniQ] <- geneAnnot

#add cn.meanlogFC
res$cn.meanlogFC <- unlist(meanlogFC)

#Get stats
for(i in 1:length(lrts1)) {
  df <- cbind(tabcom[[i]], de = ifelse(tabcom[[i]]$FDR <= 0.05, 1, 0))
  colnames(df) <- paste(names(lrts1)[i], colnames(df), sep=".")
  res <- cbind(res, df)
}

#save(res, file ="unsorted.res.olddesign_predFC.RData") 

#sort
pvals <- res[, grep("PValue", colnames(res))]
rm <- rowMins(as.matrix(pvals))
o <- order(rm)
res <- res[o,] #322551


#save txt and bed files
write.table(res[res$cn.de == 1,1:3], file = "../Data/MBD_csaw_CRCVsNORM_DMRs.bed",
            col.names = FALSE,
            row.names=FALSE,
            quote=FALSE, sep="\t")
save(res, file = "../Data/MBD_csaw_allregions.RData")
write.table(res, "../Data/MBD_csaw_allregions.csv", row.names=FALSE, 
            quote=FALSE, sep="\t")

