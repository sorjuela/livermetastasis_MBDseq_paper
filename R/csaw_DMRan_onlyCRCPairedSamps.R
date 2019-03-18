#!/usr/bin/env Rscript
# chmod +x 
# run as [R < scriptName.R --no-save]

#########################################################################################
# R script to repeat differential methylation using the csaw package.
# removing the non-paired samples
#
# MBD-seq with 7 metastatis untreated, 5 metastasis treated, 6 CRC, 3 normal mucosa
# Collaboration with Mirco Menigatti and Giancarlo Marra
#
#
# Stephany Orjuela, March 2019
#########################################################################################

### screen -r cpgdensity 
library("parallel")
library("csaw")
library("edgeR")

samples <- read.table("../metadata.txt", header =T)
samples <- samples[order(samples$sample),] 
samples$file <- paste0("../bam_dedup/dedups_local/", list.files("../bam_dedup/dedups_local/", "_s.bam$"))

samples <- samples[c(1,2,17:20),]

#### Counting Number of Reads in Each Bin ####-----------------------------------------------

# param <- readParam(restrict = paste0("chr",seq(1:22)))
# # maximizeCcf(x)
# # #180
# 
# #window count
# data <- windowCounts(as.character(samples$file), ext = 180,
#                      width = 200, spacing = 10, shift = 0,
#                      param=param, filter=30)
# 
# colnames(data) <- samples$sample
# 
# save(data, file = "data_nomapqfilt_paired.RData")
load("data_nomapqfilt_paired.RData")

#### Filter ####--------------------------------------------------------------------

#Independent filtering for count data
abundances <- aveLogCPM(asDGEList(data))
keep.simple <- abundances > -1
data2 <- data[keep.simple,]

rm(data)

#merging region before diff to see difference after

#merged <- mergeWindows(rowRanges(data2), tol=50L) #679.854, after 322551

#### Normalize ####--------------------------------------------------------------

#Deal with trended bias
data2 <- normOffsets(data2, type = "loess", se.out=T)

data.off <- assay(data2, "offset")

png("myFigs/pre_norm_verify_paired.png")
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

save(data2, file = "csaw_data_nomapqfilt_abundfilt_paied.RData")

#remove non paired samples
#y <- asDGEList(data2[,-c(3:16,21)])
y <- asDGEList(data2)

#Make model
#grp <- as.character(samples$condition[-c(3:16,21)])
grp <- as.character(samples$condition)
grp <- as.factor(grp)
#grp2 <- factor(ifelse(grp == "metastasis_untreated" | grp == "metastasis_treated", "metastasis", grp))
#samps <- factor(c(1,1,rep(0,12), 2,2,3,3))
samps <- factor(c(1,1,2,2,3,3))

mm <- model.matrix(~ 0 + grp + samps)

#Estimate dispersions

y <- estimateDisp(y, mm) 
save(y, file = "csaw_DGEwithDisp_nomapqfilt_paired2.RData")

fit <- glmQLFit(y, mm, robust=TRUE)

png("myFigs/AveLogCPM_paired2.png")
o <- order(y$AveLogCPM)
plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
     ylim=c(0, 1), xlab=expression("Ave."~Log[2]~"CPM"),
     ylab=("Biological coefficient of variation"))
plotQLDisp(fit)
dev.off()

#Testing for each contrast
#treated and untreated together

mc <- makeContrasts(cn = "colorectal_cancer-normal_mucosa",
                    levels=levels(grp))


lrts1  <- glmQLFTest(fit)

merged <- mergeWindows(rowRanges(data2), tol=50L)

tables <- combineTests(merged$id, lrts1$table)

#calculate mean logFC per merged window for cancer Vs normal
cntable <- lrts1$table$logFC
uniqueid <- unique(sort(merged$id))
meanlogFC <- lapply(uniqueid, function(x){
  w <- which(merged$id == x)
  out <- mean(cntable[w])
  return(out)
})


res <- as.data.frame(merged$region)
df <- cbind(tables, de = ifelse(tables$FDR <= 0.05, 1, 0))
res <- cbind(res, df)

#add cn.meanlogFC

res$meanlogFC <- unlist(meanlogFC)

save(res, file ="unsorted.res.3Vs3crc.from0.RData")
#sort
# pvals <- res[, grep("PValue", colnames(res))]
# rm <- rowMins(as.matrix(pvals))
# o <- order(rm)
# res <- res[o,]


