

#### Upset plot within islands ####

islands <- grep("island|shore|shelf", res$CpG)
res_islands <- res$CpG[islands]

test <- res_islands[1:10]
isleMat <- matrix(0, 28691, 3)
rownames(isleMat) <- c(1:28691)
colnames(isleMat) <- c("Island","Shore","Shelf")


sapply(test, function(c){
  isHit <- stringr::str_extract_all(test[1], "\\/*island:[0-9]+\\/*")[[1]]
  
  if(isHit){
    isNum <- as.numeric(stringr::str_extract_all(isHit, "[0-9]+")[[1]])
  }
  
  shoHit <- stringr::str_extract_all(test[1], "\\/*shore:[0-9]+\\/*")[[1]]
  shoNum <- as.numeric(stringr::str_extract_all(shoHit, "[0-9]+")[[1]])
  isID <- shoNum - floor(shoNum/2)
  
  sheHit <- stringr::str_extract_all(test[1], "\\/*shelf:[0-9]+\\/*")[[1]]
  sheNum <- as.numeric(stringr::str_extract_all(sheHit, "[0-9]+")[[1]])
  isID <- sheNum - floor(sheNum/2)
})

x <- 15634
isleID <- x - floor(x/2)


#### Count number of single island ids
library(annotatr)
annotscpg = c("hg19_cpg_islands", "hg19_cpg_shores",
              "hg19_cpg_shelves", "hg19_cpg_inter")
annotations = build_annotations(genome = 'hg19', annotations = annotscpg)

shores <- annotations$id[grep("shore", annotations$id)]
DMR <- GRanges("chr6", IRanges(133562000,133563686))
subsetByOverlaps(annotations, DMR)
#conclusion: island 2 has shores 3 and 4 --->  3|3|3|2|4|4|4
#other example:   1|1|1|1|2|2|2

#x <- limma::strsplit2(annotations$id, ":")
#x <- paste0(":",unique(sort(x)),"/") #51918?, 28691 CpG islands...where do the rest come from??
x <- paste0(":",1:28691,"/") #just the island IDs

pres.islands <- sapply(x, function(c){
  hit <- grep(c, res_islands)
  return(hit[1])
})

sum(!is.na(pres.islands)) #21738 #32226??

#### Upset plot within genic regions in open sea ####

islands <- grep("^inter:[0-9]+$", res$CpG)
res_islands <- res[islands,] #246860
res_isleGR <- GRanges(res_islands$seqnames, IRanges(res_islands$start,
                                                    res_islands$end),
                      mcols = res_islands[,6:7])

UTR3 <- grep("3UTR", res_isleGR$mcols.gene)
res_3UTR <- res_isleGR[UTR3]

res_no3UTR <- res_isleGR[-UTR3]
proms <- grep("(promoter|1to5kb|5UTR)", res_no3UTR$mcols.gene, perl =T)
res_proms <- res_no3UTR[proms]

res_no3UTR_noproms <- res_no3UTR[-proms]
gene_body <- grep("(exon|intron|CDS)", res_no3UTR_noproms$mcols.gene, perl =T)
res_genebody <- res_no3UTR_noproms[gene_body]

#### islands without regions ####
annoshores <- annotations[grep("shore", annotations$id)]

resGR <- GRanges(res$seqnames, IRanges(res$start,
                                       res$end),
                 mcols = res[,6:7])
resGR <- sortSeqlevels(resGR)
resGR <- sort(resGR)
s <- as.data.frame(resGR)
write.table(s[,c(1:3,6:7)], "csaw_regions.bed", col.names = F, row.names=F, quote=FALSE, sep="\t")

islesNoRegs <- subsetByOverlaps(annoisles,resGR, invert = T) #Why the FFF???
s <- as.data.frame(islesNoRegs)
write.table(s[,1:6], "Islands_without_MBDregions.csv", row.names=FALSE, quote=FALSE, sep="\t")

#### Number of islands (either the island, or its shore, or shelf) "covered"
annoisles <- annotations[grep("island", annotations$id)] #28691
annoisles <- annoisles[seqnames(annoisles) %in% paste0("chr", c(1:23, "X", "Y"))] #27718
sloppedisles <- GenomicRanges::promoters(annoisles, upstream=4000, downstream=4000)

resGR <- GRanges(res$seqnames, IRanges(res$start,
                                       res$end),
                 mcols = res[,6])

over <- findOverlaps(sloppedisles, resGR)
length(sort(unique(queryHits(over)))) #24644

#### number of ssIss covered hyperDMRs by Roche or MBD ####
library(annotatr)
annotscpg = c("hg19_cpg_islands", "hg19_cpg_shores",
              "hg19_cpg_shelves", "hg19_cpg_inter")
annotations = build_annotations(genome = 'hg19', annotations = annotscpg)
annoisles <- annotations[grep("island", annotations$id)] #28691
annoisles <- annoisles[seqnames(annoisles) %in% paste0("chr", c(1:23, "X", "Y"))] #27718
sloppedisles <- GenomicRanges::promoters(annoisles, upstream=4000, downstream=4000)

load("MBD_csaw_verify_mod_3grps2_nomapqfilt.RData")
resGR <- csawGR <-  GRanges(res$seqnames, IRanges(res$start, res$end), direction = res$cn.direction)
resGR <- resGR[res$cn.de == 1] #2155
resGRup <- resGR[resGR$direction == "up"]

overmbd <- subsetByOverlaps(sloppedisles, resGRup)
length(sort(unique(queryHits(over)))) #2488

load("../../../../../sorjuela/serrated_pathway_paper/BiSulf/Data/nonCIMP.allclusters.trimmed.RData")
DMRs <- BiSeq::findDMRs(allClusters.trimmed, max.dist = 200, diff.dir = F) #19,478
probes <- read.table(file="../../../../../Shared_taupo/steph/reference/130912_HG19_CpGiant_4M_EPI.bed") #Available online from Roche
probes <- GRanges(probes$V1, IRanges(start = probes$V2, end = probes$V3)) #240.131
probes <- probes[!seqnames(probes) %in% c("chrX", "chrY", "chrM")]
DMRs <- subsetByOverlaps(DMRs,probes) #10,608
DMRsfilt <- DMRs[width(DMRs) >= 3]
DMRsfilthyper <- DMRsfilt[DMRsfilt$median.meth.diff > 0]

overprobes <- subsetByOverlaps(sloppedisles, DMRsfilthyper)
length(sort(unique(queryHits(over)))) #1069

over <- findOverlaps(overmbd, overprobes) #1521, qH:733, sH:658
