#!/usr/bin/env Rscript

#########################################################################################
# R script to compare metastasis data from other sources to ours
#
# MBD-seq with 7 metastatis untreated, 5 metastasis treated, 6 CRC, 3 normal mucosa
# Collaboration with Mirco Menigatti and Giancarlo Marra
#
#
# Stephany Orjuela, January 2020
#########################################################################################

library(dplyr)
#library(ComplexHeatmap)
library(limma)
library(ggplot2)

#### Timp data ####

f <- list.files("../Data/Validation_data/timp/", ".txt")
samp <- gsub("GSM[0-9]+\\-[0-9]+_colon_[a-z]+_","",f, perl = TRUE) %>%
        gsub(pattern = "\\.txt", replacement = "", perl = TRUE)
tissue <- gsub("GSM[0-9]+\\-[0-9]+_colon_","",f, perl = TRUE) %>%
  gsub(pattern = "_[0-9]{1}\\.txt", replacement = "", perl = TRUE)

metadata <- data.frame(files = paste0("../Data/Validation_data/timp/",f),
                       patient = samp,
                       state = tissue,
                       stringsAsFactors = FALSE,
                       sample = paste0(samp, "_",tissue))

#read all files into list
tablist <- lapply(1:18, function(f){
  read.table(metadata$files[f], header = TRUE)
}
)

all_keys <- lapply(tablist, function(f){
  f$ID_REF
})

key <- as.character(unique(unlist(all_keys)))

#join all samples in a table
betas <- mapply( function(df,k){
  m <- match(key, k)
  df$VALUE[m]
}, tablist, all_keys)

colnames(betas) <- metadata$sample
rownames(betas) <- key

save(betas, file = "../Data/Validation_data/timp/betas_matrix.RData")

load("../Data/Validation_data/timp/betas_matrix.RData")

#betas <- betas[,order(colnames(betas))]
dat <- reshape2::melt(betas)
dat$Tissue <- metadata$state[match(dat$Var2,metadata$sample)]
#ggplot(dat) + geom_density(aes(value, color = Tissue, group = Var2)) + theme_bw()

ggplot(dat) + 
  geom_violin(aes(Var2,value, fill= Tissue, color = Tissue)) +
  theme_bw()

#MDS
#mval <- log2(betas/(1-betas))

x <- limma::plotMDS(betas, top = 1000, plot = FALSE, dim.plot = c(1,2))
mdsscale <- data.frame(Dim1 = x$cmdscale.out[,1], 
                       Dim2 = x$cmdscale.out[,2],
                       Samples = rownames(x$cmdscale.out),
                       Tissue = metadata$state)
timpmds <- ggplot(mdsscale) +
  geom_point(aes(Dim1,Dim2, color=Tissue), size = 6) +
  geom_text(aes(Dim1,Dim2-0.01, label =Samples), size = 3.5) + 
  labs(x="dimension 1", y = "dimension 2") +
  theme_bw()
#ggplot2::ggsave("myFigs/Timpdat_MDS.png")
timpmds

#diagnostic plots
means <- rowMeans(betas)
vars <- matrixStats::rowVars(betas)
plot(vars,means)

#Make plot to show similarity between metastasis and cancer

x <- data.frame(
  Met = rowMeans(betas[,metadata$state == "metastasis"]),
  CRC = rowMeans(betas[,metadata$state == "cancer"]),
  NM = rowMeans(betas[,metadata$state == "normal"]),
  loc = "Genome"
)

#Annotate probes
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Manifest)
data(Islands.UCSC)
#data(Locations)
data(Other)

y <- match(rownames(betas),Manifest$Name)

anno <- data.frame(probe = Manifest$Name[y], 
                   island = Islands.UCSC$Relation_to_Island[y],
                   gene = Other$UCSC_RefGene_Name[y],
                   annot = Other$UCSC_RefGene_Group[y]
)
betas_isles <- betas[anno$island == "Island",]

xisles <- data.frame(
  Met = rowMeans(betas_isles[,metadata$state == "metastasis"]),
  CRC = rowMeans(betas_isles[,metadata$state == "cancer"]),
  NM = rowMeans(betas_isles[,metadata$state == "normal"]),
  loc = "Islands"
)

x <- rbind(x,xisles)

a <- ggplot(x, aes(NM, Met)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') +
  geom_abline() +
  facet_wrap(~loc, nrow = 2) + 
  theme_bw()+
  labs(x="mean beta value NM", y="mean beta value Met")+
  theme(legend.position = "none",panel.spacing = unit(0, "lines"))

b <- ggplot(x, aes(NM, CRC)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') +
  geom_abline() +
  facet_wrap(~loc, nrow = 2) +
  theme_bw() +
  labs(x="mean beta value NM", y="mean beta value CRC")+
  theme(legend.position = "none", panel.spacing = unit(0, "lines"))

c <- ggplot(x, aes(Met, CRC)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') +
  geom_abline() +
  facet_wrap(~loc, nrow = 2) +
  theme_bw() +
  labs(x="mean beta value Met", y="mean beta value CRC")+
  theme(panel.spacing = unit(0, "lines"))
  

mid <- cowplot::plot_grid(a,b,c, ncol = 3, rel_widths = c(0.75,0.75,1))
#ggplot2::ggsave("myFigs/Timpdat_betaScatters.png", width = 10, height = 5)


#Make Figure4C-D from Timp paper
x <- data.frame(meth = c(colMeans(betas),colMeans(betas_isles)),
                loc = c(rep("Genome",18),rep("Islands",18)),
                Samples = rep(gsub("_[a-z]+","",colnames(betas)),2),
                Tissue= factor(rep(metadata$state,2), levels=c("normal",
                                                        "cancer",
                                                        "metastasis")))

mean_mean <- aggregate(meth ~ loc+Tissue, mean, data=x)
cols <- scales::hue_pal()(3)[c(3,1,2)]

jitt <- ggplot(x,aes(Tissue,meth)) +
  geom_jitter(aes(color = Tissue), size = 3, width = 0.1) +
  #geom_line(aes(group = Samples), color = "gray") +
  #geom_point(aes(color = Tissue), size = 3) +
  theme_bw() +
  geom_crossbar(data=mean_mean, aes(ymin = meth, ymax = meth, color = Tissue),
                width = 0.5) +
  geom_text(data=mean_mean, aes(label = round(meth, digits = 3)),
            vjust = 1, hjust = 0) +
  facet_wrap(~loc,nrow = 2, scales = "free") +
  ylab("mean beta value per sample") +
  scale_color_manual(values = cols)

jitt

cowplot::plot_grid(timpmds,mid,jitt, ncol = 1, nrow = 3, labels = "AUTO")
ggplot2::ggsave("myFigs/Timpdat_suppfig.png", width = 8, height = 11)


#### Martinez-C data ####

#Use MD (TS) samples
f <- list.files("../Data/Validation_data/martinez/", ".txt")
samp <- gsub("GSM[0-9]+\\-[0-9]+_col","",f, perl = TRUE) %>%
  gsub(pattern = "[A-Z]+\\.txt", replacement = "", perl = TRUE)
tissue <- gsub("GSM[0-9]+\\-[0-9]+_col[0-9]{2}","",f, perl = TRUE) %>%
  gsub(pattern = "\\.txt", replacement = "", perl = TRUE)

metadata <- data.frame(files = paste0("../Data/Validation_data/martinez/",f),
                       patient = samp,
                       state = tissue,
                       stringsAsFactors = FALSE,
                       sample = paste0(samp, "_",tissue))

#read all files into list
tablist <- lapply(1:12, function(f){
  read.table(metadata$files[f], header = TRUE, sep = "\t")
}
)

all_keys <- lapply(tablist, function(f){
  f$ID_REF
})

key <- as.character(unique(unlist(all_keys)))

#join all samples in a table
betas <- mapply( function(df,k){
  m <- match(key, k)
  df$VALUE[m]
}, tablist, all_keys)

colnames(betas) <- metadata$sample
rownames(betas) <- key

save(betas, file = "../Data/Validation_data/martinez/betas_matrix.RData")

#Transform beta vals?
#betas_trans <- asin(2*betas-1)
#bad <- BiocGenerics::rowSums(is.finite(betas_trans)) < ncol(betas_trans)
#if(any(bad)) betas_trans <- betas_trans[!bad,,drop=FALSE]

# Get mvals?
#mval <- log2(betas/(1-betas))

# MDS
x <- limma::plotMDS(betas, top = 1000, plot = FALSE, dim.plot = c(1,2))
mdsscale <- data.frame(Dim1 = x$cmdscale.out[,1], 
                       Dim2 = x$cmdscale.out[,2],
                       Samples = rownames(x$cmdscale.out),
                       Tissue = factor(metadata$state,
                                       levels = (c("MD","S","N"))))
mds <- ggplot(mdsscale) +
  geom_point(aes(Dim1,Dim2, color=Tissue), size = 6) +
  geom_text(aes(Dim1,Dim2-0.02, label =Samples)) + 
  scale_color_discrete(name = "Tissue") +
  theme_bw()
mds
#ggplot2::ggsave("myFigs/Martinezdat_MDS.png")

#Plots
x <- data.frame(
  Met = rowMeans(betas[,metadata$state == "S"]),
  CRC = rowMeans(betas[,metadata$state == "MD"]),
  NM = rowMeans(betas[,metadata$state == "N"]),
  loc = "Genome"
)

y <- match(rownames(betas),Manifest$Name)

anno <- data.frame(probe = Manifest$Name[y], 
                   island = Islands.UCSC$Relation_to_Island[y],
                   gene = Other$UCSC_RefGene_Name[y],
                   annot = Other$UCSC_RefGene_Group[y]
)
betas_isles <- betas[anno$island == "Island",]

xisles <- data.frame(
  Met = rowMeans(betas_isles[,metadata$state == "S"]),
  CRC = rowMeans(betas_isles[,metadata$state == "MD"]),
  NM = rowMeans(betas_isles[,metadata$state == "N"]),
  loc = "Islands"
)

x <- rbind(x,xisles)

a <- ggplot(x, aes(NM, Met)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') +
  geom_abline() +
  facet_wrap(~loc, nrow = 2) + 
  theme_bw()+
  theme(legend.position = "bottom",
        legend.box = "vertical")

b <- ggplot(x, aes(NM, CRC)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') +
  geom_abline() +
  facet_wrap(~loc, nrow = 2) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical")

c <- ggplot(x, aes(Met, CRC)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') +
  geom_abline() +
  facet_wrap(~loc, nrow = 2) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical")


mid <- cowplot::plot_grid(a,b,c, ncol = 3)
mid

#Make Figure4C-D from Timp paper
x <- data.frame(meth = c(colMeans(betas, na.rm = TRUE),
                         colMeans(betas_isles, na.rm = TRUE)),
                loc = c(rep("Genome",12),rep("Islands",12)),
                Samples = rep(colnames(betas),2),
                Tissue= factor(rep(metadata$state,2), levels=c("N",
                                                               "MD",
                                                               "S")))

mean_mean <- aggregate(meth ~ loc+Tissue, mean, data=x)
cols <- scales::hue_pal()(3)[c(3,1,2)]

jitt <- ggplot(x,aes(Tissue,meth,color = Tissue)) +
  geom_jitter(width = 0.15) +
  theme_bw() +
  geom_crossbar(data=mean_mean, aes(ymin = meth, ymax = meth, color = Tissue),
                width = 0.5) +
  facet_wrap(~loc,nrow = 2, scales = "free") +
  ylab("Average methylation") +
  scale_color_manual(values = cols)

jitt

cowplot::plot_grid(mds,mid,jitt, ncol = 1, nrow = 3, labels = "AUTO",
                   axis = "r", align ="v", rel_widths = c(1,2,1))
ggplot2::ggsave("myFigs/Martinezdat_suppfig.png", width = 8, height = 11)


