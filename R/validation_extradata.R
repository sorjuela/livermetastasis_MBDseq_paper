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
library(ComplexHeatmap)
library(limma)

#### Timp data ####

f <- list.files("../Data/Validation_data/timp/", "txt")
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

dat <- reshape2::melt(mval)
dat$Tissue <- metadata$state[match(dat$Var2,metadata$sample)]
ggplot(dat) + geom_density(aes(value, color = Tissue, group = Var2)) + theme_bw()

#MDS
mval <- log2(betas/(1-betas))

limma::plotMDS(mval, top = 1000, col=as.integer(factor(metadata$state)))

#diagnostic plots
means <- rowMeans(betas)
vars <- matrixStats::rowVars(betas)
plot(vars,means)

#Do limma with betas
design <- model.matrix(~0+state + patient, data = metadata)
design
cont <- limma::makeContrasts(statemetastasis-statenormal, 
                             statecancer-statenormal, 
                             statemetastasis-statecancer,
                             levels = design)
cont

fit <- limma::lmFit(mval,design)
fit.cont <- limma::contrasts.fit(fit, contrasts = cont)
fit2 <- limma::eBayes(fit.cont)

testmat <- limma::decideTests(fit2)
summary(testmat)


#Make upset plot

## Hyper
#filter for all 0 rows, and rows with any negative, and NA rows
hypermat <- testmat@.Data
hypermat <- ifelse(hypermat == -1, 0,hypermat)

hypermat <- hypermat[!rowSums(hypermat) == 0 &
                       !is.na(rowSums(hypermat)),]

dim(hypermat) #7,123
head(hypermat)
colnames(hypermat) <- c("Met Vs NM", "CRC Vs NM", 
                        "Met Vs CRC")

upset_withnumbers <- function(m, color){
  col_size = comb_size(m)
  row_size = set_size(m)
  
  ht = UpSet(m, pt_size = unit(5, "mm"), lwd = 3, 
             #set_order = order_vec,
             top_annotation = 
               HeatmapAnnotation("No. probes" = 
                                   anno_barplot(comb_size(m), 
                                                border = FALSE, 
                                                gp = gpar(fill = color), 
                                                height = unit(7, "cm"),
                                                ylim = c(0, max(col_size)*1.1)
                                   )),
             right_annotation = upset_right_annotation(m,
                                                       width = unit(4, "cm"),
                                                       ylim = c(0, max(row_size)*1.1)))
  ht = draw(ht)
  
  col_od = column_order(ht)
  row_od = row_order(ht)
  
  decorate_annotation("No. probes", {
    grid.text(col_size[col_od], 
              seq_len(length(col_size)), 
              unit(col_size[col_od], "native") + unit(2, "mm"), 
              default.units = "native", just = "bottom",
              gp = gpar(fontsize = 8))
  })
  decorate_annotation("Set size", {
    grid.text(row_size[row_od], 
              unit(row_size[row_od], "native") + unit(2, "mm"), 
              rev(seq_len(length(row_size))), 
              default.units = "native", just = "bottom", rot = -90,
              gp = gpar(fontsize = 8))
  })
}

m <- make_comb_mat(hypermat)
upset_withnumbers(m, "#e1cf22")

## Hypo
#filter for all 0 rows, and rows with any negative, and NA rows
hypomat <- testmat@.Data
hypomat <- ifelse(hypomat == 1, 0,hypomat)
hypomat <- abs(hypomat)

hypomat <- hypomat[!rowSums(hypomat) == 0 &
                     !is.na(rowSums(hypomat)),]


dim(hypomat) #14,634
head(hypomat)
colnames(hypomat) <- c("Met Vs NM", "CRC Vs NM", 
                       "Met Vs CRC")
m <- make_comb_mat(hypomat)
upset_withnumbers(m, "#3b56d8")



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

#Transform beta vals
betas_trans <- asin(2*betas-1)
bad <- BiocGenerics::rowSums(is.finite(betas_trans)) < ncol(betas_trans)
if(any(bad)) betas_trans <- betas_trans[!bad,,drop=FALSE]

# Get mvals?
mval <- log2(betas/(1-betas))

# MDS
limma::plotMDS(mval, top = 1000, col=as.integer(factor(metadata$state)))

#Do limma with mvals
design <- model.matrix(~0+state + patient, data = metadata)
design
cont <- limma::makeContrasts(stateS-stateN, 
                             stateMD-stateN, 
                             stateS-stateMD,
                             levels = design)
cont

#fit <- limma::lmFit(mval,design)
fit <- limma::lmFit(betas,design)
fit.cont <- limma::contrasts.fit(fit, contrasts = cont)
fit2 <- limma::eBayes(fit.cont)
testmat <- limma::decideTests(fit2)
summary(testmat)


#Make upset plot

## Hyper
#filter for all 0 rows, and rows with any negative, and NA rows
hypermat <- testmat@.Data
hypermat <- ifelse(hypermat == -1, 0,hypermat)

hypermat <- hypermat[!rowSums(hypermat) == 0 &
                       !is.na(rowSums(hypermat)),]

dim(hypermat) #72,547
head(hypermat)
colnames(hypermat) <- c("Met Vs NM", "CRC Vs NM", 
                        "Met Vs CRC")

m <- make_comb_mat(hypermat)
upset_withnumbers(m, "#e1cf22")

## Hypo
#filter for all 0 rows, and rows with any negative, and NA rows
hypomat <- testmat@.Data
hypomat <- ifelse(hypomat == 1, 0,hypomat)
hypomat <- abs(hypomat)

hypomat <- hypomat[!rowSums(hypomat) == 0 &
                       !is.na(rowSums(hypomat)),]


dim(hypomat) #27,797
head(hypomat)
colnames(hypomat) <- c("Met Vs NM", "CRC Vs NM", 
                        "Met Vs CRC")
m <- make_comb_mat(hypomat)
upset_withnumbers(m, "#3b56d8")
