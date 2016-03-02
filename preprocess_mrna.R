rm(list=ls())
library(magrittr)

######################################
## Illumina data
######################################

illumina <- read.table("~/Dropbox/NatureCommentData/Illumina/Illumina_FPKM_allsamples.txt", header=TRUE, sep="\t", row.names=1)
illumina[illumina==0] <- NA
illumina <- log10(illumina)

cnms <- colnames(illumina)

tissues <- strsplit(cnms, split="\\.|_", fixed=FALSE) %>%
    sapply(FUN = function(x) x[2])

## normalize mrna's against mrna[,1]
for(i in 2:ncol(illumina)){

    illumina[,i] <- illumina[,i] + median( illumina[,1] - illumina[,i], na.rm=TRUE)

}

## Aggregate over samples
illumina <- t(aggregate(t(illumina), by=list(tissues), function(x) mean(x, na.rm=TRUE)))
illumina <- illumina[-1, ]
colnames(illumina) <- sort(unique(tissues))
colnames(illumina)[colnames(illumina)=="adrenal"] <- "adrenal.gland"
colnames(illumina)[colnames(illumina)=="thyroid"] <- "thyroid.gland"
colnames(illumina)[colnames(illumina)=="salivarygland"] <- "salivary.gland"
colnames(illumina)[colnames(illumina)=="bonemarrow"] <- "bone.marrow"
colnames(illumina)[colnames(illumina)=="smallintestine"] <- "small.intestine"
colnames(illumina)[colnames(illumina)=="urinarybladder"] <- "urinary.bladder"
colnames(illumina)[colnames(illumina)=="gallbladder"] <- "gall.bladder"

write.csv(illumina, "data/mrna_illumina.csv", quote=FALSE)

######################################
## MCP data
######################################

mcp <- read.csv("~/Dropbox/NatureCommentData/MCP/mcp.M113.csv", row.names=1)
mcp <- mcp[, -ncol(mcp)]
mcp[mcp==0] <- NA
mcp <- log10(mcp)

## normalize mrna's against mrna[,1]
for(i in 2:ncol(mcp)){

    mcp[,i] <- mcp[,i] + median( mcp[,1] - mcp[,i], na.rm=TRUE)

}

write.csv(mcp, "data/mrna_mcp.csv", quote=FALSE)

######################################
## ProteinAtlas data
######################################

pa <- read.table("~/Dropbox/NatureCommentData/ProteinAtlas/transcript_rna_tissue.tsv", sep="\t", header=TRUE)
pa[pa==0] <- NA
pa <- log10(pa)

pa_agg <- aggregate(pa[, -(1:2)], by=list(pa[, 1]), function(x) sum(x, na.rm=T))

rnms <- pa_agg[, 1]
pa_agg <- pa_agg[, -1]
rownames(pa_agg) <- rnms

cnms <- colnames(pa_agg)
tissues <- strsplit(cnms, split="\\.|_", fixed=FALSE) %>%
    sapply(FUN = function(x) x[1])

## normalize mrna's against mrna[,1]
for(i in 2:ncol(illumina)){

    pa_agg[,i] <- pa_agg[,i] + median( pa_agg[,1] - pa_agg[,i], na.rm=TRUE)

}

## Aggregate over samples
pa_agg <- t(aggregate(t(pa_agg), by=list(tissues), function(x) mean(x, na.rm=TRUE)))

pa_agg <- pa_agg[-1, ]
colnames(pa_agg) <- sort(unique(tissues))
colnames(pa_agg)[colnames(pa_agg)=="adrenal"] <- "adrenal.gland"
colnames(pa_agg)[colnames(pa_agg)=="thyroid"] <- "thyroid.gland"
colnames(pa_agg)[colnames(pa_agg)=="salivary"] <- "salivary.gland"
colnames(pa_agg)[colnames(pa_agg)=="bone"] <- "bone.marrow"
colnames(pa_agg)[colnames(pa_agg)=="urinary"] <- "urinary.bladder"
colnames(pa_agg)[colnames(pa_agg)=="small"] <- "small.intestine"
colnames(pa_agg)[colnames(pa_agg)=="lymph"] <- "lymph.node"
colnames(pa_agg)[colnames(pa_agg)=="gallbladder"] <- "gall.bladder"

write.csv(pa_agg, "data/mrna_pa.csv", quote=FALSE)
