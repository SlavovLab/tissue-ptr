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
illumina_orig <- illumina

## Aggregate over samples
illumina <- t(aggregate(t(illumina_orig), by=list(tissues), function(x) mean(x, na.rm=TRUE)))
colnames(illumina) <- illumina[1, ]
illumina <- illumina[-1, ]
colnames(illumina)[colnames(illumina)=="adrenal"] <- "adrenal.gland"
colnames(illumina)[colnames(illumina)=="thyroid"] <- "thyroid.gland"
colnames(illumina)[colnames(illumina)=="salivarygland"] <- "salivary.gland"
colnames(illumina)[colnames(illumina)=="bonemarrow"] <- "bone.marrow"
colnames(illumina)[colnames(illumina)=="smallintestine"] <- "small.intestine"
colnames(illumina)[colnames(illumina)=="urinarybladder"] <- "urinary.bladder"
colnames(illumina)[colnames(illumina)=="gallbladder"] <- "gall.bladder"

## split data into two sets to measure internal reliability
split1 <- unlist(sapply(unique(tissues), function(tissue) {

    len <- sum(tissues==tissue)
    (1:len) <= len/2

}))
split2 <- !split1

illumina_split1 <- illumina_orig[, split1]
illumina_split2 <- illumina_orig[, split2]
illumina_split1 <- t(aggregate(t(illumina_split1),
                               by=list(tissues[split1]),
                               function(x) mean(x, na.rm=TRUE)))
illumina_split2 <- t(aggregate(t(illumina_split2),
                               by=list(tissues[split2]),
                               function(x) mean(x, na.rm=TRUE)))
illumina_split1 <- illumina_split1[-1, ]
illumina_split2 <- illumina_split2[-1, ]
colnames(illumina_split1) <- colnames(illumina)
colnames(illumina_split2) <- colnames(illumina)

illumina <- illumina[, order(colnames(illumina))]
illumina_split1 <- illumina_split1[, order(colnames(illumina_split1))]
illumina_split2 <- illumina_split2[, order(colnames(illumina_split2))]

write.csv(illumina, "data/mrna_illumina.csv", quote=FALSE)
write.csv(illumina_split1, "data/mrna_illumina_split1.csv", quote=FALSE)
write.csv(illumina_split2, "data/mrna_illumina_split2.csv", quote=FALSE)

######################################
## ProteinAtlas data
######################################

pa <- read.table("~/Dropbox/NatureCommentData/ProteinAtlas/transcript_rna_tissue.tsv", sep="\t", header=TRUE)

pa_agg <- aggregate(pa[, -(1:2)], by=list(pa[, 1]), function(x) sum(x, na.rm=T))
pa_agg[pa_agg==0] <- NA
pa_agg[, -(1:2)] <- log10(pa_agg[, -(1:2)])

rnms <- pa_agg[, 1]
pa_agg <- pa_agg[, -1]
rownames(pa_agg) <- rnms

cnms <- colnames(pa_agg)
tissues <- strsplit(cnms, split="\\.|_", fixed=FALSE) %>%
    sapply(FUN = function(x) x[1])

## normalize mrna's against mrna[,1]
for(i in 2:ncol(pa_agg)){

    pa_agg[,i] <- pa_agg[,i] + median( pa_agg[,1] - pa_agg[,i], na.rm=TRUE)

}


## Aggregate over samples
pa_orig <- pa_agg
pa_agg <- t(aggregate(t(pa_agg), by=list(tissues), function(x) mean(x, na.rm=TRUE)))

colnames(pa_agg) <- pa_agg[1, ]
pa_agg <- pa_agg[-1, ]
colnames(pa_agg)[colnames(pa_agg)=="adrenal"] <- "adrenal.gland"
colnames(pa_agg)[colnames(pa_agg)=="thyroid"] <- "thyroid.gland"
colnames(pa_agg)[colnames(pa_agg)=="salivary"] <- "salivary.gland"
colnames(pa_agg)[colnames(pa_agg)=="bone"] <- "bone.marrow"
colnames(pa_agg)[colnames(pa_agg)=="urinary"] <- "urinary.bladder"
colnames(pa_agg)[colnames(pa_agg)=="small"] <- "small.intestine"
colnames(pa_agg)[colnames(pa_agg)=="lymph"] <- "lymph.node"
colnames(pa_agg)[colnames(pa_agg)=="gallbladder"] <- "gall.bladder"



## split data into two sets to measure internal reliability
split1 <- unlist(sapply(unique(tissues), function(tissue) {

    len <- sum(tissues==tissue)
    (1:len) <= len/2

}))
split2 <- !split1

pa_split1 <- pa_orig[, split1]
pa_split2 <- pa_orig[, split2]
pa_split1 <- t(aggregate(t(pa_split1),
                               by=list(tissues[split1]),
                               function(x) mean(x, na.rm=TRUE)))
pa_split2 <- t(aggregate(t(pa_split2),
                               by=list(tissues[split2]),
                               function(x) mean(x, na.rm=TRUE)))
pa_split1 <- pa_split1[-1, ]
pa_split2 <- pa_split2[-1, ]
colnames(pa_split1) <- colnames(pa_agg)
colnames(pa_split2) <- colnames(pa_agg)

pa_agg <- pa_agg[, order(colnames(pa_agg))]
pa_split1 <- pa_split1[, order(colnames(pa_split1))]
pa_split2 <- pa_split2[, order(colnames(pa_split2))]

write.csv(pa_split1, "data/mrna_pa_split1.csv", quote=FALSE)
write.csv(pa_split2, "data/mrna_pa_split2.csv", quote=FALSE)
write.csv(pa_agg, "data/mrna_pa.csv", quote=FALSE)

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
