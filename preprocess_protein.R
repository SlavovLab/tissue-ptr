rm(list=ls())
library(plyr)
source("util_functions.R")

protein_datasets <- c("Kim", "Wilhelm")
for(nm in protein_datasets) {

    protein.cur <- matrix(ncol=1)
    colnames(protein.cur) <- "Gene.names"
    tissues <- c("adrenal.gland", "esophagus", "testis", "kidney",
                 "ovary", "pancreas", "prostate", "liver", "uterus",
                 "stomach", "thryoid.gland", "salivary.gland", "spleen",
                 "colon", "heart", "lung")
    for(tissue in tissues){
        fname <- sprintf("~/Dropbox/NatureCommentData/%s_Search_Results_Data/%s_proteinGroups.txt",nm, tissue)
        if(file.exists(fname)) {
            tissue.tab <- read.table(fname,stringsAsFactors=FALSE,sep="\t",header=TRUE)
            tissue.tab <- tissue.tab[,c("Gene.names","iBAQ")]
            
            colnames(tissue.tab)[2] <- tissue
            tissue.tab <- tissue.tab[tissue.tab$Gene.names!="",]
            protein.cur <- merge(protein.cur,tissue.tab,by="Gene.names",all.x=TRUE,all.y=TRUE)
            protein.cur <- ddply(protein.cur,~Gene.names,function(x) colSums(x[,2:ncol(x),drop=FALSE],na.rm=TRUE))
        }
    }

    ens <- unlist(sapply(protein.cur[,"Gene.names"],function(x) {
        x
        
        split.lst <- strsplit(x,";")
        elst <- sapply(split.lst, n2e)

        if(all(is.na(elst))) {
            NA
        } else {
            if(length(elst)>1)
                NA
            else
                elst[!is.na(elst)]
        }
        
    }))
    not.na <- !duplicated(ens)&!is.na(ens)
    protein.cur <- protein.cur[not.na,]
    rownames(protein.cur) <- ens[not.na]
    protein.cur[protein.cur==0] <- NA
    protein.cur <- as.matrix(log10(protein.cur[,2:ncol(protein.cur)]))

    ## Normalize proteins against particular tissue
    norm.tissue <- colnames(protein.cur)[1]
    for(tissue in setdiff(colnames(protein.cur),norm.tissue)){

        protein.cur[,tissue] <- protein.cur[,tissue] +
            median(protein.cur[,norm.tissue] - protein.cur[,tissue],na.rm=TRUE)
    }

    write.csv(protein.cur,file=sprintf("data/protein_%s.csv", nm), quote=FALSE)

}
