rm(list=ls())
library(ggplot2)
library(scales)
library(corrplot)
library(fdrtool)
library(xtable)
library(Matrix)
source("util_functions.R")
load("Rdata/associations.RData")

proteinDataSource <- "Kim"
mrnaDataSource <- "illumina"

## --------------------------
proteinDataSource <- "Wilhelm"
mrnaDataSource <- "pa"

protein <- as.matrix(read.csv(sprintf("data/protein_%s.csv", proteinDataSource), row.names=1))
mrna <- as.matrix(read.csv(sprintf("data/mrna_%s.csv", mrnaDataSource),
                           row.names=1))

commonCols <- intersect(colnames(protein), colnames(mrna))
commonRows <- intersect(rownames(protein), rownames(mrna))
commonCols <- sort(commonCols)

protein <- protein[commonRows, commonCols]
mrna <- mrna[commonRows, commonCols]

########################################
## Normalize proteins and mrna
########################################

mrna[, 1] <- mrna[, 1] - median(mrna[,1], na.rm=TRUE)
protein[, 1] <- protein[, 1] - median(protein[,1], na.rm=TRUE)
for(i in 2:ncol(protein)){
    mrna[, i] <- mrna[,i] + median( mrna[,1] - mrna[,i], na.rm=TRUE)
    protein[, i] <- protein[,i] + median( protein[,1] - protein[,i], na.rm=TRUE)
}

tissue.ratios <- (10^protein)/(10^mrna)

cvs <- apply(tissue.ratios, 1, function(x) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE))

median.ratios <- apply(tissue.ratios, 1, function(x) median(x,na.rm=TRUE))

prediction.ratios <- log10((10^protein)/((10^mrna)*median.ratios))

protein.ratios <- t(apply(10^protein,1,function(x) x[combn(1:12,2)[1,]]/x[combn(1:12,2)[2,]]))
mrna.ratios <- t(apply(10^mrna,1,function(x) x[combn(1:12,2)[1,]]/x[combn(1:12,2)[2,]]))

predicted.raw <- (10^mrna)*median.ratios
tissues <- colnames(protein)
pretty.tissue.names <- sapply(tissues, strFormat)

allGroups <- unique(unlist(getGroupsForProtein(rownames(protein), ensembl=TRUE)))

## create rptr matrices
rptr.mat <- tissue.ratios
for(tissue in tissues){
    median.cur <- apply(tissue.ratios[, setdiff(tissues, tissue)], 1, function(x) median(x, na.rm=TRUE))
    rptr.mat[, tissue] <- log10(tissue.ratios[, tissue] / median.cur)
}

rptr.mat <- rptr.mat[apply(rptr.mat, 1, function(x) sum(is.na(x))) < (ncol(rptr.mat) - 1), ]

write.csv(rptr.mat,
          file=sprintf("data/rptr-%s-%s.csv", proteinDataSource, mrnaDataSource))

## normalized (by tissue) rptr matrices
transformed.mat <- rptr.mat
for(tissue in tissues){
    df <- ecdf(rptr.mat[, tissue])
    transformed.mat[, tissue] <- qnorm(df(rptr.mat[, tissue]))
}

write.csv(transformed.mat,
          file=sprintf("data/rptr-normalized-%s-%s.csv", proteinDataSource, mrnaDataSource))

test.type <- "ks"

for(tissue in tissues){
    print(sprintf("-------- %s --------",tissue))
    measured.prots <- rownames(transformed.mat)[!is.na(transformed.mat[,tissue])]
    pvals <- group.stat <- group.median <- group.mean <- group.median.protein <- rep(NA,length(allGroups))
    gene.counts <- numeric(length(allGroups))
    names(pvals) <- names(gene.counts) <- names(group.median) <- names(group.mean) <- names(group.stat) <- names(group.median.protein) <- allGroups
    count <- 0
    for(group in allGroups){
        prots <- intersect(n2e(group2proteins[[group]]),measured.prots)
        if(length(prots)>=1){

            if(test.type=="ks"){
                res <- ks.test(transformed.mat[prots, tissue],
                               transformed.mat[prots, setdiff(tissues,tissue)])
            } else if(test.type=="wilcox"){
                res <- wilcox.test(abs(transformed.mat[prots,setdiff(tissues,tissue)]),abs(transformed.mat[prots,tissue]),alternative="less")
            } else {
                stop("Must specify test type")
            }

            pvals[group] <- res$p.value
            group.stat[group] <- res$statistic
            group.median[group] <- median(rptr.mat[prots, tissue], na.rm=TRUE)
            group.mean[group] <- mean(rptr.mat[prots, tissue], na.rm=TRUE)
            if(is.na(group.mean[group]))
                browser()
            group.median.protein[group] <- median(protein[prots, tissue],na.rm=TRUE)
        }
        gene.counts[group] <- length(prots)
        count <- count+1
        if(count%%1000==0){
            print(count)
        }
    }

    save(pvals, gene.counts, group.stat, group.median, group.mean,
         file=sprintf("Rdata/%s-assoc-%s-%s-%s.RData", tissue, test.type,
                      proteinDataSource, mrnaDataSource))
}

################## Analysis and Figs ###############

significant.groups <- mean.groups <- matrix(nrow=length(allGroups),
                                            ncol=length(tissues),
                                            dimnames=list(allGroups, tissues))

for(tissue in tissues){
    load(sprintf("Rdata/%s-assoc-%s-%s-%s.RData", tissue, test.type,
                 proteinDataSource, mrnaDataSource))
    sorted.pvals <- sort(pvals)
    M <- length(sorted.pvals)
    qvals <- sorted.pvals/((1:M)/M)
    idx <- intersect(allGroups, names(qvals))
    significant.groups[idx, tissue] <- qvals[idx]
    mean.groups[idx, tissue] <- group.mean[idx]
}

tokeep <- apply(mean.groups, 1, function(x) !all(is.na(x)))
mean.groups <- mean.groups[tokeep,]
significant.groups <- significant.groups[tokeep,]

fdr.thresh <- 0.01
logical.mat <- significant.groups < fdr.thresh
logical.mat[is.na(logical.mat)] <- FALSE
hits <- sort(apply(logical.mat, 1, function(x) sum(x, na.rm=TRUE)),
             decreasing=TRUE)
hits <- hits[hits > 0]

clustered.groups <- (mean.groups[names(hits),])[hclust(dist(mean.groups[names(hits),]))$order,]
scales <- list(x=list(at=1:12,rot=90,lab=tissues),
               y=list(at=1:40, lab=sapply(
                                   go.names[rownames(clustered.groups),3],
                                   function(x) substring(x,1,30))))

unrefined.mat <- clustered.groups[rownames(clustered.groups),]

save(unrefined.mat, clustered.groups, transformed.mat,
     logical.mat, mean.groups, significant.groups,
     file=sprintf("matrix-%s-%s.RData", proteinDataSource, mrnaDataSource))

####

## pvals (deal with 0 and 1 pvals)
## pass into qnorm
## then run in locfdr

##################################
## Vertical mark plots
##################################


excluded.groups <- c("GO:0006120","GO:0005747","GO:0006120","GO:0005605",
                     "GO:0003954","GO:0071685","GO:0015078","GO:0022857",
                     "GO:0008137","GO:0015935","GO:0004129","GO:0019083",
                     "GO:0042776")

glen <- sapply(group2proteins,function(x) length(unlist(x)))
sub.grps <- names(glen[glen>=5 & glen<=100])

## These are manually picked from list of significant at 1%fdr 
toshow.list <- list()
toshow.list[["stomach"]] <- c("GO:0006415", "GO:0022904", "GO:0017119", "GO:0006418", "GO:0000049", "GO:0008536", "GO:0006418", "GO:0006811", "GO:0004364", "GO:0005840", "GO:0010827")
toshow.list[["kidney"]] <- c("GO:0022627", "GO:0006414", "GO:0022904", "GO:0048365", "GO:0006090", "GO:0031966", "GO:0003995", "GO:0006105", "GO:0006099", "GO:0016651")

for(tissue in tissues){
    load(sprintf("Rdata/%s-assoc-%s-%s-%s.RData", tissue, test.type,
                 proteinDataSource, mrnaDataSource))
    pdf(sprintf("Figs/within-marks-%s-%spercent-%s-%s-%s.pdf",
                tissue, round(100*fdr.thresh), test.type,
                proteinDataSource, mrnaDataSource))
    par(mar = c(0,14,4,8))
    pvals <- pvals[sub.grps]
    pvals <- pvals[setdiff(names(pvals),excluded.groups)]
    sorted.pvals <- sort(pvals)
    M <- length(sorted.pvals)
    qvals <- sorted.pvals/((1:M)/M)

    ## For FDR and pep calculations
    tmp.pvals <- sorted.pvals
    tmp.pvals <- tmp.pvals[!is.na(tmp.pvals)]
    tmp.pvals[tmp.pvals==0] <- 10e-6
    tmp.pvals[tmp.pvals==1] <- 0.99999
    
    pep <- fdrtool(tmp.pvals,statistic="pvalue",plot=FALSE)$lfdr
    names(pep) <- names(tmp.pvals)

    ## Only take groups that meet qval criteria
    qvals <- qvals[qvals < fdr.thresh]
    pep <- pep[qvals < fdr.thresh]
    print(go.names[names(qvals), 3, drop=FALSE])

    if(tissue %in% names(toshow.list)){
        grps <- intersect(toshow.list[[tissue]], names(qvals))
    } else {
        grps <- names(qvals)
        grps.counts <- sapply(grps,function(x) length(unlist(group2proteins[x])))
        grps <- cullGroups(grps,10,0)
        grps <- names(sort(abs(group.mean[grps]))[1:(min(10,length(grps)))])
    }
    
    grps <- names(sort(group.mean[grps]))
    
    vec <- rptr.mat[, tissue]
    vec <- vec[!is.na(vec)]
    sorted.vec <- sort(vec)
    col.grad <- colorRampPalette(c("blue","grey","red"))(100)
    layout(matrix(c(1,2),ncol=1),heights=c(1.3,2))
    plot.new()
    plot.window(xlim=c(1,length(sorted.vec)),ylim=c(-4,4),xpd=TRUE)
    col.all <- col.grad[findInterval(sorted.vec,seq(-1,1,length.out=99))+1]
    segments(1:length(sorted.vec),0,1:length(sorted.vec),sorted.vec,col=col.all,lwd=0.1)
    axis(side=2, at=seq(-4,4,2), labels=seq(-4,4,2),cex.axis=1.5,font=2)
    mtext("log rPTR ratio",side=2,line=3,font=2,cex=1.4)
    par(mar=c(4,14,2,8))
    plot.new()
    plot.window(xlim=c(1,length(sorted.vec)),ylim=c(0,20),xpd=TRUE)
    par(xpd=TRUE)
    text(1.25*length(sorted.vec),20,label="P-Value",adj = c(0.5,1),cex=1.5)
    count <- 17
    for(grp in grps){ 
        
        gns <- n2e(unlist(group2proteins[grp]))
        locs <- match(gns,names(sorted.vec))
        segments(locs,count,locs,count+1,col=col.grad[findInterval(sorted.vec[locs],seq(-1,1,length.out=99))+1])

        lab <- go.names[match(grp,go.names[,1]),3]
        if(nchar(lab)>30)
            lab <- gsub('^(.{20}\\S+)(.*)$', '\\1\n\\2', lab)
        text(x=-length(sorted.vec),y=count+.5,lab,cex=0.8,adj = c(0,0.5),font=2)
        text(x=1.25*length(sorted.vec),y=count+.5, lab=sciNotation(pvals[grp]), cex=1.5,adj = c(0.5,0.5), font=2)
        count <- count-2
    }
    dev.off()
} 

## select down to n groups by minimizing overlaps
## use greedy algo
cullGroups <- function(grps, n=10, overlap.thresh=0){

    overlap.mat <- matrix(0,nrow=length(grps),ncol=length(grps))
    colnames(overlap.mat) <- rownames(overlap.mat) <- grps
    
    for(grp in grps){
        gns <- n2e(unlist(group2proteins[grp]))
        for(grp2 in setdiff(grps,grp)){
            gns2 <- n2e(unlist(group2proteins[grp2]))
            score <- length(intersect(gns, gns2)) / min(length(gns), length(gns2))
            overlap.mat[grp,grp2] <- overlap.mat[grp2,grp] <- score
        }
    }

    while(nrow(overlap.mat)>n & max(as.vector(overlap.mat)) > overlap.thresh){
        maximum <- which.max(apply(overlap.mat,1,max))
        overlap.mat <- overlap.mat[-maximum,-maximum]
    }

    rownames(overlap.mat)
}

##############################################
## GO correlation analysis
## Find groups that have higher or lower median correlation
## than median of all genes 
##############################################

pvals <- rep(NA,length(allGroups))
med.z <- rep(NA,length(allGroups))
med.cor <- rep(NA,length(allGroups))
med.sd.mrna <- med.sd.prot <- rep(NA,length(allGroups))
                                  
names(pvals) <- names(med.z) <- names(med.cor) <- names(med.sd.mrna) <- names(med.sd.prot) <- allGroups


proteinConsensus <- as.matrix(read.csv(file="data/protein_consensus.csv",
                             row.names=1))
proteinConsensus <- proteinConsensus[rownames(mrna), colnames(mrna)]

mrnaProtInfo <- getZscores(proteinConsensus, mrna)
within.cors <- mrnaProtInfo$within.cors
sd.prots <- mrnaProtInfo$sd.mat1
sd.mrnas <- mrnaProtInfo$sd.mat2

for(grp in allGroups) {
    prots <- intersect(n2e(group2proteins[[grp]]), names(z.score))
    if(length(prots)>5 & length(prots) < 200 ) {
        res <- ks.test(z.score[prots], z.score[setdiff(names(z.score),prots)])
        pvals[grp] <- res$p.value
        med.z[grp] <- median(z.score[prots],na.rm=TRUE)
        med.cor[grp] <- median(within.cors[prots],na.rm=TRUE)
        med.sd.mrna[grp] <- median(sd.mrnas[prots],na.rm=TRUE)
        med.sd.prot[grp] <- median(sd.prots[prots],na.rm=TRUE)
    }
}

pvals <- pvals[!is.na(pvals)]
sorted.pvals <- sort(pvals)
M <- length(sorted.pvals)
out <- fdrtool(pvals,"pvalue")
qvals <- sorted.pvals / ((1:M)/M)
qvals <- qvals[qvals < fdr.thresh]

med.z <- med.z[names(qvals)]
med.cor <- med.cor[names(qvals)]
med.sd.mrna <- med.sd.mrna[names(qvals)]
med.sd.prot <- med.sd.prot[names(qvals)]

nms <- sapply(names(qvals),function(x) substr(go.names[match(x,go.names[,1]),3],1,30))

table <- data.frame(Name=nms,"Median Z-Score"=med.z,"Median Correlation"=med.cor,"Median mRNA SD"=med.sd.mrna,"Median protein SD"=med.sd.prot,"Q-value"=qvals,"P-value"=sorted.pvals[names(nms)],"Number of Genes"=glen[names(med.z)])
colnames(table) <- c("Name","Median Z-Score","Median Correlation","Median mRNA SD","Median Protein SD","Q-value","P-value","Numer of Genes")

print(
    xtable(
        table[order(table[,3]),c(1:3,6:8)],
        display=c("s","s","f","f","g","g","d"),
        digits=c(0,0,2,2,2,2,2)),
    file="~/Dropbox/NatureComment/corr_table.tex",
    floating.environment = 'sidewaystable',
    hline.after=c(0,nrow(table)-sum(med.cor>pop.cor))
)

