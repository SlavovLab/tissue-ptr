library(plyr)
load("Rdata/associations.RData")

protein <- as.matrix(read.csv("protein.csv",row.names=1))
tab.mart <- read.csv("mart_export.txt",header=TRUE,stringsAsFactors=FALSE)
tab.mart <- tab.mart[!duplicated(tab.mart[,1]),]
rownames(tab.mart) <- tab.mart[,1]
tab.mart <- tab.mart[rownames(protein), ]
go.names <- read.table("Go_names",comment.char="!",sep="\t",stringsAsFactors=FALSE,quote = "")
rownames(go.names) <- go.names[,1]

scientific_10 <- function(x) {
  parse(text=paste("10^", log10(x)))
}

sciNotation <- function(x, digits = 0) {
    if(x==0){
        return(expression("<" ~ 10^-15))
    }
    if(is.na(x)){
        warning("x is NA")
        return(1)
    }
    if (length(x) > 1) {
        return(append(sciNotation(x[1]), sciNotation(x[-1])))
    }
    if (!x) return(0)
    exponent <- floor(log10(x))
    base <- round(x / 10^exponent, digits)
    if(base==10){
        base <- 1
        exponent <- exponent-1
    }
    as.expression(substitute(base %*% 10^exponent, 
			list(base = base, exponent = exponent)))
}

e2n <- function(ensembl){
    tab.mart[match(toupper(ensembl),tab.mart$Ensembl.Gene.ID),"Associated.Gene.Name"]
}

n2e <- function(nm,first=TRUE){
    if(first){
        tab.mart[match(toupper(nm),tab.mart$Associated.Gene.Name),"Ensembl.Gene.ID"]
    }else{
        tab.mart[tab.mart$Associated.Gene.Name==toupper(nm),"Ensembl.Gene.ID"]
    }
}

e2r <- function(ensembl){
    tab.mart[match(toupper(ensembl),tab.mart$Ensembl.Gene.ID),"RefSeq.Protein.ID"]
}

r2e <- function(refseq){
    tab.mart[match(toupper(refseq),tab.mart$RefSeq.Protein.ID),"Ensembl.Gene.ID"]
}


getProteinsInGroup <- function(go_id,ensembl=FALSE){
    prots <- group2proteins[go_id]
    if(ensembl){ prots <- n2e(prots)}
    prots
}

getGroupsForProtein <- function(nm,ensembl=FALSE){
    if(ensembl){nm <- e2n(nm)}
    sapply(nm,function(x) unique(protein2groups[[x]]))
}


avg.colors <- function(col1,col2){
    rgb1 <- sapply(c(2,4,6),function(i) as.double(paste("0x",substr(col1,i,i+1),sep="")))
    rgb2 <- sapply(c(2,4,6),function(i) as.double(paste("0x",substr(col2,i,i+1),sep="")))

    rgb.avg <- round((rgb1+rgb2)/2)

    rgb(rgb.avg[1],rgb.avg[2],rgb.avg[3],maxColorValue=255)
}

strFormat <- function(x) {
  s <- strsplit(x, "\\.")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}

fisher.transform <- function(r){ 1/2*log((1+r)/(1-r)) }

getZscores <- function(mat1, mat2, min.pairs=4){


    ## Remove rows which show no variation or are all NA
    sd.mat1 <- apply(mat1, 1, function(x) sd(x,na.rm=TRUE))
    sd.mat2 <- apply(mat2, 1, function(x) sd(x,na.rm=TRUE))

    toRemove <- which(sd.mat1==0 | sd.mat2==0 | is.na(sd.mat1) | is.na(sd.mat2))
    if(length(toRemove) > 0) {
        sd.mat1 <- sd.mat1[-toRemove]
        sd.mat2 <- sd.mat2[-toRemove]
        mat1 <- mat1[-toRemove, ]
        mat2 <- mat2[-toRemove, ]
    }
    
    if(!all(dim(mat1)==dim(mat2)))
        stop("non-conformable arrays")

    n.pairwise <- rowSums(!is.na(mat1*mat2))
    
    cors <- sapply(which(n.pairwise >= min.pairs), function(i) {
        if(sd(mat1[i, ],na.rm=TRUE) == 0 | sd(mat2[i, ], na.rm=TRUE) == 0)
            NA
        else
            cor(mat1[i, ], mat2[i, ], use="pairwise.complete.obs")
    })

    n.pairwise <- n.pairwise[n.pairwise >= min.pairs]

    ft <- fisher.transform(cors)

    wts <- 1/(n.pairwise-3)
    wts <- wts/sum(wts)
    z <- sum(wts*ft)
    pop.cor <- (exp(2*z)-1)/(exp(2*z)+1)

    z.score <- (ft-z)*sqrt(n.pairwise-3)
    
    list(z.score=z.score,within.cors=cors,n.pairwise=n.pairwise,sd.mat1=sd.mat1,sd.mat2=sd.mat2,pop.cor=pop.cor)
    
}


##############################################
## GO correlation analysis!
##############################################

find_go_groups <- function(z.scores, fdr.thresh, within.cors, sd.mrnas, sd.prots) {
    glen <- sapply(group2proteins, function(x) length(unlist(x)))
    
    pvals <- rep(NA,length(allGroups))
    med.z <- rep(NA,length(allGroups))
    med.cor <- rep(NA,length(allGroups))
    med.sd.mrna <- med.sd.prot <- rep(NA,length(allGroups))

    names(pvals) <- names(med.z) <- names(med.cor) <- names(med.sd.mrna) <- names(med.sd.prot) <- allGroups

    prots.lst <- list()
    for(grp in allGroups) {
        prots <- intersect(n2e(group2proteins[[grp]]),names(z.scores))
        if(length(prots)>5 & length(prots) < 200 ) {

            prots.lst[[grp]] <- prots
            res <- ks.test(z.scores[prots],z.scores[setdiff(names(z.scores),prots)])
            pvals[grp] <- res$p.value
            med.z[grp] <- median(z.scores[prots],na.rm=TRUE)
            med.cor[grp] <- median(within.cors[prots],na.rm=TRUE)
            med.sd.mrna[grp] <- median(sd.mrnas[prots],na.rm=TRUE)
            med.sd.prot[grp] <- median(sd.prots[prots],na.rm=TRUE)

        }
    }
    pvals <- pvals[!is.na(pvals)]

    sorted.pvals <- sort(pvals)
    M <- length(sorted.pvals)


    qvals <- sorted.pvals/((1:M)/M)
    qvals <- qvals[qvals<fdr.thresh]
    
    med.z <- med.z[names(qvals)]
    med.cor <- med.cor[names(qvals)]
    med.sd.mrna <- med.sd.mrna[names(qvals)]
    med.sd.prot <- med.sd.prot[names(qvals)]

    nms <- sapply(names(qvals),function(x) substr(go.names[match(x,go.names[,1]),3],1,30))
    
    table <- data.frame(Name=nms,
                        "Median Z-Score"=med.z,
                        "Median Correlation"=med.cor,
                        "Median mRNA SD"=med.sd.mrna,
                        "Median protein SD"=med.sd.prot,
                        "Q-value"=qvals,"P-value"=sorted.pvals[names(nms)],
                        "Number of Genes"=glen[names(qvals)])
    
    colnames(table) <- c("Name","Median Z-Score",
                         "Median Correlation",
                         "Median mRNA SD","Median Protein SD",
                         "Q-value","P-value","Numer of Genes")

    table <- table[order(table[,"Median Z-Score"]),]
    table
}

## Create RData files with gene names GO association mappings
if(FALSE) {

    tab <- read.csv("mart_export.txt",header=TRUE,stringsAsFactors=FALSE)
    save(tab, file="gene_names.RData")

    assoc.tab <- read.table("gene_association.goa_human",
                            comment.char="!", sep="\t",
                            stringsAsFactors=FALSE, quote = "")
    protein <- as.matrix(read.csv("protein.csv", row.names=1))
    assoc.tab <- assoc.tab[assoc.tab[, 3] %in% e2n(rownames(protein)),]


    group2proteins <- list()
    protein2groups <- list()
    for(prot in rownames(protein)){
        nm <- e2n(prot)
        indices <- which(assoc.tab[,3]%in%nm)
        protein2groups[[nm]] <- unique(assoc.tab[indices,5])
        print(prot)
    }
    for(group in allGroups){
        indices <- which(assoc.tab[,5]%in%group)
        group2proteins[[group]] <- unique(assoc.tab[indices,3])
        print(group)
    }

    save(assoc.tab, group2proteins, protein2groups, file="associations.RData")

}
                                  



