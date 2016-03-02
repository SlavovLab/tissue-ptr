library(plyr)

protein <- as.matrix(read.csv("protein.csv",row.names=1))
tab.mart <- read.csv("mart_export.txt",header=TRUE,stringsAsFactors=FALSE)
tab.mart <- tab.mart[!duplicated(tab.mart[,1]),]
rownames(tab.mart) <- tab.mart[,1]
tab.mart <- tab.mart[rownames(protein), ]
go.names <- read.table("Go_names",comment.char="!",sep="\t",stringsAsFactors=FALSE,quote = "")
rownames(go.names) <- go.names[,1]

generate_csvs <- function(directory) {

######### Load mrna data #########################
    mrna.metadata <- read.csv("mrna_metadata.txt",sep="\t")

    mrna.all <- read.csv("FPKM_allsamples.txt",sep="\t",stringsAsFactors=FALSE)
    tissues <- c("kidney","liver","pancreas","prostate","stomach","adrenal","testis","thyroid","esophagus","ovary","salivarygland","spleen")

    split1 <- split2 <- mrna <- data.frame(gene=mrna.all$ensembl.gene.ID)

    for(tissue in tissues) {

        indices <- sample(grep(tissue,colnames(mrna.all)))
        nsubj <- length(indices)
        split1[,tissue] <- rowMeans(mrna.all[,indices[1:floor(nsubj/2)],drop=FALSE])
        split2[,tissue] <- rowMeans(mrna.all[,indices[(floor(nsubj/2)+1):nsubj],drop=FALSE])
        mrna[,tissue] <- rowMeans(mrna.all[,indices,drop=FALSE])
        
    }


    ## 0's should be NA
    split1[split1==0] <- NA
    split2[split2==0] <- NA
    mrna[mrna==0] <- NA

    ## remove gene names column and work with matrix of values
    rownames(split1) <- rownames(split2) <- rownames(mrna) <- mrna[,1]
    split1 <- split1[,2:ncol(split1)]
    split2 <- split2[,2:ncol(split2)]
    mrna <- mrna[,2:ncol(mrna)]

    ## Work with log data
    split1 <- log10(split1)
    split2 <- log10(split2)
    mrna <- log10(mrna)

    split1 <- as.matrix(split1)
    split2 <- as.matrix(split2)
    colnames(split1) <- colnames(split2) <- colnames(mrna) <- c("kidney","liver","pancreas","prostate","stomach","adrenal.gland","testis","thyroid","esophagus","ovary","salivary.gland","spleen")

    ## normalize mrna's against mrna[,1]
    for(i in 2:length(tissues)){

        split1[,i] <- split1[,i] + median(split1[,1]-split1[,i],na.rm=TRUE)
        split2[,i] <- split2[,i] + median(split2[,1]-split2[,i],na.rm=TRUE)
        mrna[,i] <- mrna[,i] + median(mrna[,1]-mrna[,i],na.rm=TRUE)
    }

    mrna.reliabilities <- sapply(1:nrow(split1),function(i) cor(split1[i,],split2[i,],use="pairwise.complete.obs"))
    names(mrna.reliabilities) <- rownames(split1)
    print(median(mrna.reliabilities,na.rm=TRUE))
    write.csv(mrna,file=sprintf("%s/mrna.csv",directory),quote=FALSE)
    write.csv(split1,file=sprintf("%s/mrna_split1.csv",directory),quote=FALSE)
    write.csv(split2,file=sprintf("%s/mrna_split2.csv",directory),quote=FALSE)

########### load protein data from Wilhelm et al and Kim et al #############

    protein_datasets <- c("Kim", "Wilhelm")
    for(nm in protein_datasets) {

        protein.cur <- matrix(ncol=1)
        colnames(protein.cur) <- "Gene.names"
        tissues <- c("adrenal.gland", "esophagus", "testis", "kidney", "ovary", "pancreas", "prostate",
                     "liver", "uterus", "stomach", "thryoid.gland", "salivary.gland", "spleen", "colon",
                     "heart", "lung")
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

            protein.cur[,tissue] <- protein.cur[,tissue]+median(protein.cur[,norm.tissue]-protein.cur[,tissue],na.rm=TRUE)
        }

        write.csv(protein.cur,file=sprintf("%s/protein_%s.csv", directory, nm),
                  quote=FALSE)

    }

}

scientific_10 <- function(x) {
  parse(text=paste("10^", log10(x)))
}

sciNotation <- function(x, digits = 0) {
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

    sd.mat1 <- sd.mat1[-toRemove]
    sd.mat2 <- sd.mat2[-toRemove]
    mat1 <- mat1[-toRemove, ]
    mat2 <- mat2[-toRemove, ]
    
    
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


                                  



