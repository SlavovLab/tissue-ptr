source("util_functions.R")
library(plyr)

rawDataDir <- "~/Dropbox/tissue-ptr data/"

## data for the primary analysis comes from Wilhelm
protnm1 <- "Wilhelm"
protnm2 <- "Kim"

protein <- as.matrix(read.csv(sprintf("data/protein_%s.csv", protnm1), row.names=1))
protein2 <- as.matrix(read.csv(sprintf("data/protein_%s.csv", protnm2), row.names=1))
proteinConsensus <- as.matrix(read.csv(sprintf("data/protein_consensus.csv", protnm2), row.names=1))

tissues <- colnames(proteinConsensus)

wilhelmPepsAll <- c()
kimPepsAll <- c()

globalTab <- data.frame(Sequence=rep("", 0), Modifications=rep("", 0), Leading.Razor.Protein=rep("", 0), Gene.Names=rep("", 0))
countsMat <- c()
for(tissue in tissues){
    print(tissue)

    ## Protein counts
    if(tissue %in% colnames(protein)) {
        protWilhelm <- sum(!is.na(protein[, tissue]))
    } else {
        protWilhem <- 0
    }

    if(tissue %in% colnames(protein2)) {
        protKim <- sum(!is.na(protein2[, tissue]))
    } else {
        protKim <- 0
    }

    protConsensus <- sum(!is.na(proteinConsensus[, tissue]))

    pepfile1 <- sprintf("%sKim_Search_Results_Data/%s_evidence.txt",
                rawDataDir, tissue)
    if(file.exists(pepfile1)) {
        peptidesTab <- read.table(pepfile1, stringsAsFactors=FALSE,
                                  sep="\t", header=TRUE)
        kimPeps <- unique(subset(peptidesTab, Modifications=="Unmodified")$Sequence)

        peptidesAgg <- aggregate(peptidesTab[, "Intensity"], by=peptidesTab[, c("Sequence", "Modifications", "Leading.Razor.Protein", "Gene.Names")], function(x) sum(as.numeric(x), na.rm=TRUE))

        globalTab <- merge(globalTab, peptidesAgg, by=c("Sequence", "Modifications", "Leading.Razor.Protein", "Gene.Names"), all=TRUE)
        colnames(globalTab)[ncol(globalTab)] <- paste("kim", tissue, sep=".")

    } else {
        kimPeps <- c()
    }
    
    pepfile2 <- sprintf("%sWilhelm_Search_Results_Data/%s_evidence.txt",
                        rawDataDir, tissue)
    if(file.exists(pepfile2)) {
        peptidesTab <- read.table(pepfile2, stringsAsFactors=FALSE,
                              sep="\t", header=TRUE)
        wilhelmPeps <- unique(subset(peptidesTab, Modifications=="Unmodified")$Sequence)

        peptidesAgg <- aggregate(peptidesTab[, "Intensity"], by=peptidesTab[, c("Sequence", "Modifications", "Leading.Razor.Protein", "Gene.Names")], function(x) sum(as.numeric(x), na.rm=TRUE))

        globalTab <- merge(globalTab, peptidesAgg, by=c("Sequence", "Modifications", "Leading.Razor.Protein", "Gene.Names"), all=TRUE)
        colnames(globalTab)[ncol(globalTab)] <- paste("wilhelm", tissue, sep=".")
        
    } else {
        wilhelmPeps <- c()
    }
    
    wilhelmPepsAll <- c(wilhelmPepsAll, wilhelmPeps)
    kimPepsAll <- c(kimPepsAll, kimPeps)

    countsMat <- rbind(countsMat,
                       c(protWilhelm, length(wilhelmPeps),
                         protKim, length(kimPeps),
                         protConsensus, length(union(wilhelmPeps, kimPeps))))
    
}

write.table(globalTab, file="allPeptides.tab", quote=FALSE, row.names=FALSE, sep="\t")

countsMat <- rbind(countsMat,  c(nrow(protein),
                    length(unique(wilhelmPepsAll)),
                    nrow(protein2),
                    length(unique(kimPepsAll)),
                    length(union(rownames(protein), rownames(protein2))),
                    length(union(unique(wilhelmPepsAll),
                           unique(kimPepsAll)))
                    ))


rownames(countsMat) <- c(sapply(tissues, strFormat), "All")
colnames(countsMat) <- c("Wilhelm Prot", "Wilhelm Pep", "Kim Prot", "Kim pep",
                         "Consensus prot", "Consensus pep")
xtable(prettyNum(countsMat, big.mark=","))
print(xtable(countsMat, digits=0), format.args=list(big.mark = ","))
