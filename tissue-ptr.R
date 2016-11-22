options(stringsAsFactors=FALSE)

library(ggplot2)
library(scales)
library(corrplot)
library(RColorBrewer)
library(mvtnorm)
library(fdrtool)
library(xtable)
library(grid)
library(fields)
source("util_functions.R")

## data for the primary analysis comes from Wilhelm
protnm1 <- "Wilhelm"
mrnanm1 <- "illumina"

protnm2 <- "Kim"
mrnanm2 <- "pa"

mrna <- as.matrix(read.csv(sprintf("data/mrna_%s.csv", mrnanm1), row.names=1))
mrna2 <- as.matrix(read.csv(sprintf("data/mrna_%s.csv", mrnanm2), row.names=1))
protein <- as.matrix(read.csv(sprintf("data/protein_%s.csv", protnm1), row.names=1))
protein2 <- as.matrix(read.csv(sprintf("data/protein_%s.csv", protnm2), row.names=1))

## Targeted MS dataset for validating quality of datasets
protVal <- read.csv("data/validation.csv", row.names=1)


tissue.names <- intersect(colnames(mrna), colnames(protein))
pretty.tissue.names <- sapply(tissue.names, strFormat)
n.tissues <- length(tissue.names)

gene.names <- intersect(rownames(mrna), rownames(protein))
n.genes <- length(gene.names)

mrna <- mrna[gene.names, tissue.names]
protein <- protein[gene.names, tissue.names]

tissue.ratios <- (10^protein) / (10^mrna)
median.ratios <- apply(tissue.ratios, 1, function(x) median(x, na.rm=TRUE))
prediction.ratios <- log10(tissue.ratios * median.ratios)
predicted.raw <- (10^mrna) * median.ratios

keep.indices <- which(log10(predicted.raw) != protein)
cor.raw <- cor(log10(predicted.raw[keep.indices]), as.numeric(protein[keep.indices]))

## number of proteins quantified (intersection and union) by tissue
commonProteins <- intersect(rownames(protein), rownames(protein2))
allProteins <- union(rownames(protein), rownames(protein2))

print(length(commonProteins))
print(length(allProteins))

abundInBoth <- sapply(intersect(colnames(protein), colnames(protein2)), function(tiss) {
    sum(!is.na(protein[commonProteins, tiss]) & !is.na(protein2[commonProteins, tiss]))
})

abundInEither <- sapply(intersect(colnames(protein), colnames(protein2)), function(tiss) {
    length(union(rownames(protein)[!is.na(protein[, tiss])], rownames(protein2)[!is.na(protein[, tiss])]))
})

print(cbind(abundInBoth, abundInEither))

#########################################################################
## 1a: fraction of total protein variance explained by scaled mRNA levels
#########################################################################

pdf("Figs/1a_mrna_v_protein.pdf", width=8, height=7)
print(

    ggplot(data=data.frame(x=predicted.raw[keep.indices], y=10^as.numeric(protein[keep.indices])),aes(x=x,y=y))+
    geom_point(colour="blue",alpha=0.03)+
    scale_x_log10(breaks=10^(3:11),label=scientific_10, limits=c(10^3, 10^11)) +
    scale_y_log10(breaks=10^(3:11),labels=scientific_10, limits=c(10^3, 10^11)) +
    labs(x="Scaled mRNA",y="Measured Protein",title="")+
    annotate("text", x = 10^4, y = 10^9,
             label=as.character(as.expression(
                 substitute(italic(R)[T]^2~"="~r2,
                            list(r2=formatC(cor.raw^2, format="f", digits=2))))),
             parse=TRUE, size=10) +
    annotate("text", x = 10^4, y = 10^10,
             label=as.character(as.expression(
                 substitute(italic(R)[T]~"="~r,
                            list(r=formatC(cor.raw, format="f", digits=2))))),
             parse=TRUE, size=10) +
    theme(axis.title.x = element_text(face="bold",size=30),
          axis.text.x  = element_text(size=30),
          axis.title.y = element_text(face="bold",size=30),
          axis.text.y  = element_text(size=30,hjust=0))
    
    )
dev.off()

#########################################################################
## 1b: simpson's paradox figure, manually pick genes to make point
#########################################################################

rsquared.vec <- c()
small.slopes <- c()
small.slopes.vals <- c()
for(i in 1:nrow(protein)){
    if( sum(!is.na(predicted.raw[i,]) & !is.na(protein[i,])) > 4){
        lm.fit <- lm(protein[i,] ~ log10(predicted.raw)[i,])
        slope <- lm.fit$coefficients[2]
        rsq <- summary(lm.fit)$r.squared
        if(rsq < 1 & sum(!is.na(predicted.raw[i, ]*protein[i,])) > 4)
            rsquared.vec <- c(rsquared.vec,rsq)
        if(slope < 0.2 & slope > -0.2 &
           max(log10(predicted.raw[i, ]), na.rm=TRUE) -
           min(log10(predicted.raw[i, ]), na.rm=TRUE) < 2) {

            small.slopes <- c(small.slopes, i)
            small.slopes.vals <- c(small.slopes.vals, slope)
        }
    }
}


## proteins which are measured in all tissues 
full.protein <- names(
    which(
        apply(protein[small.slopes, ], 1, function(x)
            sum(!is.na(x))) ==n.tissues &
        apply(predicted.raw[small.slopes, ], 1, function(x)
            sum(!is.na(x)))==n.tissues))

focus.prot <- "ENSG00000163631"
    
mrna.vec <- as.numeric(t(log10(predicted.raw)))
protein.vec <- as.numeric(t(protein))
id <- rep(rownames(mrna), each=n.tissues)
df <- data.frame(x=mrna.vec, y=protein.vec, id=id)

mean.prot <- rowMeans(protein[small.slopes, ], na.rm=TRUE)
qts <- quantile(mean.prot, c(0.05, 0.95))
## manually add in some proteins with high abundances
large.prots <- sample(names(mean.prot[mean.prot > qts[2]]), 10)
small.prots <- sample(names(mean.prot[mean.prot < qts[1]]), 10)

## sample proteins 
sub.proteins <- unique(c(large.prots, small.prots,
                         sample(unique(id)[small.slopes], 80)))
       
df.sub <- df[df$id %in% sub.proteins,]
cor.segs <- cor(df.sub$x, df.sub$y, use="complete.obs")


pdf("Figs/1b_mrna_v_protein_segs_horiz.pdf", width=8, height=7)
print(
    ggplot(data=df.sub, aes(x=x, y=y)) +
    geom_smooth(method=lm, se=FALSE, aes(group=id), colour=alpha("black",0.5)) +
    scale_x_continuous(breaks=3:11, label=function(x) parse(text=paste("10^", x)), limits=c(3,11)) +
    scale_y_continuous(breaks=3:11,labels=function(x) parse(text=paste("10^", x)), limits=c(3,11)) +
    labs(x="Scaled mRNA",y="Measured Protein")+
    annotate("text", x = 4, y = 9,
             label=as.character(as.expression(substitute(
                 italic(R)[T]^2~"="~r2,
                 list(r2=formatC(cor.segs^2, format="f", digits=2))))),
             parse=TRUE,size=10) +
    annotate("text", x = 4, y = 10,
             label=as.character(as.expression(substitute(
                 italic(R)[T]~"="~r,
                 list(r=formatC(cor.segs, format="f", digits=2))))),
             parse=TRUE,size=10)+
    theme(axis.title.x = element_text(face="bold",size=30),
          axis.text.x  = element_text(size=30),
          axis.title.y = element_text(face="bold",size=30),
          axis.text.y  = element_text(size=30,hjust=0)) +
    geom_point(data=data.frame(x=log10(predicted.raw[focus.prot,]),y=protein[focus.prot,],Tissue=tissue.names),aes(x=x,y=y,label=length(tissue.names):1,color=Tissue),size=3)+
    theme(legend.text=element_text(size=15),legend.title.align=0.5,
          legend.background=element_rect(fill=alpha("gray", 0.75)),legend.key=element_blank(),legend.title=element_text(size=14)) +
    guides(colour = guide_legend(override.aes = list(size=5)))+
    geom_smooth(data=data.frame(x=log10(predicted.raw[focus.prot,]),y=protein[focus.prot,]),method=lm,se=FALSE,colour="black",size=1.5)+theme(legend.justification=c(1,0), legend.position=c(1,0)),
)
dev.off()

within.indices <- which(apply(!is.na(protein * predicted.raw), 1, sum) > 3)
within.cors <- sapply(within.indices,function(i) cor(protein[i,], log(predicted.raw)[i,], use="pairwise.complete.obs"))
between.cors <- sapply(1:ncol(protein),function(i) cor(protein[,i],log10(predicted.raw)[,i],use="pairwise.complete.obs"))
between.cors.raw <- sapply(1:ncol(protein),function(i) cor(protein[,i],mrna[,i],use="pairwise.complete.obs"))

#########################################################################
## Within and between gene correlations
#########################################################################

pdf("Figs/appendix_within_hist.pdf")
ggplot(data.frame(x=within.cors))+geom_histogram(aes(x=x,fill="within"),colour="black",binwidth=0.05,size=0.1)+labs(x=expression(paste("Correlation, ", R[P])),y="",title="Across Tissue, Empirical")+
    theme(title = element_text(size=20),
          axis.title.x = element_text(face="bold",size=25),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16,hjust=0)) +
    scale_y_continuous(limits=c(0,300))+
    guides(fill=FALSE)
dev.off()

ggplot(data.frame(x=within.cors))+geom_histogram(aes(x=x,fill="within"),colour="black",binwidth=0.05,size=0.1)+labs(x=expression(paste("Correlation, ", R[P])),y="",title="Across Tissue, Empirical")+
    theme(title = element_text(size=20),
          axis.title.x = element_text(face="bold",size=25),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16,hjust=0)) +
    scale_y_continuous(limits=c(0,300))+
    guides(fill=FALSE)

pdf("Figs/appendix_between_hists.pdf")
ggplot(data.frame(x=between.cors, x2=between.cors.raw, labels=tissue.names)) + 
    geom_histogram(aes(x=x, fill="Scaled mRNA Levels"), binwidth=0.03, size=0.1) + 
    geom_histogram(aes(x=x2, fill="Measured mRNA Levels"), binwidth=0.03, size=0.1) + 
    scale_x_continuous(limits=c(-1, 1)) + 
    scale_y_continuous(limits=c(0, 5)) + 
    labs(x=expression(paste("Correlation, ", R[T])), y="", 
         title="Within Tissue,  Empirical") +
    theme(title = element_text(size=20),
          axis.title.x = element_text(face="bold",size=25),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16,hjust=0),
          legend.text=element_text(size = 16, face = "bold"),
          legend.justification=c(1,1), legend.position=c(1,1),
          legend.background = element_rect(fill = "transparent")) +
    guides(fill=guide_legend(title=NULL))
dev.off()

#########################################################################
### Histogram of R-squared
#########################################################################

pdf("Figs/histogram-rsq.pdf")
rsq.list <- data.frame(rsq=rsquared.vec)
ggplot(rsq.list,aes(x=rsq))+geom_histogram(binwidth=.05, colour="black", fill="grey")+labs(x="R-Squared",y="Count",title="Across Tissue, Empirical")+
    theme(axis.title.x = element_text(face="bold",size=30),
          axis.text.x  = element_text(size=16),
          axis.title.y = element_text(face="bold",size=30),
          axis.text.y  = element_text(size=16,hjust=0))
dev.off()


mrnaProtInfo <- getZscores(protein, mrna)
pop.cor <- mrnaProtInfo$pop.cor

print(sprintf("Estimated population correlation: %.3f", pop.cor))

## Q- Value analysis
z.score<- mrnaProtInfo$z.score
fd.out <- fdrtool(z.score,plot=FALSE)

head(sort(fd.out$qval),n=20)
print(sprintf("Number of q-values < 0.01: %i",sum(fd.out$qval<0.01)))
print(sprintf("Number of q-values < 0.1: %i",sum(fd.out$qval<0.1)))

## simulate correlations
sim.within.cors <- sapply(1:length(z.score), function(i) {
    cor(rmvnorm(n=mrnaProtInfo$n.pairwise[i], 
                mean=c(0,0), 
                sigma=matrix(c(1,pop.cor,pop.cor,1),nrow=2)))[1,2] 
})

pdf("Figs/appendix_sim_within_hist.pdf")
ggplot(data.frame(x=sim.within.cors))+geom_histogram(aes(x=x,fill="within"),colour="black",binwidth=0.05,size=0.1)+labs(x=expression(paste("Correlation, ", R[P])),y="",title=expression(paste("Across Tissue, Simulated (",rho,"=0.35)",sep="")))+
    theme(title = element_text(size=20),
          axis.title.x = element_text(face="bold",size=25),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16,hjust=0)) +
    scale_y_continuous(limits=c(0,300))+
    guides(fill=FALSE)
dev.off()


pdf("Figs/sim_and_within_cors.pdf")
pooled.cors <- c(sim.within.cors, within.cors)
label <- factor(rep(c("Simulated", "Empirical"), 
                    levels=c("Simulated", "Empirical"), 
                    each=length(within.cors)), ordered=TRUE)

ggplot(data.frame(cors=pooled.cors, cut=label), aes(x=cors, fill=label)) + 
    geom_bar(position=position_dodge(width=0.05), binwidth=0.1, size=0,
             aes(order=c(rep(c(1,2), each=length(within.cors))))) +
    labs(x=expression(paste("Correlation, ", R[P])), 
         y="", title=expression("Across Tissue")) +
    theme(title = element_text(size=20),
          axis.title.x = element_text(face="bold",size=25),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16,hjust=0),
          legend.title=element_blank(),
          legend.text = element_text(size = 16, face = "bold"),
          legend.justification=c(1,1), legend.position=c(1,1),
          legend.background = element_rect(fill = "transparent")) +
    scale_y_continuous(limits=c(0,500)) +
    geom_vline(xintercept = mrnaProtInfo$pop.cor, linetype=2)
dev.off()

#########################################################################
## Fold error / reliability simulation
#########################################################################

n.prot <- rowSums(!is.na(protein))
sd.prot.err <- uniroot(function(x) pnorm(log10(1.5),mean=0,sd=x)-0.75,lower=0,upper=10)$root
sd.prot.tot <- sqrt(sum((n.prot-1)*apply(protein,1,function(x) var(x,na.rm=TRUE)),na.rm=TRUE)/(sum(n.prot,na.rm=TRUE)-sum(n.prot>1)))
#sd.prot.tot <- median(apply(protein,1,function(x) sd(x,na.rm=TRUE)),na.rm=TRUE)
Rprot <- 1-sd.prot.err^2/sd.prot.tot^2

n.mrna <- rowSums(!is.na(mrna))
sd.mrna.err <- uniroot(function(x) pnorm(log10(1.2),mean=0,sd=x)-0.75,lower=0,upper=10)$root
## pooled estimator
sd.mrna.tot <- sqrt(sum((n.mrna-1)*apply(mrna,1,function(x) var(x,na.rm=TRUE)),na.rm=TRUE)/(sum(n.mrna,na.rm=TRUE)-sum(n.mrna>1)))
Rmrna <- 1-sd.mrna.err^2/sd.mrna.tot^2

## 
mar2 <- par()$mar
mar2[2] <- 5.1
par(mar=mar2)

curve(sapply(x,function(fc){ 1-(uniroot(function(s) pnorm(log10(fc),mean=0,sd=s)-0.75,lower=0,upper=1000)$root)^2/sd.prot.tot^2}),from=1.00,to=1.5,
      ylab="Reliability",xlab="Median Fold Error (%)",
      lwd=3,main="",
      cex.lab=2,ylim=c(0,1),col="red",xaxt="n",cex.axis=2,cex.main=2)

grid()
text(1.37,0.6,"mRNA",col="blue",cex=1.5)
text(1.40,0.9,"protein",col="red",cex=1.5)
axis(1,at=seq(1,1.5,by=0.1),labels=c(0,10,20,30,40,50),cex.axis=2)

segments(1.2,0,1.2,Rmrna,col="dark grey",lwd=3)
segments(0,Rmrna,1.2,Rmrna,col="dark grey",lwd=3)
segments(1.5,0,1.5,Rprot,col="dark grey",lwd=3)
segments(0,Rprot,1.5,Rprot,col="dark grey",lwd=3)

curve(sapply(x,function(fc){ 1-(uniroot(function(s) pnorm(log10(fc),mean=0,sd=s)-0.75,lower=0,upper=1000)$root)^2/sd.mrna.tot^2}),from=1.00,to=1.5,cex.lab=2,lwd=3,add=TRUE,col="blue")

dev.off()

#########################################################################
## R-squared as a function of reliabilities
## assuming true mrna/prot correlation = pop.cor
#########################################################################

mrnaProtInfo <- getZscores(protein, mrna)
pop.cor <- mrnaProtInfo$pop.cor

corrected.cor <- pop.cor/(Rprot*Rmrna)
corrected.rsq <- corrected.cor^2

## Compute mrna reliabilities
mrnaRows <- intersect(rownames(mrna), rownames(mrna2))
mrnaCols <- intersect(colnames(mrna), colnames(mrna2))

mrnaRelInfo <- getZscores(mrna[mrnaRows, mrnaCols], mrna2[mrnaRows, mrnaCols])
mrna.reliabilities <- mrnaRelInfo$within.cors
mrnaRel <- mrnaRelInfo$pop.cor

## Compute protein reliabilities
pcols <- intersect(colnames(protein), colnames(protein2))
prows <- intersect(rownames(protein), rownames(protein2))
cts <- apply(protein[prows, pcols] * protein2[prows, pcols], 1,
             function(x) {sum(!is.na(x))})

ppScores <- getZscores(protein[prows, pcols], protein2[prows, pcols])
protein.reliabilities <- ppScores$within.cors
protRel <- ppScores$pop.cor

#########################################################################
## Histogram of across tissue mRNA reliabilities
#########################################################################

pdf("Figs/mrna_reliabilities.pdf")
ggplot(data.frame(x=mrna.reliabilities))+geom_histogram(aes(x=x,fill="mRNA"),colour="black",binwidth=0.05,size=0.1)+labs(x="Reliability",y="",title=expression(paste("mRNA Reliabilities, Median = ", med.mrna, sep="")))+
    theme(title = element_text(size=20),
          axis.title.x = element_text(face="bold",size=25),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16,hjust=0)) +
    scale_x_continuous(limits=c(-1,1))+
    guides(fill=FALSE)
dev.off()

#########################################################################
## Histogram of across tissue protein reliabilities
#########################################################################

pdf("Figs/protein_reliabilities.pdf")
ggplot(data.frame(x=protein.reliabilities))+geom_histogram(aes(x=x,fill="mRNA"),colour="black",binwidth=0.05,size=0.1)+labs(x="Reliability",y="",title=expression(paste("Protein Reliabilities, Median = ",med.prot,sep="")))+
    theme(title = element_text(size=20),
          axis.title.x = element_text(face="bold",size=25),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16,hjust=0)) +
    scale_x_continuous(limits=c(-1,1))+
    guides(fill=FALSE)
dev.off()

#########################################################################
## ### Compare Kim and Wilhelm and create consensus dataset
#########################################################################

pcols <- Reduce(intersect,
                list(colnames(protein),
                     colnames(protein2),
                     colnames(mrna),
                     colnames(mrna2)))
prows <- Reduce(intersect,
                list(rownames(protein),
                     rownames(protein2),
                     rownames(mrna),
                     rownames(mrna2)))

## To facilitate fair comparison, only look at genes 
## measured in both data sets.
NAMask <- protein[prows, pcols] * protein2[prows, pcols]
NAMask[!is.na(NAMask)] <- 1

p1m1scores <- getZscores(protein[prows, pcols]*NAMask, mrna[prows ,pcols])
p2m1scores <- getZscores(protein2[prows, pcols]*NAMask, mrna[prows, pcols])
p1m2scores <- getZscores(protein[prows, pcols]*NAMask, mrna2[prows ,pcols])
p2m2scores <- getZscores(protein2[prows, pcols]*NAMask, mrna2[prows, pcols])

p1m1cor <- p1m1scores$pop.cor
p2m1cor <- p2m1scores$pop.cor
p1m2cor <- p1m2scores$pop.cor
p2m2cor <- p2m2scores$pop.cor

mean(!is.na(protein[prows, pcols]*mrna[prows, pcols]))
mean(!is.na(protein[prows, pcols]*mrna2[prows, pcols]))
mean(!is.na(protein2[prows, pcols]*mrna[prows, pcols]))
mean(!is.na(protein2[prows, pcols]*mrna2[prows, pcols]))

protConsensus <- protein2[prows, pcols]
for(i in 1:length(pcols)) {
    notNA <- which(!is.na(protein[prows, pcols[i]] * protein2[prows, pcols[i]]))
    p1NA <- which(is.na(protein[prows, pcols[i]]) & !is.na(protein2[prows, pcols[i]]))
    p2NA <- which(!is.na(protein[prows, pcols[i]]) & is.na(protein2[prows, pcols[i]]))                  

    
    p1m1 <- cor(protein[prows[notNA], pcols[i]], mrna[prows[notNA], pcols[i]], use="pairwise.complete.obs")
    p2m1 <- cor(protein2[prows[notNA], pcols[i]], mrna[prows[notNA], pcols[i]], use="pairwise.complete.obs")
    p1p2 <- cor(protein[prows, pcols[i]], protein2[prows, pcols[i]], use="pairwise.complete.obs")

    p1rel1 <- p1p2*p1m1/p2m1
    p2rel1 <- p1p2*p2m1/p1m1
    
    p1m2 <- cor(protein[prows, pcols[i]], mrna2[prows, pcols[i]], use="pairwise.complete.obs")
    p2m2 <- cor(protein2[prows, pcols[i]], mrna2[prows, pcols[i]], use="pairwise.complete.obs")
    p1rel2 <- p1p2*p1m2/p2m2
    p2rel2 <- p1p2*p2m2/p1m2

    ## average independent reliability estimates
    p1rel <- mean(c(p1rel1, p1rel2))
    p2rel <- mean(c(p2rel1, p2rel2))

    w1 <- 1 / (1 + 1/(p1rel/p2rel))
    w2 <- 1 - w1

    print(sprintf("%s weights: %f, %f", pcols[i], w1, w2))

    ## if measurements observed in both value is weighted average
    protConsensus[prows[notNA], pcols[i]] <-
        w1 * protein[prows[notNA], pcols[i]] +
        w2 * protein2[prows[notNA], pcols[i]]

    ## if observed in one but not the other
    protConsensus[prows[p1NA], pcols[i]] <- protein2[prows[p1NA], pcols[i]]
    protConsensus[prows[p2NA], pcols[i]] <- protein[prows[p2NA], pcols[i]]
}

## TODO
pcm1scores <- getZscores(protConsensus[prows, pcols]*NAMask, mrna[prows, pcols])
pcm1cor <- pcm1scores$pop.cor
pcm2scores <- getZscores(protConsensus[prows, pcols]*NAMask, mrna2[prows, pcols])
pcm2cor <- pcm2scores$pop.cor


## consensus dataset more correlated with mrna datasets
## than either the kim or wilhelm data individually:
cat(sprintf("\nConsensus / %s cor: %f,\n%s / %s cor: %f,\n%s / %s cor: %f\n",
              mrnanm1, pcm1cor,
              protnm1, mrnanm1, p1m1cor,
            protnm2, mrnanm1, p2m1cor))
cat(sprintf("\nConsensus / %s cor: %f,\n%s / %s cor: %f,\n%s / %s cor: %f\n",
              mrnanm2, pcm2cor,
              protnm1, mrnanm2, p1m2cor,
              protnm2, mrnanm2, p2m2cor))


#########################################################################
## ### Compare Kim, Wilhelm and consensus vs empirical prot rel
#########################################################################

prelCuts <- cut(protein.reliabilities,
                breaks=quantile(protein.reliabilities, seq(0, 1, length.out=6)))

medianPrel <- aggregate(protein.reliabilities,
                        list(lev=prelCuts), function(x) median(x, na.rm=TRUE))$x
p1m1relcor <- aggregate(p1m1scores$within.cors[names(protein.reliabilities)],
                 list(lev=prelCuts), function(x) median(x, na.rm=TRUE))$x
p2m1relcor <- aggregate(p2m1scores$within.cors[names(protein.reliabilities)],
          list(lev=prelCuts), function(x) median(x, na.rm=TRUE))$x
pcm1relcor <- aggregate(pcm1scores$within.cors[names(protein.reliabilities)],
                        list(lev=prelCuts), function(x) median(x, na.rm=TRUE))$x
p1m2relcor <- aggregate(p1m2scores$within.cors[names(protein.reliabilities)],
                 list(lev=prelCuts), function(x) median(x, na.rm=TRUE))$x
p2m2relcor <- aggregate(p2m2scores$within.cors[names(protein.reliabilities)],
          list(lev=prelCuts), function(x) median(x, na.rm=TRUE))$x
pcm2relcor <- aggregate(pcm2scores$within.cors[names(protein.reliabilities)],
                        list(lev=prelCuts), function(x) median(x, na.rm=TRUE))$x

pdf("Figs/RelVsCor_fagerberg.pdf", width=5, height=5)
plot(x=medianPrel, y=p1m1relcor, type="l", ylim=c(0,0.5), lwd=2,
     xlab="Median Protein/Protein Correlation", ylab="Median mRNA/Protein Correlation",
     main="Median across-tissue correlations\nversus protein reliability  (Fagerberg et al)", cex.axis=1.2, cex.lab=1.3, col="blue")
lines(x=medianPrel, y=p2m1relcor, type="l", ylim=c(0,1), lwd=2, col="black")
lines(x=medianPrel, y=pcm1relcor, type="l", ylim=c(0,1), lwd=2, col="red")
legend("topleft", legend=c("Consensus", "Wilhelm et al.", "Kim et al."), col=c("red", "blue", "black"), lty=1, lwd=2, box.lty=0)
dev.off()

pdf("Figs/RelVsCor_encode.pdf", width=5, height=5)
plot(x=medianPrel, y=p1m2relcor, type="l", ylim=c(0,0.5), lwd=2,
     xlab="Median Protein/Protein Correlation", ylab="Median mRNA/Protein Correlation",
     main="Median across-tissue correlations\nversus protein reliability  (ENCODE)", cex.axis=1.2, cex.lab=1.3, col="blue")
lines(x=medianPrel, y=p2m2relcor, type="l", ylim=c(0,1), lwd=2, col="black")
lines(x=medianPrel, y=pcm2relcor, type="l", ylim=c(0,1), lwd=2, col="red")
legend("topleft", legend=c("Consensus", "Wilhelm et al.", "Kim et al."), col=c("red", "blue", "black"), lty=1, lwd=2, box.lty=0)
dev.off()

#######################################################
## Write complete consensus dataset to file, 
## include genes or tissues missing in one datset but not the other
#######################################################

protConsensus <- rbind(protConsensus,
                       protein[setdiff(rownames(protein), prows), pcols])
protConsensus <- rbind(protConsensus,
                       protein2[setdiff(rownames(protein2),
                                        rownames(protConsensus)), pcols])
otherCols <- setdiff(union(colnames(protein), colnames(protein2)),
                     pcols)

## Add tissues that were in one dataset but not the other
## for data set with more coverage
protConsensus <- cbind(protConsensus, matrix(NA,
                                             nrow=nrow(protConsensus),
                                             ncol=length(otherCols),
                                             dimnames=list(
                                                 rownames(protConsensus),
                                                 otherCols)))
protConsensus[rownames(protein), setdiff(colnames(protein), pcols)] <-
    protein[rownames(protein), setdiff(colnames(protein), pcols)]
protConsensus[rownames(protein2), setdiff(colnames(protein2), pcols)] <-
    protein2[rownames(protein2), setdiff(colnames(protein2), pcols)]


write.csv(protConsensus, file="data/protein_consensus.csv")

#######################################################
## Compute "MSE" relative to validation data.
## Use gene / tissue pairs that are observed in all data
#######################################################

pvCols <- Reduce(intersect,
                 list(colnames(protVal), colnames(protein), colnames(protein2)))
pvRows <- Reduce(intersect,
                 list(rownames(protVal), rownames(protein), rownames(protein2)))

## Find prot-tissue measurements that are observed in all datasets
comparisonGenes <- which(!is.na(protein[pvRows, pvCols]) & !is.na(protein2[pvRows, pvCols]) & !is.na(protVal[pvRows, pvCols]))

## compute differences after removing mean differences
diffP1 <- scale(
    as.vector(protein[pvRows, pvCols][comparisonGenes]) -
    as.vector(unlist(protVal[pvRows, pvCols]))[comparisonGenes],
    scale=FALSE)
    
msep1 <- sum(diffP1^2, na.rm=TRUE) / length(diffP1)


diffP2 <- scale(
    as.vector(protein2[pvRows, pvCols][comparisonGenes]) -
    as.vector(unlist(protVal[pvRows, pvCols]))[comparisonGenes],
    scale=FALSE)
    
msep2 <- sum(diffP2^2, na.rm=TRUE) / length(diffP2)

diffPC <- scale(
    as.vector(protConsensus[pvRows, pvCols][comparisonGenes]) -
    as.vector(unlist(protVal[pvRows, pvCols]))[comparisonGenes],
    scale=FALSE)
mseConsensus <- sum(diffPC^2, na.rm=TRUE) / length(diffPC)

errorMat <- rbind(c(msep1, msep2, mseConsensus),
                  c(median(abs(diffP1)),
                    median(abs(diffP2)),
                    median(abs(diffPC))),
                  c(sum(!is.na(protein)),
                    sum(!is.na(protein2)),
                    sum(!is.na(protConsensus))))

colnames(errorMat) <- c(protnm1, protnm2, "Consensus")
rownames(errorMat) <- c("MSE", "MAD", "Coverage")
xtable(errorMat)

## Show improvement is significant

## Bootstrap interval for MSE improvement over Kim dataset
quantile(sapply(1:1000, function(x) {
    idx <- sample(length(diffPC), replace=TRUE)
    mean(diffPC[idx]^2) - mean(diffP2[idx]^2)
}), c(0.025, 0.975))

## Bootstrap interval for MAD improvement over Kim dataset
quantile(sapply(1:1000, function(x) {
    idx <- sample(length(diffPC), replace=TRUE)
    median(abs(diffPC[idx])) - median(abs(diffP2[idx]))
}), c(0.025, 0.975))


## Compare MSE by tissue
for(tis in pvCols) {
    print(tis)
    cgTis <- which(!is.na(protein[pvRows, tis]) & !is.na(protein2[pvRows, tis]) & !is.na(protVal[pvRows, tis]))

    diffP2 <- scale(
        as.vector(protein2[pvRows, tis][cgTis]) -
        as.vector(unlist(protVal[pvRows, tis]))[cgTis],
        scale=FALSE)
    
    print(sum(diffP2^2, na.rm=TRUE) / length(diffP2))

    diffPC <- scale(
        as.vector(protConsensus[pvRows, tis][cgTis]) -
        as.vector(unlist(protVal[pvRows, tis]))[cgTis],
        scale=FALSE)
    print(sum(diffPC^2, na.rm=TRUE) / length(diffPC))
    
}



#######################################################
## Figure 2e: R^2 vs reliability for empirical mrna/prot cor= 0.29
#######################################################

pdf("Figs/5e-noise-correction.pdf")
par(xpd=FALSE)
par(mar=c(6.1,5.1,4.1,6.1))
cols <- c("#08306B","#08519C","#2171B5","#4292C6","#9ECAE1")
cols2 <- rep(cols, each=2)
cols <- sapply(1:9, function(i) avg.colors(cols2[i], cols2[i+1]))
cols <- c(cols, "#D8EAF3")

plot(0, xlim=c(0, 1), ylim=c(0, 1), cex=0,
     xlab="Reliability of mRNA Measurements",
     ylab="Reliability of Protein Measurements",
     main=paste0("Fraction of Across-Tissue Protein Variance\n ",
                 "Explained By Transcript Levels\n"),
     cex.lab=2, cex.axis=1.5, cex.main=1.5, xaxs="i", yaxs="i")

out.prev <- curve(p2m1cor^2/x, from=p2m1cor^2, to=1, col="red", lwd=3, n=1000, add=TRUE)
polygon(c(0, 0, out.prev$x,1), c(0, 1, out.prev$y,0), density = c(10, 20), angle = c(45, -45))
count <- 1
for(cur in seq(0.9,0.1,by=-0.1)){
    col <- cols[count]
    out.cur <- curve(p2m1cor^2/(cur*x), from=p2m1cor^2, to=1,
                     lwd=1, lty=2,col="black", n=1000, add=TRUE)
    polygon(c(out.prev$x, rev(out.cur$x)),
            c(out.prev$y, rev(out.cur$y)), col=col, lty=0)
    out.prev <- out.cur
    count <- count+1
}
polygon(c(out.prev$x, rev(out.prev$x)),
        c(out.prev$y, rep(1, length(out.prev$y))),
        col="#D8EAF3", lty=0) 
null.line <- curve(p2m1cor^2/x, from=p2m1cor^2, to=1,
                   col="red",lwd=3,n=1000,add=TRUE)
text(mrnaRel, protRel, "X", col="red", cex=2)
text(median(mrna.reliabilities), median(protein.reliabilities), "X", col="red", cex=2)

## colorbar
par(xpd=TRUE)
ybot <- 0.15
xleft <- 1.05
rectHeight <- 0.075
rectWidth <- 0.1

for(i in 1:length(cols)) {
    rect(xleft, ybot + rectHeight*(i-1),
         xleft + rectWidth,  ybot + rectHeight*i, col=rev(cols)[i], lwd=0)
    text(xleft + 3/2*rectWidth, ybot + rectHeight*(i-1/2), paste0(i*10, "%"))
}
par(xpd=FALSE)

dev.off() 

#######################################################
## Look at variation in expression and abundance across tissues
#######################################################

## noise corrected 
sd.mrnas <- apply(mrna, 1, function(x) sd(x,na.rm=TRUE))
sd.prots <- apply(protein, 1, function(x) sd(x,na.rm=TRUE))
## temporary
sd.mrnas.corrected <- sqrt(sd.mrnas^2 * mrnaRel)
sd.prots.corrected <- sqrt(sd.prots^2 * p1.reliability)
sd.diffs.corrected <- sd.prots.corrected - sd.mrnas.corrected

pdf("Figs/physiological_sds.pdf")
ggplot(data.frame(x=sd.mrnas.corrected,x2=sd.prots.corrected),labels=c("mRNA","protein"))+geom_histogram(aes(x=x,fill="mRNA"),binwidth=0.03,colour="black",size=0.1,alpha=1)+geom_histogram(aes(x=x2,fill="protein"),binwidth=0.03,colour="black",size=0.1,alpha=0.5)+scale_x_continuous(limits=c(0,1.5))+labs(x="Standard Deviation",y="Count")+
    theme(axis.title.x = element_text(face="bold",size=25),
          axis.text.x  = element_text(size=16),
          axis.title.y = element_text(face="bold",size=25),
          axis.text.y  = element_text(size=16,hjust=0),
          legend.text=element_text(size = 20, face = "bold"),
          legend.position=c(0.8,0.9))+
         guides(fill=guide_legend(title=NULL))
dev.off()

pdf("Figs/sds-diff.pdf")
ggplot(data.frame(x=(sd.prots.corrected-sd.mrnas.corrected)))+geom_histogram(aes(x=x,fill="protein-mRNA"),binwidth=0.05,colour="black",size=0.1,alpha=1)+
    scale_x_continuous(limits=c(-1.5,1.5))+
labs(x="Difference in Standard Deviation",y="Count")+
    theme(axis.title.x = element_text(face="bold",size=25),
          axis.text.x  = element_text(size=16),
          axis.title.y = element_text(face="bold",size=25),
          axis.text.y  = element_text(size=16,hjust=0),
          legend.text=element_text(size = 20, face = "bold"),
          legend.position=c(0.8,0.9))+
     guides(fill=guide_legend(title=NULL))+geom_vline(xintercept=0,size=1.5)
dev.off()

median(sd.prots.corrected/sd.mrnas.corrected,na.rm=TRUE)


################################################
### Fold Change Figs in supplement: Choose 3 tissues
################################################

fold.comparison.tissues <- c("prostate","kidney","liver")
fold.combos <- combn(fold.comparison.tissues, 2)

for(i in 1:ncol(fold.combos)){

    t1 <- fold.combos[1,i]
    t2 <- fold.combos[2,i]

    pdf(sprintf("Figs/2_top_%s_%s_fold_cor.pdf",t1,t2))

    fold.cor <- cor(mrna[, t1] - mrna[, t2],
                    protein[, t1] - protein[, t2], use="complete.obs")

    print(ggplot(data=data.frame(mrna=10^(mrna[,t1]-mrna[,t2]),
                                 protein=10^(protein[,t1]-protein[,t2]))) +
          geom_point(aes(x=mrna,y=protein),colour="blue",alpha=0.3) +
          scale_x_log10(breaks=10^(-2:2), label=scientific_10,
                        limits=c(10^-2,10^2)) +
          scale_y_log10(breaks=10^(-3:3), labels=scientific_10,
                        limits=c(10^-3,10^3)) +
          labs(x="mRNA Fold Change",y="Protein Fold Change")+
          annotate("text", x = 10^-1.2, y = 10^2.5,label=as.character(as.expression(substitute(italic(R)~"="~r,list(r=round(fold.cor,digits=3))))),parse=TRUE,size=10)+
          annotate("text", x = 10^-1.2, y = 10^2,label=as.character(as.expression(substitute(italic(R)^2~"="~r2,list(r2=round(fold.cor^2,digits=3))))),parse=TRUE,size=10)+
          theme(axis.title.x = element_text(face="bold",size=30),
                axis.text.x  = element_text(size=30),
                axis.title.y = element_text(face="bold",size=30),
                axis.text.y  = element_text(size=30,hjust=0)))
    
    dev.off()

    pdf(sprintf("Figs/2_%s_scaledby_%s.pdf", t1, t2))
    
    pred <- mrna[, t1] + (protein[, t2] - mrna[, t2])
    raw <- protein[, t1]
    cor.raw <- cor(pred, raw, use="complete.obs")
    print(
        ggplot(data=data.frame(mrna=pred, protein=raw)) +
        geom_point(aes(x=mrna, y=protein), colour="blue",alpha=0.3) +
        theme(plot.margin=unit(c(5.5, 12.5, 5.5, 5.5), "points")) +
        scale_x_continuous(breaks=1:10,labels=function(x) parse(text=paste("10^", x)),limits=c(3,10))+
        scale_y_continuous(breaks=1:10,labels=function(x) parse(text=paste("10^", x)),limits=c(3,10))+
        labs(x=sprintf("%s mRNA scaled by %s PTR", strFormat(t1), strFormat(t2)),
             y=sprintf("Measured Protein in %s",strFormat(t1)))+
        annotate("text", x = 4, y = 9.5,label=as.character(as.expression(substitute(italic(R)~"="~r,list(r=format(cor.raw, digits=2))))),parse=TRUE,size=10)+
        annotate("text", x = 4, y = 9,label=as.character(as.expression(substitute(italic(R)^2~"="~r2,list(r2=format(cor.raw^2, digits=2))))),parse=TRUE,size=10)+
        theme(axis.title.x = element_text(face="bold",size=22),
              axis.text.x  = element_text(size=22),
              axis.title.y = element_text(face="bold",size=30),
              axis.text.y  = element_text(size=30,hjust=0))
        )
    dev.off()
    
    
    
}

########################################
## Raw Correlations 
########################################

for(tis in fold.comparison.tissues){

    pdf(sprintf("Figs/2_%s_raw_cor.pdf", tis))

    pred <- log10(predicted.raw[, tis])
    raw <- protein[,tis]
    keep <- which(pred!=raw)
    pred <- pred[keep]
    raw <- raw[keep]
    cor.raw <- cor(pred,raw,use="complete.obs")^2
    print(
        ggplot(data=data.frame(mrna=pred,protein=raw))+geom_point(aes(x=mrna,y=protein),colour="blue",alpha=0.3) +
        theme(plot.margin=unit(c(5.5, 12.5, 5.5, 5.5), "points")) +
        scale_x_continuous(breaks=1:10,labels=function(x) parse(text=paste("10^", x)),limits=c(3,10))+
        scale_y_continuous(breaks=1:10,labels=function(x) parse(text=paste("10^", x)),limits=c(3,10))+
        labs(x=sprintf("Scaled mRNA in %s",strFormat(tis)),
             y=sprintf("Measured Protein in %s",strFormat(tis)))+
        annotate("text", x = 4, y = 9.5,label=as.character(as.expression(substitute(italic(R)~"="~r,list(r=format(cor.raw,digits=2))))),parse=TRUE,size=10)+
        annotate("text", x = 4, y = 9,label=as.character(as.expression(substitute(italic(R)^2~"="~r2,list(r2=format(cor.raw^2,digits=2))))),parse=TRUE,size=10)+
        theme(axis.title.x = element_text(face="bold",size=30),
              axis.text.x  = element_text(size=30),
              axis.title.y = element_text(face="bold",size=30),
              axis.text.y  = element_text(size=30,hjust=0))
        )
    dev.off()
    
}


#############################################################################
############ Raw Correlations using single tissue PTR (not median) ##########
#############################################################################

for(i in 1:ncol(fold.combos)){

    t1 <- fold.combos[1,i]
    t2 <- fold.combos[2,i]

    pdf(sprintf("Figs/2_bottom_%s_from_%s.pdf",t1,t2))
    pred <- log10(tissue.ratios[,t2]*10^mrna[,t1])
    raw <- protein[,t1]
    keep <- which(pred!=raw)
    pred <- pred[keep]
    raw <- raw[keep]
    cor.raw <- cor(pred,raw,use="complete.obs")^2
    print(
        ggplot(data=data.frame(mrna=pred,protein=raw))+geom_point(aes(x=mrna,y=protein),colour="blue",alpha=0.3)+
        scale_x_continuous(breaks=1:10,labels=function(x) parse(text=paste("10^", x)),limits=c(2,8))+
        scale_y_continuous(breaks=1:10,labels=function(x) parse(text=paste("10^", x)),limits=c(2,8))+
        labs(x=sprintf("%s mRNA scaled by %s PTR",strFormat(t1),strFormat(t2)),
             y=sprintf("Measured protein in %s",strFormat(t1)))+
        annotate("text", x = 3, y = 7.5,label=as.character(as.expression(substitute(italic(R)~"="~r,list(r=format(cor.raw,digits=2))))),parse=TRUE,size=10)+
        annotate("text", x = 3, y = 7,label=as.character(as.expression(substitute(italic(R)^2~"="~r2,list(r2=format(cor.raw^2,digits=2))))),parse=TRUE,size=10)+
        theme(axis.title.x = element_text(face="bold",size=22),
              axis.text.x  = element_text(size=30),
              axis.title.y = element_text(face="bold",size=30),
              axis.text.y  = element_text(size=30,hjust=0))
        )
    dev.off()
}
























## look in associations.R for group2proteins and allGroups
fdr.thresh <- 0.01

atable <- find_go_groups(z.score, fdr.thresh, within.cors, sd.mrnas, sd.prots)
med.cor <- atable[,"Median Correlation"]
xtab <- xtable(atable[order(atable[,3]),c(1:3,6:8)],display=c("s","s","f","f","g","g","d"),digits=c(0,0,2,2,2,2,2))

print(xtab,type='latex',floating.environment = 'sidewaystable',hline.after=c(0,nrow(atable)-sum(med.cor>pop.cor)))

cors.tmp <- cat.tmp <- c()

## show just a sample of the groups
toinclude <- rownames(atable)[c(sort(sample(which(med.cor<pop.cor),3)),which(med.cor>pop.cor))]

for(grp in toinclude) {
     prots <- intersect(n2e(group2proteins[[grp]]),names(z.score))
     cors.tmp <- c(cors.tmp, within.cors[prots])
     cat.tmp <- c(cat.tmp, rep(as.character(atable[grp, "Name"]), length(prots)))
}

df <- data.frame(list(cors=cors.tmp,cat=factor(cat.tmp,levels=unique(cat.tmp),ordered=TRUE)))

pdf("Figs/go_corrs.pdf", width=5, height=5)
print(
    ggplot(df, aes(x = cat, y = cors)) + 
    geom_boxplot(aes(x=cat, y=cors), outlier.size=0) + 
    scale_y_continuous(limits=c(-.5,0.8), 
                       breaks=c(-.5,0, round(pop.cor,2),.5,0.8)) + 
    labs(x="", y="Correlation") + 
    theme(axis.text.x = element_text(angle = 45,hjust=1),
          plot.margin=unit(c(1, 2, 1/2, 1), "cm")) + 
    geom_hline(yintercept=round(pop.cor,2),colour="red",lwd=1,lty=2)
)
dev.off()


prel.info <- getZscores(protein[prows, pcols], protein2[prows, pcols])
prel.groups <- find_go_groups(prel.info$z.score, fdr.thresh, prel.info$within.cors, prel.info$sd.mat1, prel.info$sd.mat2)

med.cor <- prel.groups[,"Median Correlation"]
xtab <- xtable(prel.groups[order(prel.groups[,3]),c(1:3,6:8)],display=c("s","s","f","f","g","g","d"),digits=c(0,0,2,2,2,2,2))

print(xtab, type='latex', floating.environment = 'sidewaystable', 
      hline.after=c(0, nrow(prel.groups) - sum(med.cor > prel.info$pop.cor)))

pdf("Figs/empirical_protein_rels.pdf")
ggplot(data.frame(x=prel.info$within.cors)) + 
    geom_histogram(aes(x=x,fill="within"), colour="black", binwidth=0.05,size=0.1)+labs(x="Correlation", y="",title="Protein Reliabilities") +
    theme(title = element_text(size=20),
          axis.title.x = element_text(face="bold",size=25),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16,hjust=0)) +
    scale_y_continuous(limits=c(0,300))+
    guides(fill=FALSE)
dev.off()



mrna_split1<- as.matrix(read.csv(sprintf("data/mrna_%s_split1.csv", mrnanm1) ,row.names=1))
mrna_split2<- as.matrix(read.csv(sprintf("data/mrna_%s_split2.csv", mrnanm1) ,row.names=1))
mrel.info <- getZscores(mrna_split1, mrna_split2)

mrna2_split1 <- as.matrix(read.csv(sprintf("data/mrna_%s_split1.csv", mrnanm2) ,row.names=1))
mrna2_split2 <- as.matrix(read.csv(sprintf("data/mrna_%s_split2.csv", mrnanm2) ,row.names=1))
mrel.info <- getZscores(mrna2_split1, mrna2_split2)

mrel.info <- getZscores(msplit1, msplit2)
mrel.groups <- find_go_groups(mrel.info$z.score, fdr.thresh,mrel.info$within.cors,mrel.info$sd.mat1, mrel.info$sd.mat2)

med.cor <- mrel.groups[,"Median Correlation"]
xtab <- xtable(mrel.groups[order(mrel.groups[,3]),c(1:3,6:8)],display=c("s","s","f","f","g","g","d"),digits=c(0,0,2,2,2,2,2))

print(xtab, type='latex', floating.environment = 'sidewaystable', 
      hline.after=c(0, nrow(mrel.groups) - sum(med.cor > mrel.info$pop.cor)))

pdf("Figs/empirical_mrna_rels.pdf")
ggplot(data.frame(x=mrel.info$within.cors)) + 
    geom_histogram(aes(x=x,fill="within"), colour="black", binwidth=0.05, size=0.1)+labs(x="Correlation", y="",title="mRNA Reliabilities") +
    theme(title = element_text(size=20),
          axis.title.x = element_text(face="bold",size=25),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16,hjust=0)) +
    scale_y_continuous(limits=c(0,1500))+
    guides(fill=FALSE)
dev.off()




 cors.mat <- miss.mat <- matrix(0,nrow=length(pcols),ncol=length(pcols))
 rownames(cors.mat) <- colnames(cors.mat) <- pcols
 rownames(miss.mat) <- colnames(miss.mat) <- pcols
 combinations <- combn(pcols,2)
 for(i in 1:ncol(combinations)){

     t1 <- combinations[1,i]
     t2 <- combinations[2,i]
    
     foldRelProt <- cor(protein[prows,t1]-protein[prows,t2],protein2[prows,t1]-protein2[prows,t2],use="pairwise.complete.obs")
     foldRelMrna <- cor(msplit1[prows,t1] - msplit1[prows,t2], msplit2[prows,t1] - msplit2[prows,t2],use="pairwise.complete.obs")
     rawCor <- cor(protein[prows,t1] - protein[prows,t2], mrna[prows,t1] - mrna[prows,t2],use="pairwise.complete.obs")
     
     ## obs.frac <- mean(!is.na((protein[prows,t1]-protein[prows,t2])*(protein2[prows,t1]-protein2[prows,t2])))
     
     if( foldRelProt < 0 || foldRelMrna < 0 | rawCor < 0) {
         cors.mat[combinations[1,i],combinations[2,i]] <- cors.mat[combinations[2,i],combinations[1,i]] <- NA
     } else {
         
         cors.mat[combinations[1,i],combinations[2,i]] <- cors.mat[combinations[2,i],combinations[1,i]] <- rawCor / sqrt( foldRelProt * foldRelMrna )
     
     }
     
##     miss.mat[combinations[1,i],combinations[2,i]] <- miss.mat[combinations[2,i],combinations[1,i]] <- obs.frac
    
     
 }
 rownames(cors.mat) <- colnames(cors.mat) <- sapply(rownames(cors.mat),strFormat)

cors.mat[cors.mat > 1] <- 1
cors.mat[cors.mat < 0] <- NA

## pdf("Figs/corr_mat.pdf",width=9,height=7)
col.pal <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F", 
                             "cyan", "#007FFF", "blue","#00007F"))
corrplot(cors.mat, method = "number",cl.lim=c(-1,1),col=col.pal(100),diag=FALSE,tl.col="black",tl.cex=1.5,cl.cex=1.2,cl.align.text="l",mar=c(5,0,0,0))
## Prettify names
## dev.off()


fold.comparison.tissues <- c("prostate","kidney","thyroid.gland")
fold.combos <- combn(fold.comparison.tissues,2)

for(i in 1:ncol(fold.combos)){

    t1 <- fold.combos[1,i]
    t2 <- fold.combos[2,i]

##    pdf(sprintf("Figs/2_top_%s_%s_fold_cor.pdf",t1,t2))

    fold.cor <- cor(mrna[prows,t1]-mrna[prows,t2],protein[prows,t1]-protein[prows,t2],use="complete.obs")
    fold.cor <- cor(protein[prows,t1]-protein[prows,t2],protein2[prows,t1]-protein2[prows,t2],use="complete.obs")

    print(ggplot(data=data.frame(protein2=10^(protein2[prows,t1]-protein2[prows,t2]),protein=10^(protein[prows,t1]-protein[prows,t2])))+geom_point(aes(x=protein2,y=protein),colour="blue",alpha=0.3)+
    scale_x_log10(breaks=10^(-3:3),label=scientific_10,limits=c(10^-3,10^3))+
    scale_y_log10(breaks=10^(-3:3),labels=scientific_10,limits=c(10^-3,10^3))+
    labs(x="Kim Fold Change",y="Wilhelm Fold Change")+
    annotate("text", x = 10^-1.2, y = 10^2.5,label=as.character(as.expression(substitute(italic(R)~"="~r,list(r=round(fold.cor,digits=3))))),parse=TRUE,size=10)+
    annotate("text", x = 10^-1.2, y = 10^2,label=as.character(as.expression(substitute(italic(R)^2~"="~r2,list(r2=round(fold.cor^2,digits=3))))),parse=TRUE,size=10)+
    theme(axis.title.x = element_text(face="bold",size=30),
          axis.text.x  = element_text(size=30),
          axis.title.y = element_text(face="bold",size=30),
          axis.text.y  = element_text(size=30,hjust=0)))
    
##    dev.off()
    
}

