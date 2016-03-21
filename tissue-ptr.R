## ----results='hide',echo=FALSE,warning=FALSE,message=FALSE,label=initialize,cache=FALSE----

options(stringsAsFactors=FALSE)
opts_chunk$set(cache.path = 'cache/',fig.path = 'Figs/')

library(ggplot2)
library(scales)
library(corrplot)
library(RColorBrewer)
library(mvtnorm)
library(fdrtool)
library(xtable)
library(grid)
source("util_functions.R")

## data for the primary analysis comes from Wilhelm
protnm1 <- "Wilhelm"
mrnanm1 <- "illumina"

protnm2 <- "Kim"
mrnanm2 <- "pa"

mrna <- as.matrix(read.csv(sprintf("data/mrna_%s.csv", mrnanm1) ,row.names=1))
mrna2 <- as.matrix(read.csv(sprintf("data/mrna_%s.csv", mrnanm2) ,row.names=1))
protein <- as.matrix(read.csv(sprintf("data/protein_%s.csv", protnm1),row.names=1))
protein2 <- as.matrix(read.csv(sprintf("data/protein_%s.csv", protnm2),row.names=1))

tissue.names <- intersect(colnames(mrna), colnames(protein))
pretty.tissue.names <- sapply(tissue.names, strFormat)
n.tissues <- length(tissue.names)

gene.names <- intersect(rownames(mrna), rownames(protein))
n.genes <- length(gene.names)

mrna <- mrna[gene.names, tissue.names]
protein <- protein[gene.names, tissue.names]

tissue.ratios <- (10^protein) / (10^mrna)
median.ratios <- apply(tissue.ratios, 1, function(x) median(x,na.rm=TRUE))
prediction.ratios <- log10(tissue.ratios * median.ratios)
predicted.raw <- (10^mrna) * median.ratios

## ----results='hide', fig.show='hide',label=figure-1,dependson=c("initialize"),echo=FALSE,warning=FALSE----

keep.indices <- which(log10(predicted.raw) != protein)
cor.raw <- cor(log10(predicted.raw[keep.indices]), as.numeric(protein[keep.indices]))
    
png("Figs/1a_mrna_v_protein.png")
print(

    ggplot(data=data.frame(x=predicted.raw[keep.indices], y=10^as.numeric(protein[keep.indices])),aes(x=x,y=y))+
    geom_point(colour="blue",alpha=0.03)+
    scale_x_log10(breaks=10^(1:10),label=scientific_10,limits=c(10^2,10^10))+
    scale_y_log10(breaks=10^(1:10),labels=scientific_10,limits=c(10^2,10^10))+
    labs(x="Scaled mRNA",y="Measured Protein",title="")+
    annotate("text", x = 10^3, y = 10^7,label=as.character(as.expression(substitute(italic(R)[T]^2~"="~r2,list(r2=round(cor.raw^2,digits=2))))),parse=TRUE,size=10)+
    annotate("text", x = 10^3, y = 10^7.6,label=as.character(as.expression(substitute(italic(R)[T]~"="~r,list(r=round(cor.raw,digits=2))))),parse=TRUE,size=10)+
    theme(axis.title.x = element_text(face="bold",size=30),
          axis.text.x  = element_text(size=30),
          axis.title.y = element_text(face="bold",size=30),
          axis.text.y  = element_text(size=30,hjust=0)),
    )
dev.off()


rsquared.vec <- c()
negative.slopes <- c()
negative.slopes.vals <- c()
for(i in 1:nrow(protein)){
    if(sum(!is.na(predicted.raw[i,])&!is.na(protein[i,]))>1){
        lm.fit <- lm(protein[i,]~log10(predicted.raw)[i,])
        slope <- lm.fit$coefficients[2]
        rsq <- summary(lm.fit)$r.squared
        if(rsq<1 & sum(!is.na(predicted.raw[i,]*protein[i,]))>4)
            rsquared.vec <- c(rsquared.vec,rsq)
        if(slope<0){
            negative.slopes <- c(negative.slopes,i)
            negative.slopes.vals <- c(negative.slopes.vals,slope)
        }
    }
}



full.protein <- names(which(apply(protein[negative.slopes,],1,function(x) sum(!is.na(x)))==n.tissues&apply(predicted.raw[negative.slopes,],1,function(x) sum(!is.na(x)))==n.tissues))

sort(rowMeans(protein[full.protein,]))
focus.prot <- "ENSG00000111640"


mrna.vec <- as.numeric(t(log10(predicted.raw)))
protein.vec <- as.numeric(t(protein))
id <- rep(rownames(mrna),each=n.tissues)
df <- data.frame(x=mrna.vec,y=protein.vec,id=id)

mean.prot <- rowMeans(protein[full.protein,])
large.prots <- sample(names(mean.prot[mean.prot>6.5]), 10)
sub.proteins <- unique(c(large.prots,sample(unique(id)[negative.slopes],100)))
       
df.sub <- df[df$id %in% sub.proteins,]
cor.segs <- cor(df.sub$x,df.sub$y,use="complete.obs")

pdf("Figs/1b_mrna_v_protein_segs_new.pdf",width=10,height=10)
print(
ggplot(data=df.sub,aes(x=x,y=y))+geom_smooth(method=lm,se=FALSE,aes(group=id),colour=alpha("black",0.5))+scale_x_continuous(breaks=1:10,label=function(x) parse(text=paste("10^", x)),limits=c(2,9))+
    scale_y_continuous(breaks=1:10,labels=function(x) parse(text=paste("10^", x)),limits=c(2,9))+
    labs(x="Scaled mRNA",y="Measured Protein")+
    annotate("text", x = 3, y = 8,label=as.character(as.expression(substitute(italic(R)[T]^2~"="~r2,list(r2=round(cor.segs^2,digits=2))))),parse=TRUE,size=10)+
    annotate("text", x = 3, y = 8.5,label=as.character(as.expression(substitute(italic(R)[T]~"="~r,list(r=round(cor.segs,digits=2))))),parse=TRUE,size=10)+
    theme(axis.title.x = element_text(face="bold",size=30),
          axis.text.x  = element_text(size=30),
          axis.title.y = element_text(face="bold",size=30),
          axis.text.y  = element_text(size=30,hjust=0))+
    geom_point(data=data.frame(x=log10(predicted.raw[focus.prot,]),y=protein[focus.prot,],Tissue=tissue.names),aes(x=x,y=y,label=length(tissue.names):1,color=Tissue),size=3)+#scale_colour_discrete(guide=FALSE)+
    theme(legend.text=element_text(size=15),legend.title.align=0.5,
          legend.background=element_rect(fill=alpha("gray", 0.5)),legend.key=element_blank(),legend.title=element_text(size=14)) +
    guides(colour = guide_legend(override.aes = list(size=5)))+
    geom_smooth(data=data.frame(x=log10(predicted.raw[focus.prot,]),y=protein[focus.prot,]),method=lm,se=FALSE,colour="black",size=1.5)+theme(legend.justification=c(1,0), legend.position=c(1,0)),
)
dev.off()

within.indices <- which(apply(!is.na(protein * predicted.raw), 1, sum) > 3)
within.cors <- sapply(within.indices,function(i) cor(protein[i,],log(predicted.raw)[i,],use="pairwise.complete.obs"))
between.cors <- sapply(1:ncol(protein),function(i) cor(protein[,i],log10(predicted.raw)[,i],use="pairwise.complete.obs"))
between.cors.raw <- sapply(1:ncol(protein),function(i) cor(protein[,i],mrna[,i],use="pairwise.complete.obs"))


## 
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

### Histogram of R-squared
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

pdf("Figs/5e-noise-correction.pdf")
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

corrected.cor <- pop.cor/(Rprot*Rmrna)
corrected.rsq <- corrected.cor^2

Rprot.grid <- seq(pop.cor^2,1,by=0.01)
Rmrna.grid <- seq(pop.cor^2,1,by=0.01)
Rgrid <- expand.grid(Rprot.grid,Rmrna.grid)
rsq.grid <- matrix(pop.cor/sqrt(Rgrid[,1]*Rgrid[,2]),nrow=length(Rprot.grid))^2

rsq.grid[rsq.grid>1] <- 1
pdf("Figs/5e-noise-correction.pdf")
par(mar=c(6.1,5.1,4.1,2.1))
cols <- c("#08306B","#08519C","#2171B5","#4292C6","#9ECAE1")
cols2 <- rep(cols,each=2)
cols <- sapply(1:9,function(i) avg.colors(cols2[i],cols2[i+1]))

plot(0,xlim=c(pop.cor^2,1),ylim=c(pop.cor^2,1),xlab="Reliability of mRNA Measurements",ylab="Reliability of Protein Measurements",main="Fraction of Across-Tissue Protein Variance\n Explained By Transcript Levels",cex.lab=2,cex.axis=1.5,cex.main=1.5,xaxs="i",yaxs="i")
out.prev <- curve(pop.cor^2/x,from=0.1,to=1,col="red",lwd=3,n=1000,add=TRUE)
polygon(c(0,out.prev$x,1),c(0,out.prev$y,0),density = c(10, 20), angle = c(45, -45))
count <- 1
for(cur in seq(0.9,0.1,by=-0.1)){
    col <- cols[count]
    out.cur <- curve(pop.cor^2/(cur*x),from=0.1,to=1,lwd=1,lty=2,col="black",n=1000,add=TRUE)
    polygon(c(out.prev$x,rev(out.cur$x)),c(out.prev$y,rev(out.cur$y)),col=col,lty=0)
    out.prev <- out.cur
    count <- count+1
}
null.line <- curve(pop.cor^2/x,from=0.1,to=1,col="red",lwd=3,n=1000,add=TRUE)
text(0.9,0.9,"10",cex=1.5,col="dark grey",srt=-45)
text(0.75,0.72,"20",cex=1.5,col="dark grey",srt=-45)
text(0.65,0.6,"30",cex=1.5,col="dark grey",srt=-45)
text(0.58,0.53,"40",cex=1.5,col="dark grey",srt=-45)
text(0.52,0.47,"50",cex=1.5,col="dark grey",srt=-45)
text(0.48,0.43,"60",cex=1.2,col="dark grey",srt=-45)
text(0.45,0.40,"70",cex=1.2,col="dark grey",srt=-45)
text(0.43,0.37,"80",cex=1,col="dark grey",srt=-45)
text(0.41,0.35,"90",cex=1,col="dark grey",srt=-45)
rect(0.18,0.22,0.35,0.32,col="white",lty=1)
text(0.27,0.27,"100%",cex=2,col="dark grey")
dev.off() 

## Compute mrna reliabilities
mrnaRows <- intersect(rownames(mrna), rownames(mrna2))
mrnaCols <- intersect(colnames(mrna), colnames(mrna2))

mrnaRelInfo <- getZscores(mrna[mrnaRows, mrnaCols], mrna2[mrnaRows, mrnaCols])
print(mrnaRelInfo$pop.cor)
mrna.reliabilities <- mrnaRelInfo$within.cors

pdf("Figs/mrna_reliabilities.pdf")
ggplot(data.frame(x=mrna.reliabilities))+geom_histogram(aes(x=x,fill="mRNA"),colour="black",binwidth=0.05,size=0.1)+labs(x="Reliability",y="",title=expression(paste("mRNA Reliabilities, Median = ",med.mrna,sep="")))+
    theme(title = element_text(size=20),
          axis.title.x = element_text(face="bold",size=25),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16,hjust=0)) +
    scale_x_continuous(limits=c(-1,1))+
    guides(fill=FALSE)
dev.off()

## Compute protein reliabilities
pcols <- intersect(colnames(protein), colnames(protein2))
prows <- intersect(rownames(protein), rownames(protein2))
cts <- apply(protein[prows, pcols] * protein2[prows, pcols], 1,
             function(x) {sum(!is.na(x))})

prows <- prows[cts > 4]

ppScores <- getZscores(protein[prows, pcols], protein2[prows, pcols])
med.prot <- ppScores$pop.cor
protein.reliabilities <- ppScores$within.cors

pdf("Figs/protein_reliabilities.pdf")
ggplot(data.frame(x=protein.reliabilities))+geom_histogram(aes(x=x,fill="mRNA"),colour="black",binwidth=0.05,size=0.1)+labs(x="Reliability",y="",title=expression(paste("Protein Reliabilities, Median = ",med.prot,sep="")))+
    theme(title = element_text(size=20),
          axis.title.x = element_text(face="bold",size=25),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16,hjust=0)) +
    scale_x_continuous(limits=c(-1,1))+
    guides(fill=FALSE)
dev.off()

### Compare Kim and Wilhelm

popcor1 <- getZscores(protein[, pcols], mrna[ ,pcols])$pop.cor
popcor2 <- getZscores(protein2[prows, pcols], mrna[prows, pcols])$pop.cor


## Estimated reliabilities
p1.reliability <- med.prot*popcor1/popcor2
p2.reliability <- med.prot*popcor2/popcor1

## Consensus Data

## ratio of error variances
eps1 <- (1-p1.reliability)/p1.reliability
eps2 <- (1-p2.reliability)/p2.reliability

w1 <- 1/eps1/(1/eps1+1/eps2)
w2 <- 1-w1

proteinConsensus <- w1*protein[prows, pcols]+w2*protein2[prows, pcols]
proteinConsensus[which(is.na(protein[prows, pcols]))] <- (protein2[prows, pcols])[which(is.na(protein[prows, pcols]))]
proteinConsensus[which(is.na(protein2[prows, pcols]))] <- (protein[prows, pcols])[which(is.na(protein2[prows, pcols]))]

pc2 <- cbind(proteinConsensus,
             protein[prows, setdiff(colnames(protein), colnames(protein2))],
             protein2[prows, setdiff(colnames(protein2), colnames(protein))])
getZscores(pc2[, colnames(protein)], mrna[prows, ])$pop.cor

getZscores(protein[prows, ], mrna[prows, ])$pop.cor
getZscores(protein[intersect(prows, rownames(mrna2)), ],
           mrna2[intersect(prows, rownames(mrna2)), colnames(protein)])$pop.cor


write.csv(pc2, file="data/protein_consensus.csv")

## ----fig.show='hide',label=figure-s1, dependson("initialize"),echo=FALSE,warning=FALSE, results='hide'----

## noise corrected 
sd.mrnas <- apply(mrna, 1, function(x) sd(x,na.rm=TRUE))
sd.prots <- apply(protein, 1, function(x) sd(x,na.rm=TRUE))
## temporary
sd.mrnas.corrected <- sqrt(sd.mrnas^2*med.mrna)
sd.prots.corrected <- sqrt(sd.prots^2*p1.reliability)
sd.diffs.corrected <- sd.prots.corrected-sd.mrnas.corrected

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

## ----results='asis',fig.show='hide',label=associations-1, dependson("initialize"),echo=FALSE,warning=FALSE----

## look in associations.R for group2proteins and allGroups
fdr.thresh <- 0.02

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


## ----results='asis',fig.show='hide',label=associations-2, dependson("initialize"),echo=FALSE,warning=FALSE----
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


## second protein dataset with mrna
mrnaProtInfo2 <- getZscores(mrna[prows, pcols], protein2[prows, pcols])
mrel.groups <- find_go_groups(mrel.info$z.score,fdr.thresh,mrel.info$within.cors,mrel.info$sd.mat1,mrel.info$sd.mat2)

    withinGroupCors <- matrix(NA,nrow=length(allGroups),ncol=2)
    rownames(withinGroupCors) <- allGroups
    colnames(withinGroupCors) <- c("Protein","mRNA")
    for( group in allGroups ) {
        prots <- n2e(group2proteins[[group]])
        prots <- intersect(rownames(protein),prots)
        pts <- protein[prots,]
        mrs <- mrna[prots,]
        if(length(prots)>3 & length(prots)<100) {
            rws <- which(apply(pts,1,function(x) sum(!is.na(x)))>6 & 
            apply(mrs,1,function(x) sum(!is.na(x)))>6)
            if(length(rws)>3){
                med.pts <- median(setdiff(cor(t(pts[rws,]),use="pairwise.complete.obs"),c(-1,1,NA)))
                med.mrs <- median(setdiff(cor(t(mrs[rws,]),use="pairwise.complete.obs"),c(-1,1,NA)))
                withinGroupCors[group,] <- c(med.pts,med.mrs)
        }
        }
    }
pdf("Figs/withing_group_corrs.pdf",width=5,height=5)
plot(withinGroupCors[,2],withinGroupCors[,1],pch=19,cex=0.4,xlab="mRNA",ylab="Protein",ylim=c(-1,1),xlim=c(-1,1))
dev.off()

## ------------------------------------------------------------------------

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

#######

Reliability.tissues <- pcols
direc2 <- "~/Dropbox/NatureComment/Kim_Search_Results_Data"
peptide.lists <- list()
for(rtissue in Reliability.tissues){

    print(rtissue)
    
    tab <- read.table(sprintf("%s/%s_evidence.txt", direc2, rtissue),
                          sep="\t",stringsAsFactors=FALSE, header=TRUE)
    ## colnames(tab) <- c("Protein","Peptide","Modification","a","b","Intensity")
    tab <- tab[,c("Leading.Razor.Protein", "Sequence", "Modifications", "Intensity")]
    colnames(tab) <- c("Protein", "Peptide", "Modification", "Intensity")
    tab <- subset(tab, Modification=="Unmodified")

    tab <- ddply(tab,.(Peptide), function(x) c(Protein=x$Protein[1],Intensity=sum(as.numeric(x$Intensity))))

    ## Order by protein name
    tab <- tab[order(tab[,"Protein"]),]
    colnames(tab)[3] <- paste("Intensity",rtissue,sep=".")
    tab[,3] <- as.numeric(tab[,3])
    peptide.lists[[rtissue]] <- tab
    
}

merged.df <- Reduce(function(x,y) merge(x,y,by=c("Protein","Peptide"),all.x=TRUE,all.y=TRUE),peptide.lists)
df <- cbind(merged.df,count=apply(merged.df[3:ncol(merged.df)],1,function(x) sum(!is.na(x))))
## df <- df[order(df$count,decreasing=TRUE),]
write.csv(df,file="peptide_intensities_all.csv",quote=FALSE)

for(i in 4:(ncol(df)-1)) {
    df[,i] <- df[,i]*exp(median(log(df[,3])-log(df[,i]),na.rm=TRUE))
}
peptide.cors <- ddply(subset(df,count>=4),~Protein,function(x){
    
    nr <- nrow(x)
    if(nr==1)
        return(NA)

    samp <- sample(1:nr)
    first.indices <- samp[1:floor(nr/2)]
    second.indices <- samp[(floor(nr/2)+1):nr]
    first <- apply(x[first.indices,3:(ncol(x)-1)],2,function(x) mean(x,na.rm=TRUE))
    second <- apply(x[second.indices,3:(ncol(x)-1)],2,function(x) mean(x,na.rm=TRUE))
    if(sum(!is.na(first)*!is.na(second))>3){
        med.cor <- cor(log(first),log(second), use="pairwise.complete.obs")
        return(med.cor)
    } else{
        return(NA)
    }
})

median.reliability <- median(peptide.cors[,2],na.rm=TRUE)^2
