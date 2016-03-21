rm(list=ls())
library(gridExtra)
source("util_functions.R")
load("matrix-Kim-illumina.RData")
kimTM <- transformed.mat
kimCG <- clustered.groups
kimGo <- unrefined.mat

load("matrix-Wilhelm-pa.RData")
wilhelmTM <- transformed.mat
wilhelmCG <- clustered.groups
wilhelmGo <- unrefined.mat

rnms <- intersect(rownames(wilhelmCG), rownames(kimCG))
cnms <- intersect(colnames(wilhelmCG), colnames(kimCG))
wilhelmCG <- wilhelmCG[rnms, cnms]
kimCG <- kimCG[rnms, cnms]

pdf("Figs/rptr-corrs.pdf")
par(mfrow=c(4, 3))
for(i in 1:ncol(wilhelmCG)) {
    tissue <- colnames(wilhelmCG)[i]
    cr <- cor.test(wilhelmCG[, i], kimCG[, i], use="pairwise.complete.obs")

    par(mar=c(2, 2, 2, 2))
    plot(wilhelmCG[, i], kimCG[,  i], pch=19, cex=0.7, col="red",
         main=strFormat(tissue), xlab="", ylab="", xlim=c(-1.2, 1.2), ylim=c(-1.2,1.2), xaxt="n", yaxt="n")
    axis(side=1, at=c(-1, 0, 1))
    axis(side=2, at=c(-1, 0, 1), las=2)
    abline(h=0)
    abline(v=0)
    textcol <- ifelse(cr$p.value < 0.01, "blue", "grey")
    legend("bottomright", sprintf("Cor = %s", format(cr$estimate, digits=2)),
           text.col=textcol, bg="white", box.lty=0, text.font=2)
    box()
}
dev.off()

## Regulation mat terms to keep:



## switch to log2 

rptrMax <- 1
rptrMin <- -1
kimCG2 <- kimCG
kimCG2[kimCG2 > rptrMax] <- rptrMax - 1e-9
kimCG2[kimCG2 < rptrMin] <- rptrMin + 1e-9

wilhelmCG2 <- wilhelmCG
wilhelmCG2[wilhelmCG2 > rptrMax] <- rptrMax - 1e-9
wilhelmCG2[wilhelmCG2 < rptrMin] <- rptrMin + 1e-9

##pdf(sprintf("Figs/regulation-mat-within-%s-%s-%s.pdf", test.type,
##            proteinDataSource, mrnaDataSource), height=15)
scales1 <- list(x=list(at=1:ncol(kimCG), rot=90,
                       lab=sapply(colnames(kimCG2), strFormat)),
               y=list(at=1:nrow(kimCG2),
                      lab=go.names[rownames(kimCG),3]))
scales2 <- list(x=list(at=1:ncol(wilhelmCG), rot=90,
                       lab=sapply(colnames(wilhelmCG), strFormat)),
               y=list(at=1:nrow(wilhelmCG2),
                      lab=rep("", 40)))

cols <- colorRampPalette(c("blue", "lightblue", "white", "orangered", "red"))(100)

mat1 <- image(Matrix(kimCG2[1:40, ]),
            at=seq(rptrMin, rptrMax, length.out=20),
            scales=scales1, xlab="",ylab="",sub="",
            col.regions=cols)


mat2 <- image(Matrix(wilhelmCG2[1:40, ]), main="Kim et. al",
            at=seq(rptrMin, rptrMax, length.out=20),
            scales=scales2, xlab="",ylab="",sub="",
            col.regions=cols)


pdf("Figs/mat1.pdf")
print(mat1)
dev.off()
pdf("Figs/mat2.pdf")
print(mat2)
dev.off()


############################################
nms <- intersect(rownames(wilhelmTM), rownames(kimTM))
cnms <- intersect(colnames(wilhelmTM), colnames(kimTM))

wilhelmTM <- wilhelmTM[rnms, cnms]
kimTM <- kimTM[rnms, cnms]

wilhelmTM[wilhelmTM==Inf] <- NA
kimTM[kimTM==Inf] <- NA

plot(wilhelmTM[, "kidney"], kimTM[, "kidney"], pch=19, cex=0.1)
for(i in 1:ncol(wilhelmTM)) {
    print(colnames(wilhelmTM)[i])
    print(cor.test(wilhelmTM[, i], kimTM[, i], use="pairwise.complete.obs"))
}
