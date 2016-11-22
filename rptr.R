library(gridExtra)
library(xtable)
source("util_functions.R")
load("matrix-Kim-illumina.RData")
kimTM <- transformed.mat
kimCG <- clustered.groups
kimGo <- unrefined.mat

load("matrix-Wilhelm-pa.RData")
wilhelmTM <- transformed.mat
wilhelmCG <- clustered.groups
wilhelmGo <- unrefined.mat

wilhelmCGorig <- wilhelmCG
kimCGorig <- kimCG

###############
## Only significant terms
###############

rnms <- intersect(rownames(wilhelmCG), rownames(kimCG))
cnms <- intersect(colnames(wilhelmCG), colnames(kimCG))
wilhelmCG <- wilhelmCG[rnms, cnms]
kimCG <- kimCG[rnms, cnms]

pdf("Figs/rptr-corrs.pdf")
par(mfrow=c(4, 3), oma=c(0, 1, 0, 0))
for(i in 1:ncol(wilhelmCG)) {
    tissue <- colnames(wilhelmCG)[i]
    cr <- cor.test(wilhelmCG[, i], kimCG[, i], use="pairwise.complete.obs")

    par(mar=c(2, 2, 2, 2))
    plot(wilhelmCG[, i], kimCG[,  i], pch=19, cex=0.7, col="red",
         main=strFormat(tissue), xlab="", ylab="", xlim=c(-1.2, 1.2), ylim=c(-1.2,1.2), xaxt="n", yaxt="n", cex.main=2)
    axis(side=1, at=c(-1, 0, 1), cex.axis=1.5)
    axis(side=2, at=c(-1, 0, 1), las=2, cex.axis=1.5)
    abline(h=0)
    abline(v=0)
    textcol <- ifelse(cr$p.value < 0.01, "blue", "grey")
    legend("bottomright", sprintf("Cor = %s", format(cr$estimate, digits=2)),
           text.col=textcol, box.lty=0, text.font=2, cex=1.4)
    box()
}
dev.off()


###############
## All terms
###############

rnms <- intersect(rownames(wilhelmTM), rownames(kimTM))
cnms <- intersect(colnames(wilhelmTM), colnames(kimTM))
wilhelmTM <- wilhelmTM[rnms, cnms]
kimTM <- kimTM[rnms, cnms]
wilhelmTM[wilhelmTM==Inf] <- NA
kimTM[kimTM==Inf] <- NA

rptrAllTable <- c()
for(i in 1:ncol(wilhelmTM)) {
    tissue <- colnames(wilhelmTM)[i]
    cr <- cor.test(wilhelmTM[, i], kimTM[, i], use="pairwise.complete.obs")

    rptrAllTable <- rbind(rptrAllTable, c(cr$estimate, cr$conf.int[1], cr$conf.int[2]))
}
rownames(rptrAllTable) <- sapply(colnames(wilhelmTM), function(x) strFormat(x))
colnames(rptrAllTable) <- c("Estimate", "Lower", "Upper")
print(xtable(t(rptrAllTable)))

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

wilhelmGoTerms <- c("GO:0006614", "GO:0006414", "GO:0005840", "GO:0006415", "GO:0006418", "GO:0000049", "GO:0006412", "GO:0006413", "GO:0019432", "GO:0022904", "GO:0005753", "GO:0005753", "GO:0042776", "GO:0031966", "GO:0006091","GO:0016651", "GO:0005743", "GO:0051539", "GO:0005739", "GO:0005758", "GO:0009055", "GO:0005759", "GO:0055114", "GO:0006635", "GO:0003995", "GO:0009083", "GO:0004029", "GO:0050660", "GO:0035338", "GO:0016491", "GO:0006633", "GO:0006766", "GO:0006099", "GO:0000096", "GO:0010873", "GO:0042157", "GO:0001523", "GO:0031093", "GO:0046951", "GO:0006695")
wilhelmHandSelected <- wilhelmCG[wilhelmGoTerms, ]
wilhelmHandSelected[wilhelmHandSelected > rptrMax] <- rptrMax - 1e-9
wilhelmHandSelected[wilhelmHandSelected < rptrMin] <- rptrMin + 1e-9

##pdf(sprintf("Figs/regulation-mat-within-%s-%s-%s.pdf", test.type,
##            proteinDataSource, mrnaDataSource), height=15)
scales1 <- list(x=list(at=1:ncol(kimCG2), rot=90,
                       lab=sapply(colnames(kimCG2), strFormat)),
               y=list(at=1:nrow(kimCG2),
                      lab=go.names[rownames(kimCG2),3]), alternating=c(1, 0))
scales2 <- list(x=list(at=1:ncol(wilhelmCG), rot=90,
                       lab=sapply(colnames(wilhelmCG), strFormat)),
               y=list(at=1:nrow(wilhelmCG2),
                      lab=rep("", 40)))

scales3 <- list(x=list(at=1:ncol(wilhelmHandSelected), rot=90,
                       lab=sapply(colnames(wilhelmHandSelected), strFormat)),
               y=list(at=1:nrow(wilhelmHandSelected),
                      lab=go.names[rownames(wilhelmHandSelected), 3]))

cols <- colorRampPalette(c("blue", "dodgerblue", "white", "orangered", "red"))(100)

mat1 <- image(Matrix(kimCG2[1:40, ]),
            at=seq(rptrMin, rptrMax, length.out=20),
            scales=scales1, xlab="",ylab="",sub="",
            col.regions=cols)


mat2 <- image(Matrix(wilhelmCG2[1:40, ]), main="Wilhelm et. al",
            at=seq(rptrMin, rptrMax, length.out=20),
            scales=scales2, xlab="",ylab="",sub="",
            col.regions=cols)

mat3 <- image(Matrix(wilhelmHandSelected), 
            at=seq(rptrMin, rptrMax, length.out=20),
            scales=scales3, xlab="",ylab="",sub="",
            col.regions=cols)


pdf("Figs/mat1.pdf", width=8, height=8)
print(mat1)
dev.off()
pdf("Figs/mat2.pdf")
print(mat2)
dev.off()
pdf("Figs/mat3.pdf")
print(mat3)
dev.off()


rownames(wilhelmCGorig) <- go.names[rownames(wilhelmCGorig),3]
rownames(kimCGorig) <- go.names[rownames(kimCGorig),3]

both <- intersect(rownames(wilhelmGo), rownames(kimGo))
wilhelmBoth <- wilhelmCGorig[both, ]
rownames(wilhelmBoth) <- go.names[both, 3]

write.csv(wilhelmCGorig, file="~/Desktop/wilhelmTable.csv")
write.csv(kimCGorig, file="~/Desktop/kimTable.csv")
write.csv(wilhelmBoth, file="~/Desktop/wilhelmRestricted.csv")

wilhelmCGorig[rnms, ]

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
