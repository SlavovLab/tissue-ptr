library(plyr)
library(mvtnorm)
source("util_functions.R")

mcp <- read.csv("~/Dropbox/NatureComment/ENCODE_mRNA_Data/mcp-mrna.csv", row.names=1)
mcp[mcp==0] <- NA

ilumina <- read.csv("~/course/tissue-ptr/mrna_ilumina.csv", row.names=1)

encode <- read.csv("~/Dropbox/NatureComment/ENCODE_mrna_Data/57epigenomes.exon.RPKM.pc.csv")
encodeLegend <- read.table("~/Dropbox/NatureComment/ENCODE_mrna_Data/57epigenomes.exon.RPKM.pc_EG.name.txt", row.names=1, stringsAsFactors=FALSE)

encodeColnames <-  tolower(sapply(colnames(encode)[-(1:2)],
                          function(nm) encodeLegend[nm, ]))
indices <- sapply(c("adult_liver", "gastric"),
                  function(x) grep(x, encodeColnames))
encodeColnames[indices] <- c("liver", "stomach")

colnames(encode)[-(1:2)] <- encodeColnames
encode[encode==0] <- NA

## average within gene 
avgEncode <- ddply(as.data.frame(encode), "gene_id", function(x) {
    colMeans(x[, -(1:2)], na.rm=TRUE)
})
rownames(avgEncode) <- avgEncode[, 1]
avgEncode <- avgEncode[, -1]

write.csv(avgEncode, "~/Dropbox/NatureComment/ENCODE_mrna_Data/encode_processed.csv")

cols <- intersect(colnames(mcp), colnames(avgEncode))
## rows <- intersect(rownames(mcp), rownames(avgEncode))
cols <- intersect(colnames(mcp), colnames(ilumina))
rows <- intersect(rownames(mcp), rownames(ilumina))

## out <- getZscores(as.matrix(mcp[rows, cols]), as.matrix(avgEncode[rows, cols]))
out <- getZscores(as.matrix(mcp[rows, cols]), as.matrix(ilumina[rows, cols]))
within.cors <- out$within.cors

sim.within.cors <- sapply(1:length(out$z.score), function(i) {
    cor(rmvnorm(n=out$n.pairwise[i], 
                mean=c(0,0), 
                sigma=matrix(c(1, out$pop.cor, out$pop.cor, 1), nrow=2)))[1,2] 
})

pdf("Figs/mrna_sim_and_within_cors.pdf")
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
    scale_y_continuous(limits=c(0, 10000)) +
    scale_x_continuous(limits=c(-1, 1)) +
    geom_vline(xintercept = out$pop.cor, linetype=2)
dev.off()


cor(avgEncode[, c("aorta", "left_ventricle", "right_ventricle", "right_atrium")], use="pairwise.complete.obs")
