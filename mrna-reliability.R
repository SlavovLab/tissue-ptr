library(plyr)
library(mvtnorm)
source("util_functions.R")

mrna_pa <- read.csv("data/mrna_pa.csv", row.names=1)
mrna_illumina <- read.csv("data/mrna_illumina.csv", row.names=1)
mrna_mcp <- read.csv("data/mrna_mcp.csv", row.names=1)

cols <- intersect(colnames(mrna_pa), colnames(mrna_illumina)) %>% 
    intersect(colnames(mrna_mcp))

rows <- intersect(rownames(mrna_pa), rownames(mrna_illumina)) %>% 
    intersect(rownames(mrna_mcp))

###################
## MCP vs Illumina
###################

out <- getZscores(as.matrix(mrna_mcp[rows, cols]),
                  as.matrix(mrna_illumina[rows, cols]))
within.cors <- out$within.cors

###################
## MCP vs PA
###################

out <- getZscores(as.matrix(mrna_mcp[rows, cols]),
                  as.matrix(mrna_pa[rows, cols]))
within.cors <- out$within.cors

###################
## Illumina vs PA
###################

out <- getZscores(as.matrix(mrna_illumina[rows, cols]),
                  as.matrix(mrna_pa[rows, cols]))
within.cors <- out$within.cors

sim.within.cors <- sapply(1:length(out$z.score), function(i) {
    cor(rmvnorm(n=out$n.pairwise[i], 
                mean=c(0,0), 
                sigma=matrix(c(1, out$pop.cor, out$pop.cor, 1), nrow=2)))[1,2] 
})

###################
## Histogram of correlations for Illumina/PA reliablity
###################

pdf("Figs/mrna_sim_and_within_cors.pdf")
pooled.cors <- c(sim.within.cors, within.cors)
label <- factor(rep(c("Simulated", "Empirical"), 
                    levels=c("Simulated", "Empirical"), 
                    each=length(within.cors)), ordered=TRUE)
ggplot(data.frame(cors=pooled.cors, cut=label), aes(x=cors, fill=label)) + 
geom_histogram(position=position_dodge(width=0.05), binwidth=0.08, size=0,
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
