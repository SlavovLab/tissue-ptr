load("gene_names.RData")
protein <- as.matrix(read.csv("protein.csv",row.names=1))
tab.mart <- read.csv("mart_export.txt",header=TRUE,stringsAsFactors=FALSE)
tab.mart <- tab.mart[!duplicated(tab.mart[,1]),]
rownames(tab.mart) <- tab.mart[,1]
tab.mart <- tab.mart[rownames(protein),]
go.names <- read.table("Go_names",comment.char="!",sep="\t",stringsAsFactors=FALSE,quote = "")
rownames(go.names) <- go.names[,1]

scientific_10 <- function(x) {
  parse(text=paste("10^", log10(x)))
}

sciNotation <- function(x, digits = 0) {
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
    tab[match(toupper(ensembl),tab$Ensembl_Gene),"RefseqRNA"]
}

r2e <- function(refseq){
    tab[match(toupper(refseq),tab$RefseqRNA),"Ensembl_Gene"]
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


myFilledContour <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
    length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
    col = color.palette(length(levels) - 1), plot.title, plot.axes, 
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
    axes = TRUE, vline=0.6, frame.plot = axes, hln=0, vln=0, ...) 
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
        stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2L]) * par("csi") * 2.54

    mar <- mar.orig
    mar[2] <- 5.1
    #mar[4L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp,cex.axis=1.5)

    .filled.contour(x, y, z, levels, col)

    abline(v=vln,col="dark grey",lwd=2)
    abline(h=hln,col="dark grey",lwd=2)
    
    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            Axis(x, side = 1, cex.axis=1.5)
            Axis(y, side = 2, cex.axis=1.5)
        }
    }
    else plot.axes
    if (frame.plot) 
        box()
    if (missing(plot.title)) 
        title(...)
    else plot.title
    invisible()
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
