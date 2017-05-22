#make AA transition figure


suppDat <- read.table("../dat/OverviewSelCoeffwPredictors.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

names(suppDat)
stotab <- table(suppDat$WTAA, suppDat$MUTAA)
allComps <- matrix(data = NA, nrow = sum(stotab > 0), ncol = 2)
for(i in 1:nrow(stotab)){
    for(j in 1:ncol(stotab)){
        if(stotab[i,j] != 0){ allComps[min(which(is.na(allComps[,1]))), ] <- c(rownames(stotab)[i], colnames(stotab)[j] ) }
    }
}

toMod <- c(0, -5, -10)
atvals <- c(toMod, toMod -14, toMod -28, toMod -42)

## Transparent colors
## modified from Mark Gardener 2015
## www.dataanalytics.org.uk
t.col <- function(color, percent = 50, name = NULL) {
    rgb.val <- col2rgb(color)
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = (100-percent)*255/100)
    return(t.col )
}


atvals
require(RColorBrewer)
oldcols <- brewer.pal(4, "Set2")
                                        #fix later
opac = 50
cols <- c(t.col(oldcols[1], opac), t.col(oldcols[2], opac),t.col(oldcols[3], opac),t.col(oldcols[4], opac))
pchs <- c(16, 17)
ylimval <- -13
pdf("../out/aachanges.pdf", height = 6, width = 12)
par(mar = c(5, 4, 1, 7))
plot(0, type = "n", xlim = c(1, sum(table(suppDat$WTAA, suppDat$MUTAA) > 0)), ylim = c(ylimval*4, 0), axes = FALSE, ylab = "Log Selection Coefficient", xlab = "Mutation")
axis(2, at = atvals, labels = rep("", length(atvals)), cex.axis = 1)
mtext(side = 2, at = atvals, text = rep(toMod, 4), line = 1, las = 2)
abline(v = 1:69, col = "grey90")
abline(h = atvals, col = "grey90")
abline(h = c(-13, -27, -41))
box()
for(i in 1:nrow(allComps)){
    relInds <- intersect(which(suppDat$WTAA == allComps[i, 1]), which(suppDat$MUTAA == allComps[i, 2]))
    relDat <- suppDat[relInds,]
relDat
    xjit <- rnorm(length(relInds), 0, .1)
    changesAA <- relDat$sameAAGroup + 1
    points(rep(i, length(relInds)) + xjit, relDat$logEstSelCoeff  + computeoffsets(relDat$TypeOfSite), pch = pchs[changesAA], col = cols[computeColInds(relDat$WTnt)], cex = .75)
#    points(rep(i, length(relInds)) + xjit, relDat$logEstSelCoeff  + computeoffsets(relDat$TypeOfSite), pch = 16, col = t.col(cols[computeColInds(relDat$WTnt)], 50), cex = .75)
}
labsize <- 1
mtext(side = 1, at = 1:69, text = paste(allComps[,2]), cex = labsize)
mtext("Mutant", 1, line = 0, at = -2, cex = labsize, adj = 1)
segments(-7 , -56.5, -2, -56.5, col="black", xpd=NA)
mtext("Ancestral", 1, line = 1, at = -2, cex = labsize, adj = 1)
for(i in unique(allComps[,1])){
    rangevals <- range(which(allComps[,1] == i))
    segments(rangevals[1] - .25, -56.5, rangevals[2] + .25, -56.5, col="black", xpd=NA)
    mtext(paste(i), 1, line = 1, at = mean(rangevals), cex = labsize)
}
mtext(4, at = c(-6, -20, -34, -47), text = c("Nonsyn", "Syn", "Stop", "Res."), line = .5)
legend(75, -10, c("A", "T", "C", "G"), col = c(cols, "black", "black"), pch = c(15, 15, 15, 15), xpd = NA, box.lty = "blank", title = "Ancestral\nnucleotide")
legend(75, -30, c("Yes", "No"), pch = pchs, xpd = NA, border = "white", box.lwd = 0, title = "Drastic AA\nchange", box.lty = "blank")
dev.off()

names(relDat)

allComps[43,]

computeColInds <- function(wtnt){
    colsToRet <- wtnt
    colsToRet[wtnt == "a"] <- 1
    colsToRet[wtnt == "t"] <- 2
    colsToRet[wtnt == "c"] <- 3
    colsToRet[wtnt == "g"] <- 4
    return(as.numeric(colsToRet))
}

computeoffsets <- function(types){
    offset <- types
    offset[types == "nonsyn"] <- 0
    offset[types == "syn"] <- -14
    offset[types == "stop"] <- -28
    offset[types == "res"] <- -42
    offset[types == "overlap"] <- -100 #do not plot overlap
    return(as.numeric(offset))
}

min(suppDat$logEstSelCoeff)











#make AA transition figure
#Just NS

suppDat <- read.table("../dat/OverviewSelCoeffwPredictors.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

suppDatNS <- suppDat[suppDat$TypeOfSite == "nonsyn",]
stotab <- table(suppDatNS$WTAA, suppDatNS$MUTAA)
allComps <- matrix(data = NA, nrow = sum(stotab > 0), ncol = 2)
for(i in 1:nrow(stotab)){
    for(j in 1:ncol(stotab)){
        if(stotab[i,j] != 0){ allComps[min(which(is.na(allComps[,1]))), ] <- c(rownames(stotab)[i], colnames(stotab)[j] ) }
    }
}


## Transparent colors
## modified from Mark Gardener 2015
## www.dataanalytics.org.uk
t.col <- function(color, percent = 50, name = NULL) {
    rgb.val <- col2rgb(color)
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = (100-percent)*255/100)
    return(t.col )
}

require(sfsmisc)
require(RColorBrewer)
oldcols <- brewer.pal(4, "Set2")
                                        #fix later
opac = 50
cols <- c(t.col(oldcols[1], opac), t.col(oldcols[2], opac),t.col(oldcols[3], opac),t.col(oldcols[4], opac))
pchs <- c(16, 17)
ylimval <- -12
pdf("../out/aachangesNS.pdf", height = 5, width = 9)
par(mar = c(5, 4, 1, 5))
plot(0, type = "n", xlim = c(1, nrow(allComps)), ylim = c(ylimval, 0), axes = FALSE, ylab = "Estimated Selection Coefficient (cost)", xlab = "Mutation")
axis(side = 2, at = seq(0, -12, by = -2), labels = expression(10^0, 10^-2, 10^-4, 10^-6, 10^-8, 10^-10, 10^-12), las = 2)
#mtext(side = 2, at = atvals, text = rep(toMod, 4), line = 1, las = 2)
abline(v = 1:nrow(allComps), col = "grey90")
abline(h = 0, col = "grey90")
box()
for(i in 1:nrow(allComps)){
    relInds <- intersect(which(suppDat$WTAA == allComps[i, 1]), which(suppDat$MUTAA == allComps[i, 2]))
    relDat <- suppDat[relInds,]
    xjit <- rnorm(length(relInds), 0, .1)
    changesAA <- relDat$sameAAGroup + 1
    points(rep(i, length(relInds)) + xjit, relDat$logEstSelCoeff  + computeoffsets(relDat$TypeOfSite), pch = pchs[changesAA], col = cols[computeColInds(relDat$WTnt)], cex = .75)
}
labsize <- 1
mtext(side = 1, at = 1:length(allComps[,2]), text = paste(allComps[,2]), cex = labsize)
mtext("Mutant", 1, line = 0, at = -1, cex = labsize, adj = 1)
segments(-6 , -13.15, -0.5, -13.15, col="black", xpd=NA)
mtext("Ancestral", 1, line = 1, at = -1, cex = labsize, adj = 1)
for(i in unique(allComps[,1])){
    rangevals <- range(which(allComps[,1] == i))
    segments(rangevals[1] - .25, -13.15, rangevals[2] + .25, -13.15, col="black", xpd=NA)
    mtext(paste(i), 1, line = 1, at = mean(rangevals), cex = labsize)
}
legend(49, -2, c("A", "T", "C", "G"), col = c(cols, "black", "black"), pch = c(15, 15, 15, 15), xpd = NA, box.lty = "blank", title = "Ancestral\nnucleotide")
legend(49, -8, c("Yes", "No"), pch = pchs, xpd = NA, border = "white", box.lwd = 0, title = "Drastic AA\nchange", box.lty = "blank")
dev.off()

names(relDat)




allComps[43,]

computeColInds <- function(wtnt){
    colsToRet <- wtnt
    colsToRet[wtnt == "a"] <- 1
    colsToRet[wtnt == "t"] <- 2
    colsToRet[wtnt == "c"] <- 3
    colsToRet[wtnt == "g"] <- 4
    return(as.numeric(colsToRet))
}

computeoffsets <- function(types){
    offset <- types
    offset[types == "nonsyn"] <- 0
    offset[types == "syn"] <- -14
    offset[types == "stop"] <- -28
    offset[types == "res"] <- -42
    offset[types == "overlap"] <- -100 #do not plot overlap
    return(as.numeric(offset))
}

min(suppDat$logEstSelCoeff)
