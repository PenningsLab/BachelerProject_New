#make AA transition figure
#Just NS
source("Rscripts/baseRscript.R")

computeColInds <- function(wtnt){
    colsToRet <- wtnt
    colsToRet[wtnt == "a"] <- 1
    colsToRet[wtnt == "t"] <- 2
    colsToRet[wtnt == "c"] <- 3
    colsToRet[wtnt == "g"] <- 4
    return(as.numeric(colsToRet))
}

Bach.dat<-read.csv("Output/OverviewSelCoeff_BachelerFilter.csv",stringsAsFactors = FALSE)

Bach.datNS <- Bach.dat[Bach.dat$TypeOfSite == "nonsyn",]
stotab <- table(Bach.datNS$WTAA, Bach.datNS$MUTAA)
allComps <- matrix(data = NA, nrow = sum(stotab > 0), ncol = 2) #Get a list of all AA substitutions in the dataset
for(i in 1:nrow(stotab)){
    for(j in 1:ncol(stotab)){
        if(stotab[i,j] != 0){ allComps[min(which(is.na(allComps[,1]))), ] <- c(rownames(stotab)[i], colnames(stotab)[j] ) }
    }
}

opac = 50
pchs <- c(16, 17)
ylimval <- -4
if (TRUE){
pdf("Output/aachangesNS_2017Nov.pdf", height = 5, width = 9)
par(mar = c(5, 4, 1, 5))
plot(0, type = "n", xlim = c(1, nrow(allComps)), ylim = c(ylimval, 0), axes = FALSE, ylab = "Estimated Selection Coefficient (cost)", xlab = "Mutation")
axis(side = 2, at = seq(0, -4, by = -1), labels = expression(10^0, 10^-1, 10^-2, 10^-3, 10^-4), las = 2)
#mtext(side = 2, at = atvals, text = rep(toMod, 4), line = 1, las = 2)
abline(v = 1:nrow(allComps), col = "grey90")
abline(h = 0, col = "grey90")
box()
for(i in 1:nrow(allComps)){
    relInds <- intersect(which(Bach.datNS$WTAA == allComps[i, 1]), which(Bach.datNS$MUTAA == allComps[i, 2]))
    relDat <- Bach.datNS[relInds,]
    xjit <- rnorm(length(relInds), 0, .1)
    changesAA <- relDat$sameAAGroup + 1
    # changed relDat$logEstSelCoeff  into log10(relDat$EstSelCoeff)
    points(rep(i, length(relInds)) + xjit, log10(relDat$EstSelCoeff), pch = pchs[relDat$bigAAChange[1]+1], col = cols[computeColInds(relDat$WTnt)], cex = .75)
}
labsize <- .9
mtext(side = 1, at = 1:length(allComps[,2]), text = paste(allComps[,2]), cex = labsize*.8)
mtext("Mutant", 1, line = 0, at = -1, cex = labsize, adj = 1)
linepos=4.4
segments(-6 , -linepos, -0.5, -linepos, col="black", xpd=NA)
mtext("Ancestral", 1, line = 1, at = -1, cex = labsize, adj = 1)
for(i in unique(allComps[,1])){
    rangevals <- range(which(allComps[,1] == i))
    segments(rangevals[1] - .25, -linepos, rangevals[2] + .25, -linepos, col="black", xpd=NA)
    mtext(paste(i), 1, line = 1, at = mean(rangevals), cex = labsize)
}
legend(49, -1, c("A", "T", "C", "G"), col = c(cols, "black", "black"), pch = c(15, 15, 15, 15), xpd = NA, box.lty = "blank", title = "Ancestral\nnucleotide")
legend(49, -2.5, c("No", "Yes"), pch = pchs, xpd = NA, border = "white", box.lwd = 0, title = "Drastic AA\nchange", box.lty = "blank")

dev.off()
}

