# Some code to look at how good the model predicts frequencies. Not very good. 

#Ok, let's the our real data values for Pleuni

prepMat <- matrix(data = NA, nrow = length(1:nrow(suppDat)), ncol = length(coef.vals))
for(i in 1:nrow(suppDat)){
    intercept <- 1
    TCG <- c(0,0,0)
    toReplace <- which(suppDat[i,]$WTnt == c("t", "c", "g"))
    if(length(toReplace) > 0){
        TCG[toReplace] <- 1
    }
    changeval <- bigChange[i]
    inRTval <- as.numeric(suppDat[i,]$unit == "RT")
    nonsynVal <- as.numeric(suppDat[i,]$TypeOfSite == "nonsyn")
    shapeVal <- suppDatShape[i,]$SHAPE    
    cpgval <- makesCpG[i]
    TCGns <- TCG*nonsynVal
    tCpG <- TCG[1] * cpgval
    nonsynCpG <- cpgval*nonsynVal
    nonsynCpGandT <- cpgval*nonsynVal*TCG[1]
    prepMat[i,] <-  c(intercept, TCG, changeval, inRTval, nonsynVal, shapeVal, cpgval, TCGns, tCpG, nonsynCpG, nonsynCpGandT)
}

colnames(prepMat) <- names(coef.vals)
t(t(prepMat[1:10,]))
names(coef.vals)

obsFreqs <- suppDat[-(1:40),]$colMeansTs0
modFreqs <- exp(prepMat %*% coef.vals)[-(1:40)]
cor.test(obsFreqs,modFreqs)

#Let's try looking at synonymous mutations at ancestrally As

layout(matrix(1:2, nrow = 1))
relInds <- intersect(which(suppDat$WTnt == "a"), which(suppDat$TypeOfSite == "syn"))
cpginds <- which(makesCpG[relInds] == 1)
colsp <- rep("black", length(relInds))
colsp[cpginds] <- "red"

plot(suppDat[relInds,]$colMeansTs0, exp(prepMat[relInds,] %*% coef.vals), col = colsp)
relInds <- intersect(which(suppDat$WTnt == "t"), which(suppDat$TypeOfSite == "syn"))
cpginds <- which(makesCpG[relInds] == 1)
colsp <- rep("black", length(relInds))
colsp[cpginds] <- "red"
plot(suppDat[relInds,]$colMeansTs0, exp(prepMat[relInds,] %*% coef.vals), col = colsp)

#Maybe nonsyn mutations (big change v no?)
require(RColorBrewer)
par(mar = c(4, 4, 1, 8))
plot(0, type = "n", xlim = c(0, .5), ylim = c(0, .041), xlab = "Observed frequency", ylab = "Modeled frequency")
abline(0, 1)
nucord <- c("a", "t", "c", "g")
pcols <- brewer.pal(8, "Paired")
for(i in  1:length(nucord)){
    relInds <- intersect(which(suppDat$WTnt == nucord[i]), which(suppDat$TypeOfSite == "nonsyn"))
    cpginds <- which(makesCpG[relInds] == 1)
    bigaachanges <- which(bigChange[relInds] == 1)
    pchp <- rep(1, length(relInds))
    pchp[cpginds] <- 2
    colsp <- rep(pcols[2*i - 1], length(relInds))
    colsp[bigaachanges] <- pcols[2*i]
    points(suppDat[relInds,]$colMeansTs0, exp(prepMat[relInds,] %*% coef.vals), col = colsp, pch = pchp, xlim = c(0, .025))
}
legend("bottomright", c("A", "A (big AA)", "T", "T (big AA)", "C", "C (big AA)", "G", "G (big AA)", "non-CpG-form", "CpG-form"), col = c(pcols, "black", "black"), pch = c(rep(16, 8), 1, 2), xpd = NA)
#legend(.043, .02, c("A", "A (big AA)", "T", "T (big AA)", "C", "C (big AA)", "G", "G (big AA)", "non-CpG-form", "CpG-form"), col = c(pcols, "black", "black"), pch = c(rep(16, 8), 1, 2), xpd = NA)

modelfits <- prepMat %*% coef.vals
outmat <- cbind(suppDat, exp(modelfits))
write.table(outmat, "../dat/dat.plus.modelfits.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

