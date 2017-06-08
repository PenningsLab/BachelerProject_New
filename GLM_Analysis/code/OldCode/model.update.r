#read prepared data
read.csv("GLM_Analysis/dat/datFitModel.csv",row.names=1)->datFitModel

fullmodel.int <- glm(minor ~ t + c + g + bigAAChange + inRT + t*nonsyn + c*nonsyn + g*nonsyn + shape + CpG + CpG*t  + CpG*nonsyn + CpG*nonsyn*t,  family = "binomial", data = datFitModel[datFitModel$res == 0 & datFitModel$stop == 0,])

sumOfModel <- summary(fullmodel.int)

#Create a table for latex 
require(xtable)
#This doesn't work yet. 
#write(xtable(sumOfModel, digits = 3), "ModelTable.txt")
xtable(sumOfModel, digits = 3)

modcoef <- sumOfModel$coef

coef.vals <- modcoef[,1]
coef.pSE.vals <- coef.vals + modcoef[,2]
coef.mSE.vals <- coef.vals - modcoef[,2]

#colnames(datFitModel)
#t:CpG
#t:nonsyn:CpG


makeDataFrameToModify <- function(NsOrNo = 0, CpGorNo = 0, bigAAChangeOrNo = 0){

    atcg.mat <- as.data.frame(matrix(data = 0, ncol = ncol(datFitModel), nrow = 4))
    #ATCG elements
    diag(atcg.mat[1:4, 1:4]) <- 1
    #Reserve the first column for intercept
    atcg.mat[,1] <- 1
    #synonymous or nonsynonymous mutation?
    atcg.mat[,7] <- NsOrNo
    #bigAA change?
    atcg.mat[,6] <- bigAAChangeOrNo * atcg.mat[,7]
    #CpG mutation or not
    atcg.mat[1:2,12] <- CpGorNo
    atcg.mat[3:4,12] <- 0
    #nonysynonymous interactions with a t c g
    atcg.mat <- cbind(atcg.mat, atcg.mat[,2:4] * atcg.mat[,7], atcg.mat[,7]*atcg.mat[,12], atcg.mat[, 7] * atcg.mat[,2] , atcg.mat[, 7] * atcg.mat[,2] * atcg.mat[, 12])

#    makeDataFrameToModify(1, 1, 1) #PSP: I removed 4 lines. Not sure what they should do. Not present in model.r
#    rownames(modcoef)[1]
#    makeDataFrameToModify(1,1,1)[,rownames(modcoef)[c(1:17)]]
#    rownames(modcoef)

    #Pleuni: I uncommented this line. 
    colnames(atcg.mat) <- rownames(modcoef)
#    grep("nonsyn", )
#    atcg.mat
#    c('(Intercept)', names[2:(length(names))], paste(c("t", "c", "g"), ":nonsyn", sep = ""), "nonsyn:CpG", "t:CpG", "t:nonsyn:CpG" )

#Pleuni: not sure what this does ...    
#    names(atcg.mat) <- c('(Intercept)', names[2:(length(names))], paste(c("t", "c", "g"), ":nonsyn", sep = ""), "nonsyn:CpG", "t:CpG", "t:nonsyn:CpG" )

    return(as.matrix(atcg.mat))
}


makePlot <- function(main){

    plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(0, .15), axes = FALSE, ylab = "Predicted frequency of the mutation", xlab = "Ancestral nucleotide", main = main)
    axis(1, at = 1:4, c("A", "T", "C", "G"))
    axis(2)

    col.par <- "gray90"
    abline(h = seq(0, .15, by = .02), col = col.par)
    
    box()
}

plotVals <- function(NsOrNo, CpGorNo, bigAAChangeOrNo, colVal){

    pchpar <- 16
    cexpar <- 2


    setUpDat <- makeDataFrameToModify(NsOrNo, CpGorNo, bigAAChangeOrNo)[,rownames(modcoef)]
    
    #compute
    pointToPlot <- exp(setUpDat %*% coef.vals)
    arrowToPlotLow <- exp(setUpDat %*% coef.pSE.vals)
    arrowToPlotHigh <- exp(setUpDat %*% coef.mSE.vals)

    #plot
    points(1:4, pointToPlot, col = colVal,  pch = pchpar, cex = cexpar)
    arrows(1:4, arrowToPlotLow, 1:4, arrowToPlotHigh, length = .1, code = 3, angle = 90, lwd = 1.5, col = colVal)
    
}


makePlot.forS <- function(main){
    require(sfsmisc)
    plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(.00015, .02), axes = FALSE, ylab = "Predicted selection coefficient", xlab = "Ancestral nucleotide", main = main, log = "y")
    axis(1, at = 1:4, c("A", "T", "C", "G"))
    eaxis(2, at = 10^c(-2, -3, -4, -5))
    col.par <- "gray90"
    abline(h = (1:9)*10^(-2), col = col.par)
    abline(h = (1:9)*10^(-3), col = col.par)
    abline(h = (1:9)*10^(-4), col = col.par)
    
    box()
}

plotVals.svals <- function(NsOrNo, CpGorNo, bigAAChangeOrNo, colVal){

    pchpar <- 16
    cexpar <- 2
    mus <- c(1.11e-05, 1.11e-05, 2.41e-05, 5.48e-05)
        
    setUpDat <- makeDataFrameToModify(NsOrNo, CpGorNo, bigAAChangeOrNo)[,rownames(modcoef)]

    #compute
    pointToPlot <- exp(setUpDat %*% coef.vals)
    arrowToPlotLow <- exp(setUpDat %*% coef.pSE.vals)
    arrowToPlotHigh <- exp(setUpDat %*% coef.mSE.vals)

    #plot
    points(1:4, mus/pointToPlot, col = colVal,  pch = pchpar, cex = cexpar)
    arrows(1:4, mus/arrowToPlotLow, 1:4, mus/arrowToPlotHigh, length = .1, code = 3, angle = 90, lwd = 1.5, col = colVal)
    
}

require(RColorBrewer)
pdf("GLM_Analysis/graphs/modeled_freqs2017May.pdf", width = 12, height = 7)

cols <- brewer.pal(3, "Set2")
layout(matrix(1:2, nrow = 1))
par(mar = c(4, 4.5, 1.5, 1))
makePlot(main = "Synonymous Sites")
plotVals(0, 1, 0, cols[1])
plotVals(0, 0, 0, cols[2])
legend("topleft", c("CpG-forming", "non-CpG-forming"), col = cols[1:2], pch = 16, bg = "white")
makePlot(main = "Non-synonymous Sites")
plotVals(1, 1, 0, cols[1])
plotVals(1, 0, 0, cols[2])
plotVals(1, 1, 1, cols[3])
legend("topleft", c("CpG-forming", "non-CpG-forming", "Major AA change"), col = cols, pch = 16, bg = "white" )

dev.off()

require(RColorBrewer)
pdf("GLM_Analysis/graphs/modeled_selcoefs2017May.pdf", width = 12, height = 7)
cols <- brewer.pal(3, "Set2")
layout(matrix(1:2, nrow = 1))
par(mar = c(4, 4.5, 1.5, 1))
makePlot.forS(main = "Synonymous Sites")
plotVals.svals(0, 1, 0, cols[1])
plotVals.svals(0, 0, 0, cols[2])
legend("topleft", c("CpG-forming", "non-CpG-forming"), col = cols[1:2], pch = 16, bg = "white")
makePlot.forS(main = "Non-synonymous Sites")
plotVals.svals(1, 1, 0, cols[1])
plotVals.svals(1, 0, 0, cols[2])
plotVals.svals(1, 1, 1, cols[3])
legend("topleft", c("CpG-forming", "non-CpG-forming", "Major AA change"), col = cols, pch = 16, bg = "white" )
dev.off()


toPlot <- exp(makeDataFrameToModify(1, 1, 0)[,rownames(modcoef)] %*% as.matrix(modcoef[,1]))
points(1:4, toPlot, col = cols[1], pch = pchpar, cex = cexpar)
toPlot <- exp(makeDataFrameToModify(1, 0, 0)[,rownames(modcoef)] %*% as.matrix(modcoef[,1]))
points(1:4, toPlot, col = cols[2], pch = pchpar, cex = cexpar)
toPlot <- exp(makeDataFrameToModify(1, 0, 1)[,rownames(modcoef)] %*% as.matrix(modcoef[,1]))
points(1:4, toPlot, col = cols[3], pch = pchpar, cex = cexpar)

#Pleuni commented this
#modcoef[,1]
#exp(as.matrix(syn.atcg[,names(modcoef)]) %*% (modcoef[,1]))
#syn.atcg[,]
#rownames(modcoef)[1]
#dim(datFitModel)


listofall <- list()
for(i in 1:length(suppDat$num)){

    print(i)
    sto <- as.data.frame(t(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "ID")))
    sto[,c(1:13)] <- as.numeric(sto[,c(1:13)])
    sto[,14] <- as.factor(sto[,14])
    names(sto) <- names

    pos <- suppDat[i,]$num
    inRT <- as.numeric(pos < 298)
    atcg <- c(0,0,0,0)
    atcg[which(nucord == suppDat[i,]$WTnt)] <- 1
    bigAAChange <- bigChange[i]    
    nonsyn <- as.numeric(regexpr("nonsyn",suppDat[i,]$TypeOfSite) > 0)
    stop <- as.numeric(regexpr("stop",suppDat[i,]$TypeOfSite) > 0)
    res <- as.numeric(regexpr("res",suppDat[i,]$TypeOfSite) > 0)
    pairingprob <- suppDat[i,]$PAIRINGPROB
    shape <- suppDat[i,]$SHAPE
    CpG <- makesCpG[i]
    pref <- c(atcg, inRT, bigAAChange, nonsyn, stop, res, shape, pairingprob, CpG)


#pref is the same for all patients. 
#Now, let's go through each of the patients and do this:
    match.is <- which(dat[,3 ] == i)
    unique(dat[,2])
    for(j in unique(dat[,2])){
        toExp <- dat[intersect(which(dat[,2] == j), match.is),]
        minorNum <- toExp$MutCount
        majorNum <- toExp$WTcount - minorNum
        toReturn <- rbind(repMe(pref, 1, minorNum), repMe(pref, 0, majorNum))
        if(!is.null(toReturn)){
            toReturn <- data.frame(toReturn, rep(j, nrow(toReturn)))
            names(toReturn) <- names
            sto <- rbind(sto, toReturn)
        }
    }
#    allDat <- rbind(allDat, sto[-1,])
    listofall[[i]] <- sto[-1,]
}



ktable(fulldat[,13])

allDat <- as.data.frame(t(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,  "ID")))
allDat[,c(1:13)] <- as.numeric(allDat[,c(1:13)])
allDat[,14] <- as.factor(allDat[,14])
names(allDat) <- names
for(i in 1:200){
    print(i)
    allDat <- rbind(allDat, listofall[[i]])
}
allDat1 <- allDat[-1,]
print("1done")
allDat <- as.data.frame(t(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,  "ID")))
allDat[,c(1:13)] <- as.numeric(allDat[,c(1:13)])
allDat[,14] <- as.factor(allDat[,14])
names(allDat) <- names
for(i in 201:400){
    print(i)
    allDat <- rbind(allDat, listofall[[i]])
}
allDat2 <- allDat[-1,]
print("2done")
allDat <- as.data.frame(t(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,  "ID")))
allDat[,c(1:13)] <- as.numeric(allDat[,c(1:13)])
allDat[,14] <- as.factor(allDat[,14])
names(allDat) <- names
for(i in 401:600){
    print(i)
    allDat <- rbind(allDat, listofall[[i]])
}
allDat3 <- allDat[-1,]
print("3done")
allDat <- as.data.frame(t(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,  "ID")))
allDat[,c(1:13)] <- as.numeric(allDat[,c(1:13)])
allDat[,14] <- as.factor(allDat[,14])
names(allDat) <- names
for(i in 601:800){
    print(i)
    allDat <- rbind(allDat, listofall[[i]])
}
allDat4 <- allDat[-1,]
print("4done")
allDat <- as.data.frame(t(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,  "ID")))
allDat[,c(1:13)] <- as.numeric(allDat[,c(1:13)])
allDat[,14] <- as.factor(allDat[,14])
names(allDat) <- names
for(i in 801:984){
    print(i)
    allDat <- rbind(allDat, listofall[[i]])
}
allDat5 <- allDat[-1,]


allDat <- rbind(allDat1, allDat2, allDat3, allDat4, allDat5)

names(allDat)
fullmodel.int <- glm(minor ~ t+c+g+ bigAAChange + inRT + t*nonsyn + c*nonsyn + g*nonsyn + shape + CpG + pairingprob,  family = "binomial", data = allDat[allDat$res == 0 & allDat$stop == 0,])
coef.save <- summary(fullmodel.int)$coef

require(xtable)
xtable(coef.save[, c(1, 2, 4)])

table(allDat[,'minor'])




names = c( "a", "t", "c", "g", "inRT", "bigAAChange", "nonsyn", "stop", "res", "shape", "pairingprob","minor", "ID")

nucord <- c("a", "t", "c", "g")

listofall <- list()
#for(i in 1:length(suppDat$num)){
for(i in 1:length(suppDat$num)){
    print(i)
    sto <- as.data.frame(t(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "ID")))
    sto[,c(1:12)] <- as.numeric(sto[,c(1:12)])
    sto[,13] <- as.factor(sto[,13])
    names(sto) <- c( "a", "t", "c", "g", "inRT", "bigAAChange", "nonsyn", "stp", "res", "shape", "pairingprob", "minor", "ID")
    pos <- suppDat[i,]$num
    inRT <- as.numeric(pos < 298)
    atcg <- c(0,0,0,0)
    atcg[which(nucord == suppDat[i,]$WTnt)] <- 1
    bigAAChange <- bigChange[i]    
    nonsyn <- as.numeric(regexpr("nonsyn",suppDat[i,]$TypeOfSite) > 0)
    stop <- as.numeric(regexpr("stop",suppDat[i,]$TypeOfSite) > 0)
    res <- as.numeric(regexpr("res",suppDat[i,]$TypeOfSite) > 0)
    pairingprob <- suppDat[i,]$PAIRINGPROB
    shape <- suppDat[i,]$SHAPE
    pref <- c(atcg, inRT, bigAAChange, nonsyn, stop, res, shape, pairingprob)
#pref is the same for all patients. 
#Now, let's go through each of the patients and do this:
    match.is <- which(dat[,3 ] == i)
    for(j in unique(dat[,2])){
        toExp <- dat[intersect(which(dat[,2] == j), match.is),]
        minorNum <- toExp$MutCount
        majorNum <- toExp$WTcount - minorNum
        toReturn <- rbind(repMe(pref, 1, minorNum), repMe(pref, 0, majorNum))
        if(!is.null(toReturn)){
            toReturn <- data.frame(toReturn, rep(j, nrow(toReturn)))
            names(toReturn) <- c( "a", "t", "c", "g", "inRT", "bigAAChange", "nonsyn", "stp", "res", "shape", "pairingprob", "minor", "ID")
            sto <- rbind(sto, toReturn)
        }
    }
#    allDat <- rbind(allDat, sto[-1,])
    listofall[[i]] <- sto[-1,]
}




allDat <- as.data.frame(t(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "ID")))
allDat[,c(1:12)] <- as.numeric(allDat[,c(1:12)])
allDat[,13] <- as.factor(allDat[,13])
names(allDat) <- c( "a", "t", "c", "g", "inRT", "bigAAChange", "nonsyn", "stp", "res", "shape", "pairingprob", "minor", "ID")
for(i in 1:250){
    print(i)
    allDat <- rbind(allDat, listofall[[i]])
}
allDat1 <- allDat[-1,]
print("1done")
allDat <- as.data.frame(t(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "ID")))
allDat[,c(1:12)] <- as.numeric(allDat[,c(1:12)])
allDat[,13] <- as.factor(allDat[,13])
names(allDat) <- c( "a", "t", "c", "g", "inRT", "bigAAChange", "nonsyn", "stp", "res", "shape", "pairingprob", "minor", "ID")
for(i in 251:500){
    print(i)
    allDat <- rbind(allDat, listofall[[i]])
}
allDat2 <- allDat[-1,]
print("2done")
allDat <- as.data.frame(t(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "ID")))
allDat[,c(1:12)] <- as.numeric(allDat[,c(1:12)])
allDat[,13] <- as.factor(allDat[,13])
names(allDat) <- c( "a", "t", "c", "g", "inRT", "bigAAChange", "nonsyn", "stp", "res", "shape", "pairingprob", "minor", "ID")
for(i in 501:750){
    print(i)
    allDat <- rbind(allDat, listofall[[i]])
}
allDat3 <- allDat[-1,]
print("3done")
allDat <- as.data.frame(t(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "ID")))
allDat[,c(1:12)] <- as.numeric(allDat[,c(1:12)])
allDat[,13] <- as.factor(allDat[,13])
names(allDat) <- c( "a", "t", "c", "g", "inRT", "bigAAChange", "nonsyn", "stp", "res", "shape", "pairingprob", "minor", "ID")
for(i in 751:984){
    print(i)
    allDat <- rbind(allDat, listofall[[i]])
}
allDat4 <- allDat[-1,]
print("4done")

allDat <- rbind(allDat1, allDat2, allDat3, allDat4)
#allDatRT2 <- allDat
#allDatRT1 <- allDat
#allDatPR <- allDat
#allDat <- rbind(allDatPR, allDatRT1, allDatRT2)

predict( tmp, allDat[1:1000,]) 

fullmodel <- glm(minor ~ t+c+g+bigAAChange + inRT + nonsyn + stp,  family = "binomial", data = allDat[allDat$res == 0,])
fullmodel.int <- glm(minor ~ t+c+g+ bigAAChange + inRT + t*nonsyn + c*nonsyn + g*nonsyn + shape,  family = "binomial", data = allDat[allDat$res == 0 & allDat$stp == 0,])
summary(fullmodel.int)

Gs <- intersect(intersect(intersect(which(allDat$res == 0), which(allDat$stp == 0)), which(allDat$nonsyn == 1)), which(allDat$g == 1))
tab <- table(allDat[Gs, ]$minor, allDat[Gs,]$ID)
hist(tab[2,]/tab[1,], breaks = seq(0, .1, by = .01))

summary(fullmodel.int)
fakeObs <- as.data.frame(t(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "ID")))
fakeObs[,c(1:13)] <- as.numeric(fakeObs[,c(1:13)])
fakeObs[,14] <- as.factor(fakeObs[,14])
fakeObs[1,c(1:13)] <- c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 1)
fakeObs[2,c(1:13)] <- c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1)
fakeObs[3,c(1:13)] <- c(0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1)
fakeObs[4,c(1:13)] <- c(0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)
fakeObs[5,c(1:13)] <- c(1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1)
fakeObs[6,c(1:13)] <- c(0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1)
fakeObs[7,c(1:13)] <- c(0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1)
fakeObs[8,c(1:13)] <- c(0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1)
fakeObs[9,c(1:13)] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
fakeObs[10,c(1:13)] <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
fakeObs[11,c(1:13)] <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
fakeObs[12,c(1:13)] <- c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1)
names(fakeObs) <- c( "a", "t", "c", "g", "inRT", "bigAAChange", "nonsyn", "stp", "res", "minor", "(Intercept)", "ID")

fakeObs[,'t:nonsyn'] <- fakeObs[,'t'] * fakeObs[,'nonsyn']
fakeObs[,'c:nonsyn'] <- fakeObs[,'c'] * fakeObs[,'nonsyn']
fakeObs[,'g:nonsyn'] <- fakeObs[,'g'] * fakeObs[,'nonsyn']



fm.coef <- summary(fullmodel)$coef
fm.coef.int <- summary(fullmodel.int)$coef
fm.coef.int

round(fm.coef.int, 3)

atcg.p <- exp(as.matrix(fakeObs[,rownames(fm.coef.int)]) %*% (fm.coef.int[,1]))
atcg.high <- exp(as.matrix(fakeObs[,rownames(fm.coef.int)]) %*% (fm.coef.int[,1] + 1.96*fm.coef.int[,2]))
atcg.low <- exp(as.matrix(fakeObs[,rownames(fm.coef.int)]) %*% (fm.coef.int[,1] - 1.96*fm.coef.int[,2]))

atcg.ns <- atcg.p[1:4]
atcg.bAAc <- atcg.p[5:8]
atcg.syn <- atcg.p[9:12]
atcg.ns.low <- atcg.low[1:4]
atcg.bAAc.low <- atcg.low[5:8]
atcg.syn.low <- atcg.low[9:12]
atcg.ns.high <- atcg.high[1:4]
atcg.bAAc.high <- atcg.high[5:8]
atcg.syn.high <- atcg.high[9:12]

require(sfsmisc)




colsNS <- brewer.pal(3, name = "Set2")
mu.rates <- c(1.11047e-05, 1.10618e-05, 2.41358e-05, 5.47833e-05)
#plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(0.0001, .013), axes = FALSE, xlab = "Reference nucleotide", ylab = "mean estimated s")
#axis(2)
plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(0.0002, .05), axes = FALSE, xlab = "Nucleotide", ylab = "mean estimated s", log = "y", bg = "grey")
eaxis(2, at = 10^c(-2, -3, -4, -5))
axis(1, 1:4, c("A", "T", "C", "G"))
box()
points(1:4, mu.rates/atcg.bAAc, pch = "-", col = colsNS[1])
arrows(1:4, mu.rates/atcg.bAAc.high, 1:4, mu.rates/atcg.bAAc.low, length = 0, col = colsNS[1], lwd = 2)
points(1:4, mu.rates/atcg.ns, pch = "-", col = colsNS[2])
arrows(1:4, mu.rates/atcg.ns.high, 1:4, mu.rates/atcg.ns.low, length = 0, col = colsNS[2], lwd = 2)
points(1:4, mu.rates/atcg.syn, pch = "-", col = colsNS[3])
arrows(1:4, mu.rates/atcg.syn.high, 1:4, mu.rates/atcg.syn.low, length = 0, col = colsNS[3], lwd = 2)
legend("topleft", c("Major AA change (Non-syn)", "No major AA change (Non-syn)", "Synonymous"), col = colsNS, lwd = 2)


mu.rates/
    atcg.syn
mu.rates

#no interactions

#predict
atcg.p <- exp(as.matrix(fakeObs[,rownames(fm.coef)]) %*% (fm.coef[,1]))
atcg.high <- exp(as.matrix(fakeObs[,rownames(fm.coef)]) %*% (fm.coef[,1] + 1.96*fm.coef[,2]))
atcg.low <- exp(as.matrix(fakeObs[,rownames(fm.coef)]) %*% (fm.coef[,1] - 1.96*fm.coef[,2]))

atcg.ns <- atcg.p[1:4]
atcg.bAAc <- atcg.p[5:8]
atcg.syn <- atcg.p[9:12]
atcg.ns.low <- atcg.low[1:4]
atcg.bAAc.low <- atcg.low[5:8]
atcg.syn.low <- atcg.low[9:12]
atcg.ns.high <- atcg.high[1:4]
atcg.bAAc.high <- atcg.high[5:8]
atcg.syn.high <- atcg.high[9:12]

table(suppDat$TSmutrate, suppDat$WTnt)

require(sfsmisc)

if (FALSE){
    pdf("GLM_Analysis/graphs/modelfits.atcg2017May.pdf", height = 5, width = 5)
    colsNS <- brewer.pal(3, name = "Set2")
    par(mar = c(4, 4.5, 1, 1))
    mu.rates <- c(1.11047e-05, 1.10618e-05, 2.41358e-05, 5.47833e-05)
    #plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(0.0001, .006), axes = FALSE, xlab = "Reference nucleotide", ylab = "mean estimated s")
    #axis(2)
    plot(0, type = "n", xlim = c(.75, 4.25), ylim = c(0.00025, .015), axes = FALSE, xlab = "Nucleotide", ylab = "mean estimated s", log = "y", cex.lab = 1.5)
    eaxis(2, at = 10^c(-2, -3, -4, -5))
    axis(1, 1:4, c("A", "T", "C", "G"))
    lwd.par <- 3
    cex.par <- 2
    col.par <- "gray90"
    abline(h = (1:9)*10^(-2), col = col.par)
    abline(h = (1:9)*10^(-3), col = col.par)
    abline(h = (1:9)*10^(-4), col = col.par)
    box()
    points(1:4, mu.rates/atcg.bAAc, pch = "-", col = colsNS[1], cex = cex.par)
    arrows(1:4, mu.rates/atcg.bAAc.high, 1:4, mu.rates/atcg.bAAc.low, length = 0, col = colsNS[1], lwd = lwd.par)
    points(1:4, mu.rates/atcg.ns, pch = "-", col = colsNS[2], cex = cex.par)
    arrows(1:4, mu.rates/atcg.ns.high, 1:4, mu.rates/atcg.ns.low, length = 0, col = colsNS[2], lwd = lwd.par)
    points(1:4, mu.rates/atcg.syn, pch = "-", col = colsNS[3], cex = cex.par)
    arrows(1:4, mu.rates/atcg.syn.high, 1:4, mu.rates/atcg.syn.low, length = 0, col = colsNS[3], lwd = lwd.par)
    legend("topleft", c("Group AA change", "No-group AA change", "Synonymous"), col = colsNS, lwd = lwd.par, bg = "white")
    dev.off()
}
    
#allDatnew <- allDat
#allDatold <- allDat
#fullalldat <- allDatpallDat <- rbind(allDatold[-1,], allDatnew[-1,])

mean(allDat[allDat$a == 1 & allDat$nonsyn == 1,]$minor)
mean(allDat[allDat$t == 1 & allDat$nonsyn == 1,]$minor)
mean(allDat[allDat$c == 1 & allDat$nonsyn == 1,]$minor)
mean(allDat[allDat$g == 1 & allDat$nonsyn == 1,]$minor)

summary(tmp)

table(allDat$minor)
table(allDat$stp)

colnames(sto) <- c( "a", "t", "c", "g", "inRT", "bigAAChange", "nonsyn", "stop", "res", "minor", "ID")
#Ok, sto is really big
processed.dat <- as.data.frame(sto[-1,])

tmp <- glm(minor ~ nonsyn + stop + res,  data = processed.dat)
plot(tmp)

table(processed.dat$nonsyn)


sto <-  expInf(1)
colnames(sto) <- c("patno", "a", "t", "c", "g", "inRT", "bigAAChange", "nonsyn", "stop", "minor")
for(i in 2:nrow(dat)){
    if(i %% 1000 == 0){ print(paste(i ,"out of ",nrow(dat))) }
    tmp <- expInf(i)
    if(!is.null(tmp)){
        sto <- rbind(sto, tmp)
    }
}


layout(matrix(1:4, nrow = 2))
hist(suppDat[suppDat$TypeOfSite == "nonsyn",]$SHAPE, main = "nonsyn", xlab = "", xlim = c(0, 2.5))
abline(v= mean(suppDat[suppDat$TypeOfSite == "nonsyn",]$SHAPE), lwd = 2, col = "red")
hist(suppDat[suppDat$TypeOfSite == "syn",]$SHAPE, main = "syn", xlab = "", xlim = c(0, 2.5))
abline(v= mean(suppDat[suppDat$TypeOfSite == "syn",]$SHAPE), lwd = 2, col = "red")
hist(suppDat[suppDat$TypeOfSite == "stop",]$SHAPE, main = "stop", xlab = "", xlim = c(0, 2.5))
abline(v= mean(suppDat[suppDat$TypeOfSite == "stop",]$SHAPE), lwd = 2, col = "red")
hist(suppDat[suppDat$TypeOfSite == "res",]$SHAPE, main = "res", xlab = "", xlim = c(0, 2.5))
abline(v= mean(suppDat[suppDat$TypeOfSite == "res",]$SHAPE), lwd = 2, col = "red")
