dat <- read.table("../dat/CountDataLehmanRT.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)
suppDat <- read.table("../dat/OverviewSelCoeffwProteinFeatures.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
suppDatShape <- read.table("../dat/OverviewSelCoeffwSHAPE.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

table(dat[,2])
summary(dat[,6] + dat[,5])


#This function is amCat old in other files, but they are functionally the same on this data.
#function to judge amino acid categories
pos <- "R|H|K"
neg <- "D|E"
unc <- "S|T|N|Q"
spe <- "C|U|G|P"
hyd <- "A|I|L|F|M|W|Y|V"
amCat <- function(AA){
    if(regexpr(pos, AA) > 0){ return(0) }
    if(regexpr(neg, AA) > 0){ return(1) }
    if(regexpr(unc, AA) > 0){ return(2) }
    if(regexpr(spe, AA) > 0){ return(3) }
    if(regexpr(hyd, AA) > 0){ return(4) }
    return(5)
}

#Which positions will have big amino acid changes if they undergo transitions
bigChange <- rep(NA, nrow(suppDat))
for(i in 1:nrow(suppDat)){
  WT <- amCat(suppDat[i,'WTAA'])
  MUT <- amCat(suppDat[i,'MUTAA'])
  if(WT == MUT){ bigChange[i] <- 0
             }else{
                 bigChange[i] <- 1
             }
}

#Which positions will create CpG sites if they undergo transitions
makesCpG <- rep(0, length = nrow(suppDat))
for(i in 1:nrow(suppDat)){
    trip <- suppDat$WTnt[c(i-1, i, i + 1)]
    if(trip[1] == "c" & trip[2] == "a" ){
        makesCpG[i] <- 1
    }
    if(trip[2] == "t" & trip[3] == "g"){
        makesCpG[i] <- 1
    }
}

#Helper function to do some matrix replication in creating the dataset
repMe <- function(pref, majormin, toRep){
    if(toRep == 0){ return() }
    toRet <- matrix(data = NA, nrow = toRep, ncol = length(pref) + 1)
    for(i in 1:toRep){
        toRet[i,] <- c(pref, majormin)
    }
    return(toRet)
}

names = c( "a", "t", "c", "g", "inRT", "bigAAChange", "nonsyn", "stop", "res", "shape", "pairingprob", "CpG", "helix", "alpha", "beta", "coil", "minor")
nucord <- c("a", "t", "c", "g")


#numRows <- sum(dat$WTcount) + sum(dat$MutCount)
datRows <- matrix(data = NA, nrow = 141405456, ncol = length(names))
colnames(datRows) <- names
relInd <- 1
dsnum <- 350
#Leave out the first 40 positions
for(i in 297:length(suppDat$num)){
    print(i)
    pos <- suppDat[i,]$num
    inRT <- as.numeric(pos < 298)
    atcg <- c(0,0,0,0)
    atcg[which(nucord == suppDat[i,]$WTnt)] <- 1
    bigAAChange <- bigChange[i]    
    nonsyn <- as.numeric(regexpr("nonsyn",suppDat[i,]$TypeOfSite) > 0)
    stop <- as.numeric(regexpr("stop",suppDat[i,]$TypeOfSite) > 0)
    res <- as.numeric(regexpr("res",suppDat[i,]$TypeOfSite) > 0)
    pairingprob <- suppDatShape[i,]$PAIRINGPROB
    shape <- suppDatShape[i,]$SHAPE
    CpG <- makesCpG[i]
    proteinfeats <- c(0,0,0,0)
    proteinfeats[which(suppDat[i,]$Protein_feature == c("3-10-helix", "alpha", "beta", "coil"))] <- 1
    pref <- c(atcg, inRT, bigAAChange, nonsyn, stop, res, shape, pairingprob, CpG, proteinfeats)
#pref is the same for all patients. 
#Now, let's go through each of the patients and do this:
    match.is <- which(dat[,3 ] == i - 296)
    for(j in unique(dat[,2])){
        toExp <- dat[intersect(which(dat[,2] == j), match.is),]
        if(nrow(toExp) == 0){
            toReturn = NULL
        }else{
            minorNum <- toExp$MutCount
            majorNum <- toExp$WTcount - minorNum
            if(minorNum + majorNum > dsnum){
                minorNumNew <- rbinom(1, dsnum, c(minorNum)/sum(minorNum, majorNum))
                majorNumNew <- dsnum - minorNumNew
                toReturn <- rbind(repMe(pref, 1, minorNumNew), repMe(pref, 0, majorNumNew))
            }else{
                toReturn <- rbind(repMe(pref, 1, minorNum), repMe(pref, 0, majorNum))
            }
        }
        if(!is.null(toReturn)){
            datRows[relInd:(relInd + nrow(toReturn) - 1),] <- toReturn
            relInd <- relInd + nrow(toReturn)
        }
    }
}


NAvals <- which(!is.na(datRows[,1]))
datFitModel <- as.data.frame(datRows[NAvals,])
dim(datFitModel)
## toExc <- length(NAvals):152735204
## #datFitModel <- as.data.frame(datRows[NAvals,])
## datFitModel.notDF <- datRows[-toExc,]
## datRows <- ""
## datFitModel <- as.data.frame(datFitModel.notDF)
## datFitModel.notDF <- ""

write.table(datFitModel, "LehmanToFitSmall.txt", quote = FALSE, row.names= FALSE, col.names = TRUE)

                                        #141405456

#Substantially more complicated model with structural elements
#fullmodel.int <- glm(minor ~ t + c + g + bigAAChange + inRT + t*nonsyn + c*nonsyn + g*nonsyn + shape + CpG + CpG*t  + CpG*nonsyn + CpG*nonsyn*t + helix*nonsyn + beta*nonsyn + coil*nonsyn,  family = "binomial", data = datFitModel[datFitModel$res == 0 & datFitModel$stop == 0,])

datFitModel <- as.data.frame(datRows)


fullmodel.int <- glm(minor ~ t + c + g + bigAAChange + t*nonsyn + c*nonsyn + g*nonsyn + CpG + CpG*t  + CpG*nonsyn + CpG*nonsyn*t,  family = "binomial", data = datFitModel[datFitModel$res == 0 & datFitModel$stop == 0,])
#fullmodel.int <- glm(minor ~ t + c + g + bigAAChange + inRT + t*nonsyn + c*nonsyn + g*nonsyn + shape + CpG + CpG*t  + CpG*nonsyn + CpG*nonsyn*t,  family = "binomial", data = datFitModel[datFitModel$res == 0 & datFitModel$stop == 0,])

sumOfModel <- summary(fullmodel.int)
modcoef <- sumOfModel$coef
coef.vals <- modcoef[,1]
coef.pSE.vals <- coef.vals + modcoef[,2]
coef.mSE.vals <- coef.vals - modcoef[,2]


makeDataFrameToModify <- function(nsOrNo = 0, CpGorNo = 0, bigAAChangeOrNo = 0){

    atcg.mat <- as.data.frame(matrix(data = 0, ncol = ncol(datFitModel), nrow = 4))
    #ATCG elements
    diag(atcg.mat[1:4, 1:4]) <- 1
    #Reserve the first column for intercept
    atcg.mat[,1] <- 1
    #synonymous or nonsynonymous mutation?
    atcg.mat[,7] <- nsOrNo
    #bigAA change?
    atcg.mat[,6] <- bigAAChangeOrNo * atcg.mat[,7]
    #CpG mutation or not
    atcg.mat[1:2,12] <- CpGorNo
    atcg.mat[3:4,12] <- 0
    #nonysynonymous interactions with a t c g
    atcg.mat <- cbind(atcg.mat, atcg.mat[,2:4] * atcg.mat[,7], atcg.mat[,7]*atcg.mat[,12], atcg.mat[, 7] * atcg.mat[,2] , atcg.mat[, 7] * atcg.mat[,2] * atcg.mat[, 12])

    names(atcg.mat) <- c('(Intercept)', names[2:(length(names))], paste(c("t", "c", "g"), ":nonsyn", sep = ""), "nonsyn:CpG", "t:CpG", "t:nonsyn:CpG" )

    return(as.matrix(atcg.mat))
}


makePlot <- function(main){

    plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(0, 1), axes = FALSE, ylab = "Predicted frequency of the mutation", xlab = "Ancestral nucleotide", main = main)
    axis(1, at = 1:4, c("A", "T", "C", "G"))
    axis(2)

    col.par <- "gray90"
    abline(h = seq(0, .15, by = .025), col = col.par)
    
    box()
}

plotVals <- function(NsOrNo, CpGorNo, bigAAChangeOrNo, colVal, offset){

    pchpar <- "-"
    cexpar <- 2

    arrowlwd <- 2.5
    arrowlen <-.15
    
    setUpDat <- makeDataFrameToModify(NsOrNo, CpGorNo, bigAAChangeOrNo)[,rownames(modcoef)]

    
    #compute
    pointToPlot <- exp(setUpDat %*% coef.vals)
    arrowToPlotLow <- exp(setUpDat %*% coef.pSE.vals)
    arrowToPlotHigh <- exp(setUpDat %*% coef.mSE.vals)

    #don't plot estimates for CpG sites at Cs or Gs
    restrict <- 1:4
    if(CpGorNo == 1){
        restrict <- 1:2
    }
    
    #plot
    points((1:4)[restrict] + offset, pointToPlot[restrict], col = colVal,  pch = pchpar, cex = cexpar)
    arrows((1:4)[restrict] + offset, arrowToPlotLow[restrict], (1:4)[restrict] + offset, arrowToPlotHigh[restrict], length = arrowlen, code = 3, angle = 90, lwd = arrowlwd, col = colVal)
    
}

## Transparent colors
## modified from Mark Gardener 2015
## www.dataanalytics.org.uk
t.col <- function(color, percent = 50, name = NULL) {
    rgb.val <- col2rgb(color)
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = (100-percent)*255/100)
    return(t.col )
}




plotDat <- function(NsOrNo, CpGorNo, bigAAChangeOrNo, colVal, offsetval){

    pchpar <- 16
    cexpar <- 1
    opacity <- 65

    nsinds <- which(suppDat$TypeOfSite == c("syn", "nonsyn")[NsOrNo + 1])
    cpginds <- which(makesCpG == CpGorNo)
    aachangeinds <- which(bigChange == bigAAChangeOrNo)

    datInds <- intersect(intersect(nsinds, cpginds), aachangeinds)

    xinds <- suppDat$WTnt[datInds]
    for(i in 1:4){
        xinds[xinds == nucord[i]] <- i
    }
    xinds <- as.numeric(xinds)
    yinds <- suppDat$colMeansTs0[datInds]
    xjit <- rnorm(length(xinds), offsetval, .05)
    
    #plot
    points(xinds + xjit, yinds, col = t.col(colVal, percent = opacity),  pch = pchpar, cex = cexpar)
    
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
    arrows(1:4, mus/arrowToPlotLow, 1:4, mus/arrowToPlotHigh, length = .2, code = 3, angle = 90, lwd = 2, col = colVal)
    
}

require(RColorBrewer)
pdf("../out/modeled_freqs_lehman.pdf", width = 12, height = 7)
cols <- brewer.pal(4, "Set2")
layout(matrix(1:2, nrow = 1))
par(mar = c(4, 4.5, 1.5, 1))
makePlot(main = "Synonymous Sites")
plotVals(0, 1, 0, cols[1], .1)
plotVals(0, 0, 0, cols[2], -.1)
plotDat(0, 1, 0, cols[1], .1)
plotDat(0, 0, 0, cols[2], -.1)
abline(v = 1:3 + .5, col = "black")
legend("topleft", c("CpG-forming", "non-CpG-forming"), col = cols[1:2], pch = 16, bg = "white")
makePlot(main = "Non-synonymous Sites")
plotDat(1, 0, 1, cols[3], .1)
plotDat(1, 0, 0, cols[2], -.3)
plotDat(1, 1, 0, cols[1], -.1)
plotDat(1, 1, 1, cols[4], .3)
plotVals(1, 0, 1, cols[3], .1)
plotVals(1, 0, 0, cols[2], -.3)
plotVals(1, 1, 0, cols[1], -.1 )
plotVals(1, 1, 1, cols[4], .3 )
abline(v = 1:3 + .5, col = "black")
legend("topleft", c("CpG-forming", "non-CpG-forming", "Major AA change (non-CpG-forming)",  "Major AA change (CpG-forming)"), col = cols, pch = 16, bg = "white" )
dev.off()


require(RColorBrewer)
pdf("../graphs/modeled_selcoefs.pdf", width = 12, height = 7)
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


#Print xtable for the model
require(xtable)
xtable(sumOfModel, digits = 3)




#length(intersect(intersect(which(makesCpG == 1), which(suppDat$TypeOfSite == "nonsyn")), which(suppDat$WTnt == "a")))
