read.csv("GLM_Analysis/dat/datFitModel.csv",row.names=1)->datFitModel

#Substantially more complicated model with structural elements
#fullmodel.int <- glm(minor ~ t + c + g + bigAAChange + inRT + t*nonsyn + c*nonsyn + g*nonsyn + shape + CpG + CpG*t  + CpG*nonsyn + CpG*nonsyn*t + helix*nonsyn + beta*nonsyn + coil*nonsyn,  family = "binomial", data = datFitModel[datFitModel$res == 0 & datFitModel$stop == 0,])

#Run model 
fullmodel.int <- glm(minor ~ t + c + g + bigAAChange + inRT + t*nonsyn + c*nonsyn + g*nonsyn + shape + CpG + CpG*t  + CpG*nonsyn + CpG*nonsyn*t,  family = "binomial", data = datFitModel[datFitModel$res == 0 & datFitModel$stop == 0,])
sumOfModel <- summary(fullmodel.int)


#fullmodel.small <- glm(minor ~ t + c + g + bigAAChange + t*nonsyn + c*nonsyn + g*nonsyn + CpG + CpG*t  + CpG*nonsyn + CpG*nonsyn*t,  family = "binomial", data = datFitModel[datFitModel$res == 0 & datFitModel$stop == 0,])
#sumOfModel <- summary(fullmodel.small)

modcoef <- sumOfModel$coef

coef.vals <- modcoef[,1]
coef.pSE.vals <- coef.vals + modcoef[,2]
coef.mSE.vals <- coef.vals - modcoef[,2]

modcoef
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
    atcg.mat <- cbind(atcg.mat, atcg.mat[,2:4] * atcg.mat[,7], atcg.mat[,7]*atcg.mat[,12], atcg.mat[, 12] * atcg.mat[,2] , atcg.mat[, 7] * atcg.mat[,2] * atcg.mat[, 12])
    names(atcg.mat) <- c('(Intercept)', names[2:(length(names))], paste(c("t", "c", "g"), ":nonsyn", sep = ""), "nonsyn:CpG", "t:CpG", "t:nonsyn:CpG" )

    return(as.matrix(atcg.mat))
}



makeDataFrameToModify.withSHAPEandinRT <- function(nsOrNo = 0, CpGorNo = 0, bigAAChangeOrNo = 0, shapeval, inRTval){


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
    atcg.mat <- cbind(atcg.mat, atcg.mat[,2:4] * atcg.mat[,7], atcg.mat[,7]*atcg.mat[,12], atcg.mat[, 12] * atcg.mat[,2] , atcg.mat[, 7] * atcg.mat[,2] * atcg.mat[, 12], rep(shapeval, 4), rep(inRTval, 4))
    names(atcg.mat) <- c('(Intercept)', names[2:(length(names))], paste(c("t", "c", "g"), ":nonsyn", sep = ""), "nonsyn:CpG", "t:CpG", "t:nonsyn:CpG",  "shape", "inRT")

    return(as.matrix(atcg.mat))
}

DataFrameOfData <- function(){

    ret.mat <- as.data.frame(matrix(data = 0, ncol = nrow(modcoef), nrow = nrow(suppDatShape)))
    ret.mat[,1] <- 1
    ret.mat[,2] <- as.numeric(suppDat$WTnt == "t")
    ret.mat[,3] <- as.numeric(suppDat$WTnt == "c")
    ret.mat[,4] <- as.numeric(suppDat$WTnt == "g")
    ret.mat[,5] <- bigChange
    ret.mat[,6] <- 0
    ret.mat[1:297,6] <- 1
    ret.mat[,7] <- as.numeric(suppDat$TypeOfSite == "nonsyn")
    ret.mat[,8] <- suppDatShape$SHAPE
    ret.mat[,9] <- makesCpG
    ret.mat[,10] <- ret.mat[,2] * ret.mat[,7]
    ret.mat[,11] <- ret.mat[,3] * ret.mat[,7]
    ret.mat[,12] <- ret.mat[,4] * ret.mat[,7]
    ret.mat[,13] <- ret.mat[,2] * ret.mat[,9]
    ret.mat[,14] <- ret.mat[,7] * ret.mat[,9]
    ret.mat[,15] <- ret.mat[,7] * ret.mat[,9] * ret.mat[,2]
    names(ret.mat) <- rownames(modcoef)
    return(as.matrix(ret.mat))
}


inPRrows <- intersect(which(suppDat$TypeOfSite != "overlap"), 1:297)
inRTrows <- intersect(which(suppDat$TypeOfSite != "overlap"), 298:984)
#mean(exp(DataFrameOfData()[inPRrows,] %*% coef.vals))
#mean(exp(DataFrameOfData()[inRTrows,] %*% coef.vals))

names(suppDat)
mean(suppDat$colMeansTs0[inPRrows])
mean(suppDat$colMeansTs0[inRTrows])

which(suppDat$TypeOfSite == "overlap")

tmpDat.0 <- DataFrameOfData()
tmpDat.0[,8] <- 0
tmpDat.1 <- DataFrameOfData()
tmpDat.1[,8] <- 1

mean(exp(tmpDat.0 %*% coef.vals))
mean(exp(tmpDat.1 %*% coef.vals))

nonsynrows <- which(suppDat$TypeOfSite == "nonsyn")
synrows <- which(suppDat$TypeOfSite == "syn")
mean(exp(DataFrameOfData()[nonsynrows,] %*% coef.vals))
mean(exp(DataFrameOfData()[synrows,] %*% coef.vals))
mean(suppDat$colMeansTs0[nonsynrows])
mean(suppDat$colMeansTs0[synrows])


makePlot <- function(main){

    plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(0.0001, 1), axes = FALSE, ylab = "Predicted frequency of the mutation", xlab = "Mutation type", main = main, log = "y")
    axis(1, at = 1:4, c("A", "T", "C", "G"))
    require(sfsmisc)
    eaxis(2)

    col.par <- "gray95"
#    abline(h = seq(0, .15, by = .025), col = col.par)
    abline(h = (1:9)*10^(-1), col = col.par)
    abline(h = (1:9)*10^(-2), col = col.par)
    abline(h = (1:9)*10^(-3), col = col.par)
    abline(h = (1:9)*10^(-4), col = col.par)
    abline(h = (1:9)*10^(-5), col = col.par)
    
    box()
}


plotVals <- function(NsOrNo, CpGorNo, bigAAChangeOrNo, colVal, offset){

    pchpar <- "-"
    cexpar <- 2

    arrowlwd <- 3.5
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
    xjit <- rnorm(length(xinds), offsetval, .03)

    for(i in 1:4){
        relinds <- which(xinds == i)
#        hist(log(ysToPlot[relinds]))
#        abline(v = mean(log(ysToPlot[relinds])))
        if(length(relinds) > 0){
#            points(i+offsetval, mean(yinds[relinds]), pch = "-", cex = 4)
#            points(i+offsetval, median(yinds[relinds]), pch = "-", cex = 4, col = "blue")
        }
    }
    yinds[yinds == 0] <- 10^(-4.4)
    #plot
    points(xinds + xjit, yinds, col = t.col(colVal, percent = opacity),  pch = pchpar, cex = cexpar)
    
}

#old
makePlot.forS <- function(main){
    require(sfsmisc)
    plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(.00015, .02), axes = FALSE, ylab = "Predicted selection coefficient", xlab = "Mutation type", main = main, log = "y")
    axislabs <- c(expression("A" %->% "G"), expression("T" %->% "C"), expression("C" %->% "T"), expression("G" %->% "A"))
    axis(1, at = 1:4, axislabs)
    eaxis(2, at = 10^c(-2, -3, -4, -5))
    col.par <- "gray95"
    abline(h = (1:9)*10^(-2), col = col.par)
    abline(h = (1:9)*10^(-3), col = col.par)
    abline(h = (1:9)*10^(-4), col = col.par)
    
    box()
}

makePlot.svals <- function(main){
    require(sfsmisc)
    
    plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(.000015, 1), axes = FALSE, ylab = "Predicted selection coefficient", xlab = "Mutation type", main = main, log = "y")

    col.par <- "gray95"
    abline(h = 1, col = col.par)
    abline(h = (1:9)*10^(-1), col = col.par)
    abline(h = (1:9)*10^(-2), col = col.par)
    abline(h = (1:9)*10^(-3), col = col.par)
    abline(h = (1:9)*10^(-4), col = col.par)
    abline(h = (1:9)*10^(-5), col = col.par)
    box()

    axislabs <- c(expression("A" %->% "G"), expression("T" %->% "C"), expression("C" %->% "T"), expression("G" %->% "A"))
    axis(1, at = 1:4, axislabs)
#    axis(1, at = 1:4, c("A", "T", "C", "G"))
    eaxis(2, at = 10^c(0, -1, -2, -3, -4, -5))
    col.par <- "gray95"
}


#Let's do some fresh analysis on this for Marion.
#What's the average shape coefficient?
avShape <- mean(suppDatShape$SHAPE)
avShape
inRTval <- 1
mus <- c(1.11e-05, 1.11e-05, 2.41e-05, 5.48e-05)

# Point 1:
#synonymous CpG forming mutations
mus/exp(makeDataFrameToModify.withSHAPEandinRT(0,1,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

#synonymous non-CpG forming mutations
mus/exp(makeDataFrameToModify.withSHAPEandinRT(0,0,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

#Magnitude changes
(mus/exp(makeDataFrameToModify.withSHAPEandinRT(0,1,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals))/(mus/exp(makeDataFrameToModify.withSHAPEandinRT(0,0,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals))

# Point 2:
mus/exp(makeDataFrameToModify.withSHAPEandinRT(0,0,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

#  Point 3:
#Non-synonymous (changes AA group versus doesn’t change AA group) - among non-CpG forming mutations. 

#Does not change AA group:
mus/exp(makeDataFrameToModify.withSHAPEandinRT(1,0,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

#Does change AA group:
mus/exp(makeDataFrameToModify.withSHAPEandinRT(1,0,1, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

#Magnitude change:
(mus/exp(makeDataFrameToModify.withSHAPEandinRT(1,0,1, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals))/(mus/exp(makeDataFrameToModify.withSHAPEandinRT(1,0,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals))
#selection coefficients
(exp(makeDataFrameToModify.withSHAPEandinRT(1,0,1, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals))
(exp(makeDataFrameToModify.withSHAPEandinRT(1,0,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals))


#  Point 4:
#Nonsynonymous, does not change AA group, does not create new CpG:
noNewCpG <- mus/exp(makeDataFrameToModify.withSHAPEandinRT(1,0,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

#Nonsynonymous, does not change AA group, does create new CpG:
NewCpG <- mus/exp(makeDataFrameToModify.withSHAPEandinRT(1,1,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

#Magnitude change
NewCpG/noNewCpG

#  Point 5
#Non synonymous, non-CpG forming, does not change AA group:
noNewCpG



plotVals.svals <- function(NsOrNo, CpGorNo, bigAAChangeOrNo, colVal, offset){

    pchpar <- 16
    cexpar <- 2
    mus <- c(1.11e-05, 1.11e-05, 2.41e-05, 5.48e-05)
    arrowlen <- .15
    arrowlwd <- 3

    ## NsOrNo <- 1
    ##  CpGorNo <- 1
    ##  bigAAChangeOrNo <- 0
    ##  colVal <- "red"
    ##  offset <- 0

    setUpDat <- makeDataFrameToModify(NsOrNo, CpGorNo, bigAAChangeOrNo)[,rownames(modcoef)]
setUpDat

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
    points((1:4)[restrict] + offset, (mus/pointToPlot)[restrict], col = colVal,  pch = pchpar, cex = cexpar)
    arrows((1:4)[restrict] + offset, (mus/arrowToPlotLow)[restrict], (1:4)[restrict] + offset, (mus/arrowToPlotHigh)[restrict], length = arrowlen, code = 3, angle = 90, lwd = arrowlwd, col = colVal)
  
}

plotDat.svals <- function(NsOrNo, CpGorNo, bigAAChangeOrNo, colVal, offsetval){

    ## NsOrNo <- 1
    ## CpGorNo <- 0
    ## bigAAChange <- 0
    ## offsetval <- -.3
    
    pchpar <- 16
    cexpar <- 1
    opacity <- 65

    mus <- c(1.11e-05, 1.11e-05, 2.41e-05, 5.48e-05)

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
    xjit <- rnorm(length(xinds), offsetval, .03)

    ysToPlot <- mus[xinds]/ yinds
    if(length(yinds == 0) > 1){
        ysToPlot[yinds == 0] <- 1
    }

    for(i in 1:4){
        relinds <- which(xinds == i)
#        hist(log(ysToPlot[relinds]))
#        abline(v = mean(log(ysToPlot[relinds])))
        if(length(relinds) > 0){
#            points(i+offsetval, median(ysToPlot[relinds]), pch = "-", cex = 4)
        }
    }

    #plot
    points(xinds + xjit, ysToPlot, col = t.col(colVal, percent = opacity),  pch = pchpar, cex = cexpar)
    
}


makePlot.axisbreak <- function(main){

    plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(0.00004, 1), axes = FALSE, ylab = "Mutation frequency", xlab = "Mutation type", main = main, log = "y")

    col.par <- "gray95"
#    abline(h = seq(0, .15, by = .025), col = col.par)
    abline(h = (1:9)*10^(-1), col = col.par)
    abline(h = (1:9)*10^(-2), col = col.par)
    abline(h = (1:9)*10^(-3), col = col.par)
    abline(h = (1:9)*10^(-4), col = col.par)
    abline(h = (1:2)*10^(-5), col = col.par)
    abline(h = (1)*10^(-4.4), col = col.par)
    box()

    axislabs <- c(expression("A" %->% "G"), expression("T" %->% "C"), expression("C" %->% "T"), expression("G" %->% "A"))
    axis(1, at = 1:4, axislabs)
    eaxis(2, at = 10^c(-1*(0:4)))
    axis(2, at = 10^(-4.4), label = c("0"), las = 2)
    axis.break(2,2*10^-(4.5),style="slash")
    

}


#geeeeee
library(plotrix)
require(RColorBrewer)
pdf("../out/modeled_freqs_May2017_2.pdf", width = 12, height = 7)
cols <- brewer.pal(4, "Set2")
layout(matrix(1:2, nrow = 1))
par(mar = c(4, 4.5, 1.5, 1))
makePlot.axisbreak(main = "Synonymous Sites")
plotVals(0, 1, 0, cols[1], .1)
plotVals(0, 0, 0, cols[2], -.1)
plotDat(0, 1, 0, cols[1], .1)
plotDat(0, 0, 0, cols[2], -.1)
abline(v = 1:3 + .5, col = "black")
#legend("topleft", c("CpG-forming", "non-CpG-forming"), col = cols[1:2], pch = 16, bg = "white")
legend("bottomright", c("No drastic AA change (non-CpG-forming)", "No drastic AA change (CpG-forming)", "Drastic AA change (non-CpG-forming)",  "Drastic AA change (CpG-forming)"), col = cols[c(2,1,3,4)], pch = 16, bg = "white" )
makePlot.axisbreak(main = "Non-synonymous Sites")
plotDat(1, 0, 1, cols[3], .1)
plotDat(1, 0, 0, cols[2], -.3)
plotDat(1, 1, 0, cols[1], -.1)
plotDat(1, 1, 1, cols[4], .3)
plotVals(1, 0, 1, cols[3], .1)
plotVals(1, 0, 0, cols[2], -.3)
plotVals(1, 1, 0, cols[1], -.1 )
plotVals(1, 1, 1, cols[4], .3 )
abline(v = 1:3 + .5, col = "black")
dev.off()

pdf("../out/modeled_sels_May2017.pdf", width = 12, height = 7)
cols <- brewer.pal(4, "Set2")
layout(matrix(1:2, nrow = 1))
par(mar = c(4, 4.5, 1.5, 1))
makePlot.svals(main = "Synonymous Sites")
plotVals.svals(0, 1, 0, cols[1], .1)
plotVals.svals(0, 0, 0, cols[2], -.1)
plotDat.svals(0, 1, 0, cols[1], .1)
plotDat.svals(0, 0, 0, cols[2], -.1)
abline(v = 1:3 + .5, col = "black")
                                        #legend("topleft", c("CpG-forming", "non-CpG-forming"), col = cols[1:2], pch = 16, bg = "white")
legend("topleft", c("No drastic AA change (non-CpG-forming)", "No drastic AA change (CpG-forming)", "Drastic AA change (non-CpG-forming)",  "Drastic AA change (CpG-forming)"), col = cols[c(2,1,3,4)], pch = 16, bg = "white" )
#legend("topleft", c("Same AA group (non-CpG-forming)", "Same AA group (CpG-forming)", "Changes AA group (non-CpG-forming)",  "Changes AA group (CpG-forming)"), col = cols[c(2,1,3,4)], pch = 16, bg = "white" )
makePlot.svals(main = "Non-synonymous Sites")
plotDat.svals(1, 0, 1, cols[3], .1)
plotDat.svals(1, 0, 0, cols[2], -.3)
plotDat.svals(1, 1, 0, cols[1], -.1)
plotDat.svals(1, 1, 1, cols[4], .3)
plotVals.svals(1, 0, 1, cols[3], .1)
plotVals.svals(1, 0, 0, cols[2], -.3)
plotVals.svals(1, 1, 0, cols[1], -.1 )
plotVals.svals(1, 1, 1, cols[4], .3 )
abline(v = 1:3 + .5, col = "black")
dev.off()



#Print xtable for the model
require(xtable)
xtable(sumOfModel, digits = 3)


#length(intersect(intersect(which(makesCpG == 1), which(suppDat$TypeOfSite == "nonsyn")), which(suppDat$WTnt == "a")))


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


