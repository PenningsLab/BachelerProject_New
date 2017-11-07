source("Colors.R")

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
    names(atcg.mat) <- c('(Intercept)', names[2:(length(names))], paste(c("t", "c", "g"), ":nonsyn", sep = ""), "CpG:nonsyn", "t:CpG", "t:CpG:nonsyn" )
    
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
    names(atcg.mat) <- c('(Intercept)', names[2:(length(names))], paste(c("t", "c", "g"), ":nonsyn", sep = ""), "CpG:nonsyn", "t:CpG", "t:CpG:nonsyn",  "shape", "inRT")
    
    return(as.matrix(atcg.mat))
}

DataFrameOfData <- function(){
    
    ret.mat <- as.data.frame(matrix(data = 0, ncol = nrow(modcoef), nrow = nrow(PolShapeData)))
    ret.mat[,1] <- 1
    ret.mat[,2] <- as.numeric(OverviewDF$WTnt == "t")
    ret.mat[,3] <- as.numeric(OverviewDF$WTnt == "c")
    ret.mat[,4] <- as.numeric(OverviewDF$WTnt == "g")
    ret.mat[,5] <- bigChange
    ret.mat[,6] <- 0
    ret.mat[1:297,6] <- 1
    ret.mat[,7] <- as.numeric(OverviewDF$TypeOfSite == "nonsyn")
    ret.mat[,8] <- PolShapeData$SHAPE
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
    
    pchpar <- 16
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

makePlot.svals <- function(main){
    require(sfsmisc)
    
    plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(.00015, 1), axes = FALSE, ylab = "Predicted selection coefficient", xlab = "Mutation type", main = main, log = "y")
    
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

plotVals.svals <- function(NsOrNo, CpGorNo, bigAAChangeOrNo, colVal, offset, mutrates){
    
    pchpar <- 16
    cexpar <- 2
    #mus <- c(1.11e-05, 1.11e-05, 2.41e-05, 5.48e-05)
    mus <- mutrates
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

plotDat.svals <- function(NsOrNo, CpGorNo, bigAAChangeOrNo, colVal, offsetval,mutrates){
    
    ## NsOrNo <- 1
    ## CpGorNo <- 0
    ## bigAAChange <- 0
    ## offsetval <- -.3
    
    pchpar <- 16
    cexpar <- 1
    opacity <- 65
    
    #mus <- c(1.11e-05, 1.11e-05, 2.41e-05, 5.48e-05)
    mus <- mutrates

    nsinds <- which(OverviewDF$TypeOfSite == c("syn", "nonsyn")[NsOrNo + 1])
    cpginds <- which(makesCpG == CpGorNo)
    aachangeinds <- which(bigChange == bigAAChangeOrNo)
    
    datInds <- intersect(intersect(nsinds, cpginds), aachangeinds)
    
    xinds <- as.character(OverviewDF$WTnt[datInds])
    for(i in 1:4){
        xinds[xinds == nucord[i]] <- i
    }
    xinds <- as.numeric(xinds)
    yinds <- OverviewDF$MeanFreq[datInds]
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
    
    nsinds <- which(OverviewDF$TypeOfSite == c("syn", "nonsyn")[NsOrNo + 1])
    cpginds <- which(makesCpG == CpGorNo)
    aachangeinds <- which(bigChange == bigAAChangeOrNo)
    
    datInds <- intersect(intersect(nsinds, cpginds), aachangeinds)
    
    xinds <- as.character(OverviewDF$WTnt[datInds])
    for(i in 1:4){
        xinds[xinds == nucord[i]] <- i
    }
    xinds <- as.numeric(xinds)
    yinds <- OverviewDF$MeanFreq[datInds]
    xjit <- rnorm(length(xinds), offsetval, .03) #jitter
    
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
if (FALSE){ #Pleuni: since Alison marked this as old, I think we don't need it. 
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
}

