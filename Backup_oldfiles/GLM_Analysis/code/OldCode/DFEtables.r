setwd("/Users/pleuni/Documents/Git/bachelerProject")
source("Rscripts/baseRscript.R")

OverviewDFBacheler <- read.table("Output/OverviewSelCoeff_BachelerFilter.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFLehman <- read.table("Output/OverviewSelCoeffLehman.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFZanini <- read.table("Output/OverviewSelCoeffZanini.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

#Read in the data file and convert the first col to rownames
#Bach.dat <- read.table("../dat/OverviewSelCoeffwSHAPE.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
#Bach.dat <- read.table("../dat/OverviewSelCoeffwProteinFeatures.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
#Lehman.dat <- read.table("../dat/OverviewSelCoeffLehman-2.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
#Zanini.dat <- read.table("../dat/OverviewSelCoeffZanini-v2.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

list.files("../dat/")
#Read in other data file with mutation rates
# PSP: next line copied from dfe.r
suppDat <- read.table("../dat/OverviewSelCoeff.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
suppDat <-OverviewDFBacheler

zanMut <- rep(0, nrow(suppDat))
zanMut[suppDat$WTnt == "a"] <- 5.4*10^(-6)
zanMut[suppDat$WTnt == "t"] <- 8*10^(-6)
zanMut[suppDat$WTnt == "c"] <- 1.1*10^(-5)
zanMut[suppDat$WTnt == "g"] <- 1.5*10^(-5)

runSub <- function(datName, subDat, plotMe = TRUE, zan = FALSE){
    if(datName == "bach"){
        freqColName <- "MeanFreq"
    }
    if(datName == "zanini"){
        freqColName <- "colMeansTsZanini"
    }
    if(datName == "lehman"){
        freqColName <- "colMeansTsLehman"
    }

    process <- function(cnames){
        toReturn <- c()
        for(i in 1:length(cnames)){
            toReturn[i] <- strsplit(cnames[i], "V")[[1]][2]
        }
        return(as.numeric(toReturn))
    }

    mus <- subDat$TSmutrate
    if(zan == TRUE){
        zanMut <- rep(0, nrow(subDat))
        zanMut[subDat$WTnt == "a"] <- 5.4*10^(-6)
        zanMut[subDat$WTnt == "t"] <- 8*10^(-6)
        zanMut[subDat$WTnt == "c"] <- 1.1*10^(-5)
        zanMut[subDat$WTnt == "g"] <- 1.5*10^(-5)
        mus <- zanMut
    }

    freqs <- subDat[,freqColName]

    ss.inf <- mus/freqs
    ss <- ss.inf
    ss[ss > 1] <- 1

#Now, we need to fit a gamma distribution.
#Gamma is not bimodal, so let's exclude the infs for now
    ss.exc <- ss[ss != 1]
    ss.exc <- ss[!is.na(ss)]

    #http://stats.stackexchange.com/questions/160304/fit-gamma-distribution-to-dataset-with-r
    GammaNLL <- function(pars, data){
        alpha <- pars[[1]]
        theta <- pars[[2]]
        return (-sum(dgamma(data, shape = alpha, scale = theta, log = TRUE)))
    }


    require(nloptr)
    Fit <- nloptr(x0 = c(1, 1), eval_f = GammaNLL, lb = c(0,0), data = ss.exc,opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1e15, xtol_rel = 1e-15,  xtol_abs = 1e-15))

    Fit <- nloptr(x0 = c(1, 1), eval_f = GammaNLL, lb = c(0,0), data = ss.exc,opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1e15, xtol_rel = 1e-10,  xtol_abs = 1e-10))

    return( list(shape = Fit$solution[1], scale = Fit$solution[2], ests = ss) )
}


layout(matrix(1:2, nrow = 1))
testseq <- seq(-.2, .5, by = .001)
plotp.1 <- c()
for(i in 1:length(testseq)){
    plotp.1[i] <- sum(dgamma(ss.exc, shape = Fit$solution[1] + testseq[i], scale = Fit$solution[2], log = TRUE))
}
plot(testseq,plotp.1)
abline(v = testseq[which(plotp.1 == max(plotp.1))])
testseq <- seq(-.2, .5, by = .001)
plotp.2 <- c()
for(i in 1:length(testseq)){
    plotp.2[i] <- sum(dgamma(ss.exc, shape = Fit$solution[1], scale = Fit$solution[2] + testseq[i], log = TRUE))
}
plot(testseq,plotp.2)
abline(v = testseq[which(plotp.2 == max(plotp.2))])

auths <- c("bach", "zanini", "lehman")
tOrF <- c(TRUE, FALSE)
toXtable <- matrix(data = "", nrow = 6, ncol  = 6)
for(i in auths){
    print(i)
    rowForObs <- 2 * which(i == auths) - 1
    if(i == "bach"){
        subDat <- OverviewDFBacheler
    }
    if(i == "zanini"){
        subDat <- OverviewDFZanini
    }
    if(i == "lehman"){
        subDat <- OverviewDFLehman
    }
(table(subDat[grep("syn|stop", subDat$TypeOfSite),]$TypeOfSite))
#(table(subDat[grep("syn|stop", subDat$TypeOfSite, invert = TRUE),]$TypeOfSite))
    for(j in c(1:2)){
#        yrange <- 1:nrow(subDat)
#        if(i == "lehman"){
#            yrange <- which(!is.na(subDat$colMeansTsLehman))
        yrange <- grep("syn|stop", subDat$TypeOfSite)
#        }
        runit <- runSub(i, subDat[yrange,], zan = tOrF[j])
        obs.p <- sum(runit$ests >= 1, na.rm = TRUE)/length(runit$ests)
        obs.scale <- runit$scale
        obs.shape <- runit$shape
        toXtable[rowForObs, 1] <- nrow(subDat[yrange,])
        toXtable[rowForObs, 2*j] <- round(obs.scale, 3)
        toXtable[rowForObs, 2*j + 1] <- round(obs.shape, 3)
        toXtable[rowForObs, 6] <- round(obs.p, 3)
        numreps <- 1000
        btstrp.vals <- matrix(data = NA, nrow = numreps, ncol = 3)
        leths <- c()
        for(k in 1:numreps){
            subDat.bs <- subDat[sample(yrange, length(yrange), replace = TRUE),]
            runit.bs <- runSub(i, subDat.bs, plotMe = FALSE, zan = tOrF[j])
            leths[k] <- (length(which(runit.bs$ests >= 1)))
            btstrp.vals[k,] <- c(runit.bs$shape, runit.bs$scale, sum(runit.bs$ests >= 1)/length(yrange))
        }
    #Create confidence intervals
        shape.int <- quantile(btstrp.vals[,1], c(.025, .975) )
        scale.int <- quantile(btstrp.vals[,2], c(.025, .975) )
is.na(btstrp.vals[,3])
        if(sum(is.na(btstrp.vals[,3]))== numreps){ p.int <- c(0,0)
                                               }else{
            p.int <- quantile(btstrp.vals[,3], c(.025, .975) )
        }
        toXtable[rowForObs + 1, 2*j] <- paste("(",paste(round(scale.int, digits = 3), collapse = ", "), ")", sep = "")
        toXtable[rowForObs + 1, 2*j + 1] <- paste("(",paste(round(shape.int, digits = 3), collapse = ", "), ")", sep = "")
        toXtable[rowForObs + 1, 6] <- paste("(",paste(round(p.int, digits = 3), collapse = ", "), ")", sep = "")
    }
}
warnings()
colnames(toXtable) <- c("Sites",  "Scale", "Shape", "Scale", "Shape", "Lethal")
rownames(toXtable) <- c("Bacheler",  " ", "Zanini", "  ", "Lehman", "   ")
require(xtable)
xtable(toXtable)

if (FALSE){

#pdf("../out/fitDFEsabr.pdf", width = 5, height = 5)
par(mar = c(4.5, 4.5, 1, 1))
scales <- as.numeric(toXtable[c(1,3, 5),2])
shapes <- as.numeric(toXtable[c(1,3, 5),3])
plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 20), ylab = "Density", xlab = "Estimated s", cex.lab = 1.5)
plotseq <- seq(0.000001, 1, by = .01)
for(i in 1:3){
    lines(plotseq, dgamma(plotseq, scale =  scales[i], shape = shapes[i]), col = cols[i], lwd = 2)
}
legend("topright", c("Bacheler", "Zanini", "Lehman"), col = cols[1:3], lwd = 2)
abline(v = 0, lty = "dashed", lwd = 2, col = "gray90")
#dev.off()


#pdf("../out/fitDFEzan.pdf", width = 5, height = 5)
par(mar = c(4.5, 4.5, 1, 1))
scales <- as.numeric(toXtable[c(1,3, 5),4])
shapes <- as.numeric(toXtable[c(1,3, 5),5])
plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 20), ylab = "Density", xlab = "Estimated s", cex.lab = 1.5)
plotseq <- seq(0.000001, 1, by = .01)
for(i in 1:3){
    lines(plotseq, dgamma(plotseq, scale =  scales[i], shape = shapes[i]), col = cols[i], lwd = 2)
}
legend("topright", c("Bacheler", "Zanini", "Lehman"), col = cols[1:3], lwd = 2)
abline(v = 0, lty = "dashed", lwd = 2, col = "gray90")
#dev.off()

}