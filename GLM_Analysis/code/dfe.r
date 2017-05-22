
#Read in the data file and convert the first col to rownames
tmpdat <- read.table("../dat/freqPatTsInclDay0-threshold0.csv", header = TRUE, sep = ",")
dat <- tmpdat[,-1]
rownames(dat) = tmpdat[,1]

#Read in other data file with mutation rates
suppDat <- read.table("../dat/OverviewSelCoeff.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

zanMut <- rep(0, nrow(suppDat))
zanMut[suppDat$WTnt == "a"] <- 5.4*10^(-6)
zanMut[suppDat$WTnt == "t"] <- 8*10^(-6)
zanMut[suppDat$WTnt == "c"] <- 1.1*10^(-5)
zanMut[suppDat$WTnt == "g"] <- 1.5*10^(-5)

suppDat$WTnt
#Given a subset of data, compute the s values
subDat <- dat

runSub <- function(subDat, plotMe = TRUE, zan = FALSE){
#First, filter sites with more than 10% non-zero vals?
#    excludeDueToRefProbs <- which(apply(subDat, 2, function(x){ sum(x > 0, na.rm = TRUE)  })/nrow(subDat) > 0.1)
    excludeDueToRefProbs <- c()
    
#Also, exclude sites that have a lot of NAs (more than 10%)
#    excludeDueToNAs <- which(apply(subDat, 2, function(x){ sum(is.na(x))  })/nrow(subDat) > 0.1) 
    excludeDueToNAs <- c()
    
    allExclusions <- sort(unique(c(excludeDueToRefProbs, excludeDueToNAs)))
    incl <- setdiff(1:ncol(subDat), allExclusions)

    process <- function(cnames){
        toReturn <- c()
        for(i in 1:length(cnames)){
            toReturn[i] <- strsplit(cnames[i], "V")[[1]][2]
        }
        return(as.numeric(toReturn))
    }

    RelInds <- process(colnames(subDat)[incl])
    mus <- suppDat[RelInds,]$TSmutrate
    if(zan == TRUE){ mus <- zanMut[RelInds] }
    
    freqs <- apply(subDat[,incl], 2, mean, na.rm = TRUE)
    ss.inf <- mus/freqs
    ss <- ss.inf
    ss[ss > 1] <- 1

#Let's do a check on our selection coefficents
    cols <- rep("black", length(RelInds))
    cols[suppDat[RelInds,]$TypeOfSite == "nonsyn"] = "red"
    cols[suppDat[RelInds,]$TypeOfSite == "syn"] = "green"
    cols[suppDat[RelInds,]$TypeOfSite == "res"] = "blue"

#I'm reasonably convinced that this worked

#Now, we need to fit a gamma distribution.

#Gamma is not bimodal, so let's exclude the infs for now
    ss.exc <- ss[ss != 1]

    #http://stats.stackexchange.com/questions/160304/fit-gamma-distribution-to-dataset-with-r
    GammaNLL <- function(pars, data){
        alpha <- pars[[1]]
        theta <- pars[[2]]
        return (-sum(dgamma(data, shape = alpha, scale = theta, log = TRUE)))
    }

    require(nloptr)
    Fit <- nloptr(x0 = c(1, 1), eval_f = GammaNLL, lb = c(0,0), data = ss.exc,opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1e15, xtol_rel = 1e-15,  xtol_abs = 1e-15))

    if(plotMe == TRUE){
        layout(matrix(1:2, nrow = 1))
        plot(sort(ss), col = cols[order(ss)])
        hist(ss.exc, freq = FALSE, breaks = 100)
        lines(seq(.002, .5, by = .001), dgamma(seq(.002, .5, by = .001), shape = Fit$solution[1], scale = Fit$solution[2]), col = "red")
    }
    return( list(shape = Fit$solution[1], scale = Fit$solution[2], ests = ss) )
}


#What are all the different conditions we want to test?

yranges <- list()
nrSites <- which(suppDat$TypeOfSite != "res")
#all sites:
yranges[['all']] <- nrSites
#just PR/RT
yranges[['PR']] <- intersect(nrSites, 1:297)
yranges[['RT']] <- intersect(nrSites, 298:ncol(dat))
#just A/T/C/G
yranges[['a']] <- intersect(nrSites, which(suppDat$WTnt == 'a'))
yranges[['t']] <- intersect(nrSites, which(suppDat$WTnt == 't'))
yranges[['c']] <- intersect(nrSites, which(suppDat$WTnt == 'c'))
yranges[['g']] <- intersect(nrSites, which(suppDat$WTnt == 'g'))

allranges <- names(yranges)
toTab <- matrix(data = NA, nrow = 7, ncol = 7)
colnames(toTab) <- c("#sites", "shape", "shape CI", "scale", "scale CI", "p(v del)", "p(v del) CI")
rownames(toTab) <- allranges
#for(j in 1:length(allranges)){
    for(j in 6:7){
    #Grab the relevant columns
    yrange <- yranges[[allranges[j]]]
    subDat <- dat[,yrange]
    #Run and store the relevant estimates
    runit <- runSub(subDat, plotMe = FALSE, zan = FALSE)
    n <- ncol(subDat)
    obs.scale <- runit$scale
    obs.shape <- runit$shape
    obs.p <- sum(runit$ests >= 1)/n
    #bootstrap:
    numreps <- 1000
    btstrp.vals <- matrix(data = NA, nrow = numreps, ncol = 3)
#there are 57 lethal locations
        leths <- c()
    for(i in 1:numreps){
        #sample the data with replacement
        subDat.bs <- dat[, sample(yrange, length(yrange), replace = TRUE)]
        runit.bs <- runSub(subDat.bs, plotMe = FALSE, zan = FALSE)
        leths[i] <- (length(which(runit.bs$ests >= 1)))
        btstrp.vals[i,] <- c(runit.bs$shape, runit.bs$scale, sum(runit.bs$ests >= 1)/n)
    }
    #Create confidence intervals
    shape.int <- quantile(btstrp.vals[,1], c(.025, .975) )
    scale.int <- quantile(btstrp.vals[,2], c(.025, .975) )
    p.int <- quantile(btstrp.vals[,3], c(.025, .975) )
    toTab[j,] <- c(n, round(obs.shape, truncval), format(shape.int), round(obs.scale, truncval), format(scale.int), round(obs.p, truncval), format(p.int))
}


retab <- matrix(data = NA, ncol = 4, nrow = 14)
for(i in 1:7){
    retab[2*i - 1, ] <- toTab[i, c(1, 2, 4, 6)]
    retab[2*i, ] <- c("", toTab[i, c(3, 5, 7)])
}
rows <- rep("", 14)
rows[seq(1, 14, by = 2)] <- allranges
printretab <- cbind(rows, retab)
colnames(printretab) <- c(" ", "Num. sites", "shape", "scale", "Prob(v. deleterious)")
print(xtable(printretab), include.rownames = F)

require(xtable)
xtable(toTab)


bothtab <- cbind(rows, abretab[,1:3], fabretab[,2:4])
colnames(bothtab) <- c(" ", "Num. sites", "shape", "scale","shape", "scale", "Prob(v. deleterious)")
print(xtable(bothtab), include.rownames = F)
                                        #abretab <- retab
#fabretab <- retab
#AbramtoTab <- toTab
#FabiotoTab <- toTab

truncval <- 3
format <- function(twovals){
    return(paste("(",paste(round(twovals,truncval), collapse = ", "), ")", sep = ""))
}
              
set.seed(234819)
params <- list()
for(nuc in c("a", "t", "c", "g")){              

nuc

    params[[nuc]] <- list(shape = NA, scale = NA, p = NA, shaperange = NA, scalerange = NA, prange = NA)
    numreps <- 1000
    nsSites <- which(suppDat$TypeOfSite == "nonsyn")
    yrange <- intersect(which(suppDat$WTnt == nuc), nsSites)
    btstrp.vals <- matrix(data = NA, nrow = numreps, ncol = 3)
    for(i in 1:numreps){
        subDat <- dat[sample(1:nrow(dat), nrow(dat), replace = TRUE), yrange]
        runit <- runSub(subDat, plotMe = FALSE)
        btstrp.vals[i,] <- c(runit$shape, runit$scale, sum(runit$ests >= 1)/length(runit$ests) )
    }
    tvals <- runSub(dat[,yrange], plotMe = FALSE)
    a1vals <- quantile(btstrp.vals[,1], c(.025, .975) )
    a2vals <- quantile(btstrp.vals[,2], c(.025, .975) )
    pvals <- quantile(btstrp.vals[,3], c(.025, .975) )
#direct estimates    
    params[[nuc]]$shape = tvals$shape
    params[[nuc]]$scale = tvals$scale
    params[[nuc]]$p = sum(tvals$ests >= 1)/length(tvals$ests)
#bootstrapped confidence ints
    params[[nuc]]$shaperange = btstrp.vals[,1]
    params[[nuc]]$scalerange = btstrp.vals[,2]
    params[[nuc]]$prange = pvals
   #print
#    print(paste(nuc,": shape: ",round(tvals[1], digits = 3) ," (", paste(round(a1vals, digits = 3), collapse = ", "), ") ", paste("scale: ",round(tvals[2], digits = 3) ," (", paste(round(a2vals, digits = 3), collapse = ", "), ")", sep = ""), sep = ""))
}



              
layout(matrix(1:4, nrow = 2))
par(mar = c(4, 4, 1, 1))
cols[which(nuc == c("a", "t", "c", "g"))]
for(nuc in c("a", "t", "c", "g")){
    plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 20), xlab = "s", ylab = "density", main = paste(toupper(nuc)))
    sc <- params[[nuc]]$scalerange
    sh <- params[[nuc]]$shaperange
    for(i in 1:length(sc)){
        toPlot <- dgamma(ran, shape = sh[i], scale = sc[i])
        toPlot[toPlot > 1000] <- 100000
        lines(ran, toPlot, col = rgb(0,0,0,.01), lwd = 2)
    }
    toPlot <-dgamma(ran, shape = params[[nuc]]$shape, scale = params[[nuc]]$scale)
    toPlot[toPlot > 1000] <- 100000
    lines(ran, toPlot, col = cols[which(nuc == c("a", "t", "c", "g"))], lwd = 2)
 }

              
plot(btstrp.vals[,1], btstrp.vals[,2])
run <- runSub(dat)

hist(log(run$ests) , breaks = 100)
ran <- seq(-10, 0, by = .1)
lines(ran, log(dgamma( exp(ran), shape = run$shape, scale = run$scale)))

require(RColorBrewer)
cols <- brewer.pal(4, "Set1")
layout(matrix(1:4, nrow = 2))
par(mar = c(4, 4, 1, 1))
for (nuc in c("a", "g", "c", "t")){
#    resSites <- which(suppDat$TypeOfSite == "res")
    nsSites <- which(suppDat$TypeOfSite == "nonsyn")
#    yrange <- setdiff(which(suppDat$WTnt == nuc), resSites)
    yrange <- intersect(which(suppDat$WTnt == nuc), nsSites)
    run <- runSub(dat[,yrange], plotMe = FALSE)
    newe <- run$ests[run$ests < 1]
#    hist(log(newe, base = 10) , breaks = seq(-5, 0, by = .2), col= cols[which(nuc == c("a", "t", "c", "g"))], freq = TRUE, xlab = "log10(s)", ylab = "Count", main =  paste(toupper(nuc)))
#    ran <- seq(-5, 0, by = .01)
#    lines(ran, dgamma(10^ran, shape = run$shape, scale = run$scale), lwd = 2, col = "grey")
#    hist(run$ests , breaks = seq(0, 0.05, by = .001), col= cols[which(nuc == c("a", "t", "c", "g"))], main = paste(toupper(nuc)), xlab = "s", ylab = "Count")
#    ran <- seq(0, 1, by = .01)
#    lines(ran, dgamma(ran, shape = run$shape, scale = run$scale), lwd = 2, col = "lightgrey", lty = "solid")
    box()
    print(paste(nuc, length(yrange)))
}


nuc <- "a"
nsSites <- which(suppDat$TypeOfSite == "nonsyn")
yrange <- intersect(which(suppDat$WTnt == nuc), nsSites)
subDat <- dat[,yrange]

              resSites <- which(suppDat$TypeOfSite == "res")
              yrange <- setdiff(which(suppDat$WTnt == nuc), resSites)
              run <- runSub(dat[,yrange], plotMe = FALSE)
              hist(log(run$ests) , breaks = 10, col= cols[which(nuc == c("a", "t", "c", "g"))], ylim = c(0, 200))
              lines(log(sort(unique(run$ests))), dgamma(sort(unique(run$ests)), shape = run$shape, scale = run$scale))

              ran <- seq(-10, 0, by = .1)
              lines(ran, log(dgamma(exp(ran), shape = run$shape, scale = run$scale)), col = "black", lwd = 2)
              lines(log(dgamma(seq(min(run$ests), max(run$ests), by = .001), shape = run$shape, scale = run$scale)), lwd = 2, col = "grey")
              box()
              
dgamma(seq(min(run$ests), max(run$ests), by = .001), shape = run$shape, scale = run$scale)

pdf("../out/nonsyn.atcg.pdf", height = 6, width = 8)
library(plotrix)
layout(matrix(1:4, nrow = 2))
par(bty="n")
gaplen <- c(.65, .95)
cols <- brewer.pal(6, "Set2")[c(1, 2, 3, 6)]
ylabs <- list(a = "Count", t = "Count", c = "", g = "")
xlabs <- list(a = "", t = "s", c = "", g = "s")
mars <- list(a = c(2, 4.5, 2, 1), t = c(4, 4.5, 0, 1), c = c(2, 2, 2, 3.5), g = c(4, 2, 0, 3.5))
for (nuc in c("a", "t", "c", "g")){
    nsSites <- which(suppDat$TypeOfSite == "nonsyn")
    yrange <- intersect(which(suppDat$WTnt == nuc), nsSites)
#    yrange <- which(suppDat$WTnt == nuc)
    run <- runSub(dat[,yrange], plotMe = FALSE)
#    hist(toPlot , breaks = seq(0, .6, by = .01), col= cols[which(nuc == c("a", "t", "c", "g"))])
    forYlim <- hist(run$ests, breaks = seq(0, 1, by = .01), plot = FALSE)
    par(mar = mars[[nuc]])
    gap.plot(0, type= "n", xlim = c(0, 1), ylim = c(0, max(forYlim$counts)), gap = gaplen, gap.axis = "x", xtics = seq(0, 1, by = .2), ytics = seq(0, max(forYlim$counts), by = 10), ylab = ylabs[[nuc]], xlab = xlabs[[nuc]], main = "", cex.lab = 1.5)
                                        #abline(v=seq(gaplen[1], gaplen[2],.001), col="white", lwd = 2)  # hiding vertical lines
    text(.35, max(forYlim$counts), paste(toupper(nuc)), pos = 1, cex = 1.5)
    abline(v = c(.65), col = "white", lwd = 2)
    abline(v = c(.665), col = "white", lwd = 2)
    axis.break(1,.65,style="slash")
    toPlot <- run$ests
    toPlot[toPlot == 1] <- .7
    hist(toPlot, breaks = seq(0, 1, by = .01), add = TRUE, col = cols[which(nuc == c("a", "t", "c", "g"))])
    abline(v = mean(run$ests), lwd = 2, col = "grey", lty = "solid")
}
dev.off()

display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, colorblindFriendly=TRUE)

pdf("../out/syn.atcg.pdf", height = 6, width = 8)
library(plotrix)
layout(matrix(1:4, nrow = 2))
par(bty="n")
gaplen <- c(.65, .95)
cols <- brewer.pal(6, "Set2")[c(1, 2, 3, 6)]
ylabs <- list(a = "Count", t = "Count", c = "", g = "")
xlabs <- list(a = "", t = "s", c = "", g = "s")
mars <- list(a = c(2, 4.5, 2, 1), t = c(4, 4.5, 0, 1), c = c(2, 2, 2, 3.5), g = c(4, 2, 0, 3.5))
for (nuc in c("a", "t", "c", "g")){
    nsSites <- which(suppDat$TypeOfSite == "syn")
    yrange <- intersect(which(suppDat$WTnt == nuc), nsSites)
#    yrange <- which(suppDat$WTnt == nuc)
    run <- runSub(dat[,yrange], plotMe = FALSE)
                                        #    hist(toPlot , breaks = seq(0, .6, by = .01), col= cols[which(nuc == c("a", "t", "c", "g"))])
    breakSto <- seq(0, .05, by = .0008)
    forYlim <- hist(run$ests, breaks = breakSto, plot = FALSE)
    par(mar = mars[[nuc]])
                                        #    gap.plot(0, type= "n", xlim = c(0, 1), ylim = c(0, max(forYlim$counts)), gap = gaplen, gap.axis = "x", xtics = seq(0, 1, by = .2), ytics = seq(0, max(forYlim$counts), by = 10), ylab = ylabs[[nuc]], xlab = xlabs[[nuc]], main = "", cex.lab = 1.5)
    plot(0, type= "n", xlim = c(0, .05), ylim = c(0, max(forYlim$counts)), ylab = ylabs[[nuc]], xlab = xlabs[[nuc]], main = "", cex.lab = 1.5)
                                        #abline(v=seq(gaplen[1], gaplen[2],.001), col="white", lwd = 2)  # hiding vertical lines
    text(.025, max(forYlim$counts), paste(toupper(nuc)), pos = 1, cex = 1.5)
    toPlot <- run$ests
    hist(toPlot, breaks = breakSto, add = TRUE, col = cols[which(nuc == c("a", "t", "c", "g"))])
    abline(v = mean(toPlot), lwd = 2, col = "grey", lty = "solid")
}
dev.off()




box()




layout(matrix(1:4, nrow = 2))
par(mar = c(4, 4, 1, 1))
for (nuc in c("a", "t", "c", "g")){
#    resSites <- which(suppDat$TypeOfSite == "res")
#    yrange <- setdiff(which(suppDat$WTnt == nuc), resSites)
    nsSites <- which(suppDat$TypeOfSite == "nonsyn")
    yrange <- intersect(which(suppDat$WTnt == nuc), nsSites)
    run <- runSub(dat[,yrange], plotMe = FALSE)
    print(paste(nuc, mean(run$ests)))
    toPlot <- run$ests
    hist(toPlot , breaks = seq(0, .6, by = .01), col= cols[which(nuc == c("a", "t", "c", "g"))])
    ran <- seq(-10, 0, by = .1)
    lines(exp(ran), dgamma( exp(ran), shape = run$shape, scale = run$scale), col = "black", lwd = 2)
    box()
}


layout(matrix(1:4, nrow =2))
runSub(dat[,which(suppDat$WTnt == "a")], plotMe = TRUE)
runSub(dat[,which(suppDat$WTnt == "t")], plotMe = TRUE)
runSub(dat[,which(suppDat$WTnt == "c")], plotMe = TRUE)
runSub(dat[,which(suppDat$WTnt == "g")], plotMe = TRUE)
              

layout(matrix(1:16, nrow =4))
runSub(dat[,which(suppDat$WTnt == "t")], plotMe = TRUE)
runSub(dat[sample(1:nrow(dat), nrow(dat), replace = TRUE),which(suppDat$WTnt == "c")], plotMe = TRUE)

            

plot(0, type = "n", xlim = c(0, 1), ylim = c(.0001, 100), log = "y")
lines(seq(0, 1, by = .01), dgamma(seq(0, 1, by = .01), shape = 0.897, scale = 0.01159, log = FALSE), col= "red")
lines(seq(0, 1, by = .01), dgamma(seq(0, 1, by = .01), shape = 0.917, scale = 0.01499, log = FALSE), col= "blue")
lines(seq(0, 1, by = .01), dgamma(seq(0, 1, by = .01), shape = 0.612, scale = 0.09300, log = FALSE), col= "green")
lines(seq(0, 1, by = .01), dgamma(seq(0, 1, by = .01), shape = 0.8177, scale = 0.03902, log = FALSE), col= "yellow")
legend("topright", c("A", "T", "C", "G"), col = c("red", "blue", "green", "yellow"), lty = "solid")
              


