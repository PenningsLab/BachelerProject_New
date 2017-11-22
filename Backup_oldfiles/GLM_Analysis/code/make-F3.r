#actually, currently figure 2

Bach.dat <- read.table("../dat/OverviewSelCoeffwSHAPE.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
Lehman.dat <- read.table("../dat/OverviewSelCoeffLehman.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
Zanini.dat <- read.table("../dat/OverviewSelCoeffZanini-v2.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)


library(plotrix)

plotter <- function(datset, typematch){


    splitPlot = 0
    if(datset == "Lehman"){
        main.dat <- Lehman.dat[Lehman.dat$TypeOfSite != "res",]
        if(typematch == "nonsyn"){
            xlim.set = c(0, 1)
        }else{
            xlim.set = c(0, 1)
        }
    }
    if(datset == "Zanini"){
        main.dat <- Zanini.dat[Zanini.dat$TypeOfSite != "exclude" & Zanini.dat$TypeOfSite != "res",]
        if(typematch == "nonsyn"){
            xlim.set = c(0, 0.65)
        }else{
            xlim.set = c(0, 0.16)
        }
    }
    if(datset == "Bacheler"){
        main.dat <- Bach.dat[Bach.dat$TypeOfSite != "res",]
        if(typematch == "nonsyn"){
            xlim.set = c(0, 1)
            splitPlot = 1 
        }else{
            xlim.set = c(0, .05)
        }
    }

    print(paste(datset, "xlims = ", paste(xlim.set, collapse = " - "), splitPlot, collapse = " "))
#    layout(matrix(1:4, nrow = 2))
    par(mfcol=c(2,2),oma = c(0, 0, 3, 0))
    par(bty="n")
    gaplen <- c(.65, .95)

    cols <- brewer.pal(6, "Set2")[c(1, 2, 3, 6)]
    ylabs <- list(a = "Count", t = "Count", c = "", g = "")
    xlabs <- list(a = "", t = "s", c = "", g = "s")
    mars <- list(a = c(2, 4.5, 2, 1), t = c(4, 4.5, 0, 1), c = c(2, 2, 2, 3.5), g = c(4, 2, 0, 3.5))

    axislabs <- c(expression("A" %->% "G"), expression("T" %->% "C"), expression("C" %->% "T"), expression("G" %->% "A"))
    for (nuc in c("a", "t", "c", "g")){

        siteType <- which(main.dat$TypeOfSite == typematch)
        yrange <- intersect(which(main.dat$WTnt == nuc), siteType)
        sToPlot <- main.dat[yrange, ]$EstSelCoeff
        meanS <- mean(sToPlot)
        plotBreaks <- seq(0, xlim.set[2], by = xlim.set[2]/30)
        forYlim <- hist(sToPlot, breaks = plotBreaks, plot = FALSE)
        par(mar = mars[[nuc]])

        cex.lab.val <- 2
        if(splitPlot == 1){
            gap.plot(0, type= "n", xlim = c(0,1), ylim = c(0, max(forYlim$counts)), gap = gaplen,cex.axis = 5, gap.axis = "x", xtics = seq(0, 1, by = .2), ytics = seq(0, max(forYlim$counts), by = 10), ylab = ylabs[[nuc]], xlab = xlabs[[nuc]], main = "", cex.lab = cex.lab.val)
            abline(v = c(.65), col = "white", lwd = 2)
            abline(v = c(.665), col = "white", lwd = 2)
            axis.break(1,.65,style="slash")
            text(.35, max(forYlim$counts), axislabs[which(nuc == c("a", "t", "c", "g"))], pos = 1, cex = cex.lab.val)
            sToPlot[sToPlot == 1] <- .7
            hist(sToPlot, breaks = plotBreaks, add = TRUE, col = cols[which(nuc == c("a", "t", "c", "g"))])
        }else{
            hist(sToPlot, breaks = plotBreaks, add = FALSE, col = cols[which(nuc == c("a", "t", "c", "g"))], main = "", xlim = xlim.set, ylim = c(0, max(forYlim$counts)), ylab = ylabs[[nuc]], xlab = xlabs[[nuc]], cex.lab =  cex.lab.val)
#            text(mean(xlim.set), max(forYlim$counts), paste(toupper(nuc)), pos = 1, cex = cex.lab.val)
            text(mean(xlim.set), max(forYlim$counts), axislabs[which(nuc == c("a", "t", "c", "g"))], pos = 1, cex = cex.lab.val)
        }
        abline(v = meanS, lwd = 2, col = "grey", lty = "solid")
    }
    if(typematch == "syn"){
        typematch.out = "Synonymous"
    }else if(typematch == "nonsyn"){
        typematch.out = "Non-synonymous"
    }
    mtext(paste(typematch.out , " sites", sep = ""), outer = TRUE, cex = 1.5)
#        mtext(paste(datset, " data, ",typematch.out , " sites", sep = ""), outer = TRUE, cex = 1.5)
}



for(dat.file in c("Lehman", "Zanini", "Bacheler")){
    for(sitetype in c("nonsyn", "syn")){
        pdf(paste("../out/F2-", dat.file, "-", sitetype, ".pdf", sep = ""), height = 8, width = 8)
        plotter(dat.file, sitetype)
        dev.off()
    }
}


#Enter the Zanini mutation rates
zanMut <- rep(0, nrow(suppDat))
zanMut[suppDat$WTnt == "a"] <- 5.4*10^(-6)
zanMut[suppDat$WTnt == "t"] <- 8*10^(-6)
zanMut[suppDat$WTnt == "c"] <- 1.1*10^(-5)
zanMut[suppDat$WTnt == "g"] <- 1.5*10^(-5)
