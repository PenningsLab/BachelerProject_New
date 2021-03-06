#Read in the data file and convert the first col to rownames
source('Rscripts/baseRscript.R')

OverviewDFBacheler <- read.table("Output/OverviewSelCoeff_BachelerFilter.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFLehman <- read.table("Output/OverviewSelCoeffLehman.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFZanini <- read.table("Output/OverviewSelCoeffZanini.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

plotter <- function(datset){
    par(mar = c(4.5, 4.5, 2, 2))
    if(datset == "Lehman"){
        main.dat <- OverviewDFLehman
        wheresthebreak <- 6
        remap = 1
    }
    if(datset == "Zanini"){
        main.dat <- OverviewDFZanini
        wheresthebreak <- 5
        remap = 0
    }
    if(datset == "Bacheler"){
        main.dat <- OverviewDFBacheler
        wheresthebreak <- 5
        remap = 1
    }
    dataset <- main.dat[intersect(intersect(which(main.dat$TypeOfSite != "res"), which(main.dat$TypeOfSite != "overlap")), which(main.dat$TypeOfSite != "exclude") ),]
    if(datset == "Zanini"){ toPlot <- dataset$colMeansTsZanini }
    if(datset == "Bacheler"){ toPlot <- dataset$MeanFreq } #MeanFreq replaces colMeansTs0
    if(datset == "Lehman"){ toPlot <- dataset$colMeansTsLehman }
    toPlot <- toPlot[!is.na(toPlot)]
    colVect <- rep(0, nrow(dataset))
    colVect[dataset$TypeOfSite == "nonsyn"] <- cols[5]
    colVect[dataset$TypeOfSite == "syn"] <- cols[3]
    colVect[dataset$TypeOfSite == "stop"] <- "black"
    plot(5, type = "n", log = "y", axes = FALSE, xlim = c(0, length(toPlot[!is.na(toPlot)])), ylim = c(10^-(wheresthebreak), max(toPlot, na.rm = TRUE)),  ylab = "Mean mutation frequency", xlab = "Mutations ordered by mean mutation frequency", cex.lab = 1.3)
    for(i in 1:wheresthebreak){
        abline(h = 1:10 * 10^(-i), col = "gray70")
    }
    abline(h = 10^-(wheresthebreak), col = "gray70")
    if(remap == 1){
        eaxis(side = 2, at = 10^((-1):(-(wheresthebreak-1))))
        axis(side = 2, at = c(1, 10^-(wheresthebreak)), label = c(1, 0), las = 2)
        box()
        axis.break(2,2*10^-(wheresthebreak),style="slash")
    }else{
        eaxis(side = 2, at = 10^((0):(-(wheresthebreak))))
        axis(side = 2, at = 1, label = 1, las =2)
        box()
    }
    cexval <- 1.5
    toPlot[toPlot == 0] <- 10^-(wheresthebreak)
    points(1:length(toPlot), sort(toPlot), col = colVect[order(toPlot)], pch = "|", cex = cexval)
    axis(1)
    legend("bottomright", c("Synonymous", "Non-synonymous", "Nonsense"), col = c(cols[3], cols[5], "black"), pch = "|", bg = "white", pt.cex = cexval)
}

for(dat.file in c("Lehman", "Zanini", "Bacheler")){
    png(paste("Output/F1-ordered-Nov2017", dat.file, "-v3.png", sep = ""), height = 6, width = 9,units="in",res=100)
    plotter(dat.file)
    dev.off()
}

