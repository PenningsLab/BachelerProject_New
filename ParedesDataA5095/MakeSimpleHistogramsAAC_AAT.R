
#Feb 2016, my goal is to make simple SFS plot for talk on Saturday (BAPGXIII)

setwd("~/Dropbox/MarionKristofBachelerProject/GitMarionKristof/bachelerProject/ParedesDataA5095")

Data<-read.csv("MeanFreqDataRanSubCohort.csv")

Freqs103AAC<-Data$K103AACFreq[!is.na(Data$K103AACFreq)]

br<-seq(-0.00001+10^-7,0.1,by=0.00001)

png("Freqs103AAC.png")
hist(Freqs103AAC,xlim=c(-0.00001+10^-7,0.00014),breaks=br, 
     ylab="", xlab = "", main="Single-Site Frequency Spectrum RT103AAC",
     xaxt="n", col="blue"
)
axis(1, at=c(round(br[2],5)-0.000005,round(br[3:16],5)),labels=round(br[2:16],5), col.axis="black", las=2)
hist(rep(0,length(which(Freqs103AAC==0))),breaks=br, add = T, col=1)
mtext("Observed frequency", side=1, line=4, cex.lab=1,las=1, col="black")
mtext("Number of patients", side=2, line=2.5, cex.lab=1,las=3, col="black")
dev.off()

Freqs103AAT<-Data$K103AATFreq[!is.na(Data$K103AATFreq)]

png("Freqs103AAT.png")
hist(Freqs103AAT,xlim=c(-0.00001+10^-7,0.00014),breaks=br, 
     ylab="", xlab = "", main="Single-Site Frequency Spectrum RT103AAT",
     xaxt="n", col="blue"
)
axis(1, at=c(round(br[2],5)-0.000005,round(br[3:16],5)),labels=round(br[2:16],5), col.axis="black", las=2)
hist(rep(0,length(which(Freqs103AAT==0))),breaks=br, add = T, col=1)
mtext("Observed frequency", side=1, line=4, cex.lab=1,las=1, col="black")
mtext("Number of patients", side=2, line=2.5, cex.lab=1,las=3, col="black")
dev.off()
