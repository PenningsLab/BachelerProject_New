#This script makes / does
# SingleSiteFrequencySpectraPRO.pdf
# SingleSiteFrequencySpectraRT.pdf
# ProteaseRanking.pdf
# ProteaseRanking_EffectFiltering.pdf (PSP 2017 June Not working currently)
# RTRanking.pdf
# RTRanking_EffectFiltering.pdf  (PSP 2017 June Not working currently)

setwd("~/Documents/Git/bachelerProject/Rscripts")
source('./baseRscript.R')
library(scales)
library(plotrix)
library(RColorBrewer)

#PSP Aug 2017 changed the names of the files read so that it includes the threshhold / filter
read.table("../Output/freqPatTs_Bacheler_Threshold1.csv",sep=",",header=TRUE,row.names=1)->freqPatTs0
read.csv("../Output/OverviewSelCoeff_BachelerFilter.csv")->OverviewDF

OverviewDFOrderedByFreq <- OverviewDF[order(OverviewDF$MeanFreq),] 
PROdata<-OverviewDFOrderedByFreq[OverviewDFOrderedByFreq$num<298,]

#Make a figure with all single site frequency spectra for Protease (and then RT)
if (FALSE){
pdf("../Output/SingleSiteFrequencySpectraPRO.pdf",width=8,height=10)
par(mfrow=c(3,3))
for (i in 1:297){
    #Site excluded bc too many non-consensus patients? Threshold = 10%
    #PSP June 2017 I removed this exclusion criterion. Code doesn't work. NumPats66Excluded is missing.
    #    excluded = ""; if (NumPats66Excluded[i]>0.1) excluded = "(site excluded)"
    resistancesite = ""; if (OverviewDFOrderedByFreq$TypeOfSite[which(OverviewDFOrderedByFreq$num==i)]=="res") resistancesite = "\n resistance codon, excluded"
    #first create empty plot with title
    hist(rep(0,70),breaks=seq(0,1,by=0.02),xlim=c(0,1),ylim=c(0,70),yaxt="n",col=0,border=0,main = paste("Protease site", i,"\n","AA",ceiling(i/3),resistancesite),xlab="Frequency",ylab="Count")
    #Next, show true height of 0 bar
    if (length(which(freqPatTs0[,i]<0.02))>=70){
        hist(rep(0,70),breaks=seq(0,1,by=0.02),xlim=c(0,1),ylim=c(0,70),yaxt="n",col=OverviewDFOrderedByFreq$color[which(OverviewDFOrderedByFreq$num==i)],add=T)}
    #next show all data (unfiltered), but only until 50 for 0 cat
    hist(c(rep(0,min(60,length(which(freqPatTs0[,i]<0.02)))),freqPatTs0[,i][which(freqPatTs0[,i]>0)]),breaks=seq(0,1,by=0.02),add=T,
         col=OverviewDFOrderedByFreq$color[which(OverviewDFOrderedByFreq$num==i)])
    axis(2,labels = c(10,20,30,40,50,max(70,length(which(freqPatTs0[,i]<0.02)))), at = c(10,20,30,40,50,70), las=1)
    if (length(which(freqPatTs0[,i]<0.02))>=70){
        axis.break(axis=2,breakpos=60,bgcol="white",breakcol="black",style="slash",brw=0.02)
        points(c(0.01,0.02),c(60,60),pch=15,cex=2.5,col="white")
    }else{axis(2,labels = 60,at=60,las=1)}
}
dev.off()}

#Make a figure with all single site frequency spectra for RT
if (FALSE){
pdf("../Output/SingleSiteFrequencySpectraRT.pdf",width=8,height=10)
par(mfrow=c(3,3))
for (i in 298:length(OverviewDFOrderedByFreq$color)){
    #Site excluded bc too many non-consensus patients? Threshold = 10%
    #excluded = ""; if (NumPats66Excluded[i]>0.1) excluded = "(site excluded)"
    resistancesite = ""; if (OverviewDFOrderedByFreq$TypeOfSite[which(OverviewDFOrderedByFreq$num==i)]=="res") resistancesite = "\n resistance codon, excluded"
    #first create empty plot with title
    hist(rep(0,70),breaks=seq(0,1,by=0.02),xlim=c(0,1),ylim=c(0,70),yaxt="n",col=0,border=0,main = paste("RT site", i-297,"\n","AA",ceiling((i-297)/3),resistancesite),xlab="Frequency",ylab="Count")
    #Next, show true height of 0 bar
    if (length(which(freqPatTs0[,i]<0.02))>=70){
        hist(rep(0,70),breaks=seq(0,1,by=0.02),xlim=c(0,1),ylim=c(0,70),yaxt="n",col=OverviewDFOrderedByFreq$color[which(OverviewDFOrderedByFreq$num==i)],add=T)}
    #next show all data (unfiltered), but only until 50 for 0 cat
    hist(c(rep(0,min(60,length(which(freqPatTs0[,i]<0.02)))),freqPatTs0[,i][which(freqPatTs0[,i]>0)]),breaks=seq(0,1,by=0.02),add=T,
         col=OverviewDFOrderedByFreq$color[which(OverviewDFOrderedByFreq$num==i)])
    axis(2,labels = c(10,20,30,40,50,max(70,length(which(freqPatTs0[,i]<0.02)))), at = c(10,20,30,40,50,70), las=1)
    if (length(which(freqPatTs0[,i]<0.02))>=70){
        axis.break(axis=2,breakpos=60,bgcol="white",breakcol="black",style="slash",brw=0.02)
        points(c(0.01,0.02),c(60,60),pch=15,cex=2.5,col="white")
    }else{axis(2,labels = 60,at=60,las=1)}
}
dev.off()}

#All of the rest of this file is now replaced by "ranking_Ordered_Figure1.R"
#PRO: Make the plots (transitions) for ranking
if (FALSE){
pdf("../Output/ProteaseRanking_Aug2017.pdf",width = 13, height = 10)
par(mfrow=c(1,1))
PROdata<-OverviewDFOrderedByFreq[OverviewDFOrderedByFreq$num<298 & OverviewDFOrderedByFreq$TypeOfSite %in%c("nonsyn","stop","syn") ,]
#remove resistance mutations
#PROdata<-PROdata[PROdata$TypeOfSite!="res",]
#remove positions with too many non-consensus patients
#which(NumPats66Excluded>0.1)->listSitesToExclude
#PROdata<-PROdata[-which(PROdata$num %in% listSitesToExclude),]
plot(log(PROdata$MeanFreq+0.001), main = "Protease mutant frequencies",
     ylim=c(log(0.001),log(0.5)),cex=1.5, pch = 16, col=alpha(PROdata$color, 1), xlab = "Nucleotides ordered by mean mutation frequency", ylab = "Mean mutation frequency" , yaxt = "n")
axis(2,labels = c(0,0.001, 0.005, 0.05, 0.1), at = log(c(0.001, 0.002, 0.006, 0.051, 0.101)),las=1)
points(10*PROdata$NumPats66Excluded+log(0.001),pch=16,cex=0.5,col="grey")
#axis(4,labels = c(0,0.05,0.1), at = 10*c(0,0.05,0.1)+log(0.001),las=1)
dev.off()

#make figure with ranked / ordered frequencies with all data

plot(5, type = "n", log = "y", axes = FALSE, xlim = c(0, length(toPlot[!is.na(toPlot)])), 
     ylim = c(10^-(wheresthebreak), max(toPlot, na.rm = TRUE)),  
     ylab = "Mean mutation frequency", xlab = "Mutations ordered by mean mutation frequency", 
     cex.lab = 1.5)


#pdf("../Output/PolRanking_Aug2017.pdf",width = 13, height = 10)
par(mfrow=c(1,1))
POLdata<-OverviewDFOrderedByFreq[OverviewDFOrderedByFreq$TypeOfSite %in%c("nonsyn","stop","syn") ,]
plot(log(POLdata$MeanFreq+0.001), main = "Protease mutant frequencies",
     ylim=c(log(0.001),log(0.5)),cex=1.5, pch = 16, col=alpha(POLdata$color, 1), 
#     xlab = "Nucleotides ordered by mean mutation frequency", ylab = "Mean mutation frequency" ,
     ylab = "Mean mutation frequency", xlab = "Mutations ordered by mean mutation frequency", 
     yaxt = "n")
axis(2,labels = c(0,0.001, 0.005, 0.05, 0.1), at = log(c(0.001, 0.002, 0.006, 0.051, 0.101)),las=1)
#dev.off()



#Show the effect of our filtering out non-consensus patients at day 0
#PSP: continue here! 

pdf("../Output/ProteaseRanking_EffectFiltering.pdf",width = 13, height = 10)
plot(log(PROdata$MeanFreq+0.001), main = "Protease mutant frequencies",
     ylim=c(log(0.001),log(0.5)),cex=1.5, pch = 16, col=alpha(PROdata$color, 1), xlab = "Nucleotides ordered by mean mutation frequency", ylab = "Mean mutation frequency" , yaxt = "n")
axis(2,labels = c(0,0.001, 0.005, 0.05, 0.1), at = log(c(0.001, 0.002, 0.006, 0.051, 0.101)),las=1)
#points(10*PROdata$NumPats66Excluded+log(0.001),pch=16,cex=0.5,col="grey")
#axis(4,labels = c(0,0.05,0.1), at = 10*c(0,0.05,0.1)+log(0.001),las=1)
points(log(PROdata$colMeansTs66[order(PROdata$colMeansTs66)]+0.001),cex=1.5, pch = 16, col=alpha(PROdata$color, 0.3))
dev.off()


#RT: Make the plots (transitions) for rankng

#pdf("../Output/RTRanking.pdf",width = 13, height = 10)
RTdata<-OverviewDFOrderedByFreq[OverviewDFOrderedByFreq$num>=298,]
RTdata<-RTdata[RTdata$TypeOfSite!="res",]
#remove positions with too many non-consensus patients
which(NumPats66Excluded>0.1)->listSitesToExclude
RTdata<RTdata[-which(RTdata$num %in% listSitesToExclude),]
plot(log(RTdata$MeanFreq+0.001), main = "RT mutant frequencies",
     ylim=c(log(0.001),log(0.5)),cex=1.5, pch = 16, col=alpha(RTdata$color, 1), xlab = "Nucleotides ordered by mean mutation 
     frequency", ylab = "Mean mutation frequency" , yaxt = "n")
axis(2,labels = c(0,0.001, 0.005, 0.05, 0.1), at = log(c(0.001, 0.002, 0.006, 0.051, 0.101)),las=1)
points(10*PROdata$NumPats66Excluded+log(0.001),pch=16,cex=0.5,col="grey")
axis(4,labels = c(0,0.05,0.1), at = 10*c(0,0.05,0.1)+log(0.001),las=1)
#dev.off()


#Show the effect of our filtering out non-consensus patients at day 0

pdf("../Output/RTRanking_EffectFiltering.pdf",width = 13, height = 10)
plot(log(RTdata$MeanFreq+0.001), main = "RT mutant frequencies",
     ylim=c(log(0.001),log(0.5)),cex=1.5, pch = 16, col=alpha(PROdata$color, 1), xlab = "Nucleotides ordered by mean mutation frequency", ylab = "Mean mutation frequency" , yaxt = "n")
axis(2,labels = c(0,0.001, 0.005, 0.05, 0.1), at = log(c(0.001, 0.002, 0.006, 0.051, 0.101)),las=1)
#points(10*PROdata$NumPats66Excluded+log(0.001),pch=16,cex=0.5,col="grey")
#axis(4,labels = c(0,0.05,0.1), at = 10*c(0,0.05,0.1)+log(0.001),las=1)
points(log(RTdata$colMeansTs66[order(RTdata$colMeansTs66)]+0.001),cex=1.5, pch = 16, col=alpha(RTdata$color, 0.3))
dev.off()
}

