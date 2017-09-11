
#Prep data
if (TRUE){
    setwd("~/Documents/Git/bachelerProject/Rscripts")
    source('./baseRscript.R')
    library(scales)
    library(plotrix)
    library(RColorBrewer)
    #July 2017 now read freqPatTs_Bacheler_Threshold05.csv  and OverviewSelCoeff_BachelerFilter.csv
    read.table("../Output/freqPatTs_Bacheler_Threshold1.csv",sep=",",header=TRUE,row.names=1)->freqPatTs0
    read.csv("../Output/OverviewSelCoeff_BachelerFilter.csv")->OverviewDF
}

par(mar=c(2,2,2,2))
#Ok, I want to know whether there is an effect of location. 
#What we do to test this is to create a sliding window average and calculate the variance of the average of the windows. 
#We'll do this first for the non-syn sites
WindowSize=20
NonSynSites<-which(OverviewDF$TypeOfSite=="nonsyn")
SynSites<-which(OverviewDF$TypeOfSite=="syn")

NonSynWindowMeans<-c()
for (i in 1:(length(NonSynSites)-(WindowSize-1))){
    Window<-NonSynSites[i:(i+(WindowSize-1))] #windows of size 10
    NonSynWindowMeans<-c(NonSynWindowMeans,mean(OverviewDF$MeanFreq[Window]))
}
RealNonSynWindowVar<-var(NonSynWindowMeans)

SynWindowMeans<-c()
for (i in 1:(length(SynSites)-(WindowSize-1))){
    Window<-SynSites[i:(i+(WindowSize-1))] #windows of size 10
    SynWindowMeans<-c(SynWindowMeans,mean(OverviewDF$MeanFreq[Window]))
}
RealSynWindowVar<-var(SynWindowMeans)

#Prepare randomizations
AsitesNonSyn<-which(consensusB=="a"&OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$bigAAChange==0)
CsitesNonSyn<-which(consensusB=="c"&OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$bigAAChange==0)
GsitesNonSyn<-which(consensusB=="g"&OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$bigAAChange==0)
TsitesNonSyn<-which(consensusB=="t"&OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$bigAAChange==0)
AsitesNonSynDr<-which(consensusB=="a"&OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$bigAAChange==1)
CsitesNonSynDr<-which(consensusB=="c"&OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$bigAAChange==1)
GsitesNonSynDr<-which(consensusB=="g"&OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$bigAAChange==1)
TsitesNonSynDr<-which(consensusB=="t"&OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$bigAAChange==1)
AsitesSyn<-which(consensusB=="a"&OverviewDF$TypeOfSite=="syn")
CsitesSyn<-which(consensusB=="c"&OverviewDF$TypeOfSite=="syn")
GsitesSyn<-which(consensusB=="g"&OverviewDF$TypeOfSite=="syn")
TsitesSyn<-which(consensusB=="t"&OverviewDF$TypeOfSite=="syn")

ListOfRandomVarsNonSynSites<-c()
ListOfRandomVarsSynSites<-c()
for (j in 1:1000){
    FreqsRandom<-OverviewDF$MeanFreq
    FreqsRandom[AsitesNonSyn]<-OverviewDF$MeanFreq[sample(AsitesNonSyn)]
    FreqsRandom[CsitesNonSyn]<-OverviewDF$MeanFreq[sample(CsitesNonSyn)]
    FreqsRandom[GsitesNonSyn]<-OverviewDF$MeanFreq[sample(GsitesNonSyn)]
    FreqsRandom[TsitesNonSyn]<-OverviewDF$MeanFreq[sample(TsitesNonSyn)]
    FreqsRandom[AsitesNonSynDr]<-OverviewDF$MeanFreq[sample(AsitesNonSynDr)]
    FreqsRandom[CsitesNonSynDr]<-OverviewDF$MeanFreq[sample(CsitesNonSynDr)]
    FreqsRandom[GsitesNonSynDr]<-OverviewDF$MeanFreq[sample(GsitesNonSynDr)]
    FreqsRandom[TsitesNonSynDr]<-OverviewDF$MeanFreq[sample(TsitesNonSynDr)]
    FreqsRandom[AsitesSyn]<-OverviewDF$MeanFreq[sample(AsitesSyn)]
    FreqsRandom[CsitesSyn]<-OverviewDF$MeanFreq[sample(CsitesSyn)]
    FreqsRandom[GsitesSyn]<-OverviewDF$MeanFreq[sample(GsitesSyn)]
    FreqsRandom[TsitesSyn]<-OverviewDF$MeanFreq[sample(TsitesSyn)]
    #NonSyn
    RanNonSynWindowMeans<-c()
    for (i in 1:(length(NonSynSites)-(WindowSize-1))){
        Window<-NonSynSites[i:(i+(WindowSize-1))] #windows of size 10
        RanNonSynWindowMeans<-c(RanNonSynWindowMeans,mean(FreqsRandom[Window]))
    }
    ListOfRandomVarsNonSynSites<-c(ListOfRandomVarsNonSynSites,var(RanNonSynWindowMeans))  
    #Synonymous
    RanSynWindowMeans<-c()
    for (i in 1:(length(SynSites)-(WindowSize-1))){
        Window<-SynSites[i:(i+(WindowSize-1))] #windows of size 10
        RanSynWindowMeans<-c(RanSynWindowMeans,mean(FreqsRandom[Window]))
    }
    ListOfRandomVarsSynSites<-c(ListOfRandomVarsSynSites,var(RanSynWindowMeans))  
    }

hist(ListOfRandomVarsSynSites,xlim=c(0,RealNonSynWindowVar*2))
abline(v=RealNonSynWindowVar,col=2,lwd=2)
hist(ListOfRandomVarsSynSites,xlim=c(0,RealSynWindowVar*2))
abline(v=RealSynWindowVar,col=2,lwd=2)

