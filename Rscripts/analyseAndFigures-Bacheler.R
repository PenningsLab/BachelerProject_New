#Script to analyse the frequency data and associate with features. 
#Using the Bacheler et al data

#This script does / makes
#Tests whether non syn muts, syn muts and nonsense muts are different in freq
#Make EstSelCoeffRT.pdf
#Make EstSelCoeffPRO.pdf
#Make SingleSiteFrequencySpectraPRO_58.pdf


#* Read the csv files 
#* Perform necessary calcuations
#* Plot results (eventually new script)

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
#Test whether non syn muts, syn muts and nonsense muts are different in freq

FreqsSyn<-OverviewDF$MeanFreq[OverviewDF$TypeOfSite=="syn"]
FreqsNonSyn<-OverviewDF$MeanFreq[OverviewDF$TypeOfSite=="nonsyn"]
FreqsStop<-OverviewDF$MeanFreq[OverviewDF$TypeOfSite=="stop"]

wilcox.test(FreqsSyn, FreqsNonSyn,alternative = "greater", paired = FALSE)
wilcox.test(FreqsNonSyn,FreqsStop,alternative = "greater", paired = FALSE)

#Make a figure with the selection coefficients across Pol
if (TRUE){
#pdf("../Output/EstSelCoeffPRO_aug2017.pdf",width=12,height=8)
png("../Output/EstSelCoeffPRO_aug2017.png",width=12,height=8,units="in",res=100)
    par(mfrow=c(1,1))
#make log scale
maxnuc=984
plot(OverviewDF$num[40:maxnuc],OverviewDF$EstSelCoeff[40:maxnuc],
     log="y", ylab="Estimated Selection Coefficient (cost)", 
     xlab = "Position in Protease                                                                                 Position in RT                                                             ", 
     xaxt="n",yaxt="n", 
     col="darkgrey",t="n",pch=".", ylim=c(5*10^-4,1))
axis(1,at=c(3*seq(15,95,by=20)-1,296+30),labels=c(seq(15,95,by=20),""))
axis(1,at=3*seq(109,349,by=20)-1,labels=seq(109-99,349-99,by=20))
axis(2,at=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),labels=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0))
abline(v=297.5,lwd=4,col="grey61")

for(i in 1:5){
    abline(h = 1:10 * 10^(-i), col = "gray61")
}

cols <- brewer.pal(6, "Set2")[c(1, 2, 3, 6)]
for (i in 40:maxnuc){
    c=0
    if (OverviewDF$TypeOfSite[i]=="stop"&OverviewDF$WTnt[i]%in%c("g","c")) {c=1;p=21}
    if (OverviewDF$TypeOfSite[i]=="syn"&OverviewDF$WTnt[i]%in%c("g","c")) {c=cols[1];p=21}
    if (OverviewDF$TypeOfSite[i]=="syn"&OverviewDF$WTnt[i]%in%c("a","t")) {c=cols[1];p=21}
    if (OverviewDF$TypeOfSite[i]=="nonsyn"&OverviewDF$WTnt[i]%in%c("c","g")) {c=cols[2];p=21}
    if (OverviewDF$TypeOfSite[i]=="nonsyn"&OverviewDF$WTnt[i]%in%c("a","t")) {c=cols[4];p=21}
    if (c!=0) points(OverviewDF$num[i],OverviewDF$EstSelCoeff[i],pch=p,col=1,
                     bg=rgb(red=col2rgb(c)[1]/255,
                           green=col2rgb(c)[2]/255,
                           blue=col2rgb(c)[3]/255,
                           maxColorValue = 1,alpha=0.8),
   cex=2)
}


rect(270*3, 0.0004, 315*3, 0.0017, density = NULL, angle = 45,col=alpha("white",0.5))

#legend
points(275*3,0.001/0.7,pch=21,bg=1,col=1,cex=2)
text(290*3,0.001/0.7,"nonsense")
points(275*3,0.001,pch=21,bg=cols[2],col=1,cex=2)
text(295*3,0.001,"non-syn, C/G")
points(275*3,0.001*0.7,pch=21,bg=cols[4],col=1,cex=2)
text(295*3,0.001*0.7,"non-syn, A/T")
points(275*3,0.001*0.5,pch=21,bg=cols[1],col=1,cex=2)
text(295*3,0.001*0.5,"synonymous")

dev.off()
}


#Make a figure with the selection coefficients across Protease

if(FALSE){
#pdf("../Output/EstSelCoeffPRO.pdf",width=12,height=8)
par(mfrow=c(1,1))
#make log scale
plot(OverviewDF$num[40:297],OverviewDF$EstSelCoeff[40:297],
     log="y", ylab="Estimated Selection Coefficient (cost)", xlab = "Position in Protease", 
     xaxt="n",yaxt="n",
     col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1))
axis(1,at=3*c(14,35,55,75,95)-1,labels=c(14,35,55,75,95))
axis(2,at=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),labels=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0))
#Add colors for different types of sites

#A #SYN #no CPG
data<-OverviewDF[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==0&OverviewDF$num<=297,]

points(data$num,data$EstSelCoeff,
       pch=21,bg=brewer.pal(11, "Spectral")[11])
a=1:length(data$num)
#PSP 2017 June use this to add conf intervals
#arrows(x0=data$num[a],y0=data$TSmutrate[a]/data$lowerConf[a],
#       x1=data$num[a],y1=data$TSmutrate[a]/data$upperConf[a],
#       code=3,length = 0.05,angle=90)

#legend
points(45,2*10^-5,pch=21,bg=brewer.pal(11, "Spectral")[11])
text(48,2*10^-5,pos=4,"A/T/C, syn, no CpG")

#T #SYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==0&OverviewDF$num<=297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==0&OverviewDF$num<=297],
       pch=21,bg=brewer.pal(11, "Spectral")[11])

#C #SYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="c"&OverviewDF$makesCpG==0&OverviewDF$num<=297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="c"&OverviewDF$makesCpG==0&OverviewDF$num<=297],
       pch=21,bg=brewer.pal(11, "Spectral")[11])

#A #SYN #CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==1&OverviewDF$num<=297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==1&OverviewDF$num<=297],
       pch=24,bg=brewer.pal(11, "Spectral")[10])

points(45,1.4*10^-5,pch=24,
       bg=brewer.pal(11, "Spectral")[10])
text(48,1.4*10^-5,pos=4,"A/T, syn, CpG")

#T #SYN #CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==1&OverviewDF$num<=297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==1&OverviewDF$num<=297],pch=24,bg=brewer.pal(11, "Spectral")[10])

#G #SYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="g"&OverviewDF$makesCpG==0&OverviewDF$num<=297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="g"&OverviewDF$makesCpG==0&OverviewDF$num<=297],pch=25,bg=brewer.pal(11, "Spectral")[8],col=brewer.pal(11, "RdYlGn")[11])

points(45,1*10^-5,pch=25,bg=brewer.pal(11, "Spectral")[8],col=brewer.pal(11, "RdYlGn")[11])
text(48,1*10^-5,pos=4,"G, syn")

#A #NONSYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==0],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==0],pch=21,col=2,bg=brewer.pal(11, "Spectral")[7])

points(100,2*10^-5,pch=21,bg=brewer.pal(11, "Spectral")[7],col=2)
text(103,2*10^-5,pos=4,"A/T, non-syn, no drastic AA change")

#T #NONSYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==0],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==0],pch=21,col=2,bg=brewer.pal(11, "Spectral")[7])

#A #NONSYN #CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==1&OverviewDF$num<=297&OverviewDF$bigAAChange==0],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==1&OverviewDF$num<=297&OverviewDF$bigAAChange==0],pch=21,col=2,bg=brewer.pal(11, "Spectral")[7])

#T #NONSYN #CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==1&OverviewDF$num<=297&OverviewDF$bigAAChange==0],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==1&OverviewDF$num<=297&OverviewDF$bigAAChange==0],pch=21,col=2,bg=brewer.pal(11, "Spectral")[7])

#With big AA change

#A #NONSYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==1],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==1],pch=22,col=2,bg=brewer.pal(11, "Spectral")[5])

points(100,1.4*10^-5,pch=22,bg=brewer.pal(11, "Spectral")[5],col=2)
text(103,1.4*10^-5,pos=4,"A/T, non-syn, drastic AA change")

#T #NONSYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==1],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==1],pch=22,col=2,bg=brewer.pal(11, "Spectral")[5])

#A #NONSYN #CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==1&OverviewDF$num<=297&OverviewDF$bigAAChange==1],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==1&OverviewDF$num<=297&OverviewDF$bigAAChange==1],pch=22,col=2,bg=brewer.pal(11, "Spectral")[5])

#T #NONSYN #CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==1&OverviewDF$num<=297&OverviewDF$bigAAChange==1],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==1&OverviewDF$num<=297&OverviewDF$bigAAChange==1],pch=22,col=2,bg=brewer.pal(11, "Spectral")[5])


#no big AA change

#C #NONSYN no big AA change
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="c"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==0],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="c"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==0],pch=21,col=2,bg=brewer.pal(11, "Spectral")[3])

points(100,1*10^-5,pch=21,bg=brewer.pal(11, "Spectral")[3],col=2)
text(103,1*10^-5,pos=4,"C, non-syn, no drastic AA change")

#G #NONSYN no big AA change
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="g"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==0],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="g"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==0],pch=21,col=2,bg=brewer.pal(11, "Spectral")[1])

points(200,2*10^-5,pch=21,bg=brewer.pal(11, "Spectral")[1],col=2)
text(203,2*10^-5,pos=4,"G, non-syn, no drastic AA change")

#C #NONSYN bigAAChange
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="c"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==1],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="c"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==1],pch=22,col=2,bg=brewer.pal(11, "Spectral")[3])

points(200,1.4*10^-5,pch=22,bg=brewer.pal(11, "Spectral")[3],col=2)
text(203,1.4*10^-5,pos=4,"C, non-syn, drastic AA change")

#G #NONSYN bigAAChange
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="g"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==1],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="g"&OverviewDF$makesCpG==0&OverviewDF$num<=297&OverviewDF$bigAAChange==1],pch=22,col=2,bg=brewer.pal(11, "Spectral")[1])

points(200,1*10^-5,pch=22,bg=brewer.pal(11, "Spectral")[1],col=2)
text(203,1*10^-5,pos=4,"G, non-syn, drastic AA change")

#nonsense
points(OverviewDF$num[OverviewDF$TypeOfSite=="stop"&OverviewDF$num<=297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="stop"&OverviewDF$num<=297],pch=22,col=2,bg=1)

points(280,2*10^-5,pch=22,bg=1,col=1)
text(283,2*10^-5,pos=4,"nonsense")

#resistance
points(OverviewDF$num[OverviewDF$TypeOfSite=="res"&OverviewDF$num<=297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="res"&OverviewDF$num<=297],pch=8,col=1,bg=1)

points(280,1.4*10^-5,pch=8,bg=1,col=1)
text(283,1.4*10^-5,pos=4,"resistance")

#dev.off()
}

#Make a figure with the selection coefficients across RT
if (FALSE){
pdf("../Output/EstSelCoeffRT.pdf",width=12,height=8)
par(mfrow=c(1,1))
#make log scale
plot(OverviewDF$num[298:984],OverviewDF$EstSelCoeff[298:984],log="y", 
     ylab="Estimated Selection Coefficient (cost)", xlab = "Position in Reverse Transcriptase", xaxt="n",yaxt="n",col="darkgrey",
     t="c",pch=".",
     ylim=c(0.9*10^-5,1))
axis(1,at=298+3*seq(20,228,by=20)-1,labels=seq(20,228,by=20))
axis(2,at=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),labels=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0))
#Add colors for different types of sites

#A #SYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==0&OverviewDF$num>297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==0&OverviewDF$num>297],pch=21,bg=brewer.pal(11, "Spectral")[11])

points(298+45,2*10^-5,pch=21,bg=brewer.pal(11, "Spectral")[11])
text(298+48,2*10^-5,pos=4,"A/T/C, syn, no CpG")

#T #SYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==0&OverviewDF$num>297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==0&OverviewDF$num>297],pch=21,bg=brewer.pal(11, "Spectral")[11])

#C #SYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="c"&OverviewDF$makesCpG==0&OverviewDF$num>297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="c"&OverviewDF$makesCpG==0&OverviewDF$num>297],pch=21,bg=brewer.pal(11, "Spectral")[11])

#A #SYN #CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==1&OverviewDF$num>297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==1&OverviewDF$num>297],pch=24,bg=brewer.pal(11, "Spectral")[10])

points(298+45,1.4*10^-5,pch=24,bg=brewer.pal(11, "Spectral")[10])
text(298+48,1.4*10^-5,pos=4,"A/T, syn, CpG")

#T #SYN #CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==1&OverviewDF$num>297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==1&OverviewDF$num>297],pch=24,bg=brewer.pal(11, "Spectral")[10])

#G #SYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="g"&OverviewDF$makesCpG==0&OverviewDF$num>297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="g"&OverviewDF$makesCpG==0&OverviewDF$num>297],pch=25,bg=brewer.pal(11, "Spectral")[8],col=brewer.pal(11, "RdYlGn")[11])

points(298+45,1*10^-5,pch=25,bg=brewer.pal(11, "Spectral")[8],col=brewer.pal(11, "RdYlGn")[11])
text(298+48,1*10^-5,pos=4,"G, syn")

#A #NONSYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==0],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==0],pch=21,col=2,bg=brewer.pal(11, "Spectral")[7])

points(298+180,2*10^-5,pch=21,bg=brewer.pal(11, "Spectral")[7],col=2)
text(298+183,2*10^-5,pos=4,"A/T, non-syn, no drastic AA change")

#T #NONSYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==0],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==0],pch=21,col=2,bg=brewer.pal(11, "Spectral")[7])

#A #NONSYN #CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==1&OverviewDF$num>297&OverviewDF$bigAAChange==0],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==1&OverviewDF$num>297&OverviewDF$bigAAChange==0],pch=21,col=2,bg=brewer.pal(11, "Spectral")[7])

#T #NONSYN #CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==1&OverviewDF$num>297&OverviewDF$bigAAChange==0],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==1&OverviewDF$num>297&OverviewDF$bigAAChange==0],pch=21,col=2,bg=brewer.pal(11, "Spectral")[7])

#With big AA change

#A #NONSYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==1],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==1],pch=22,col=2,bg=brewer.pal(11, "Spectral")[5])

points(298+180,1.4*10^-5,pch=22,bg=brewer.pal(11, "Spectral")[5],col=2)
text(298+183,1.4*10^-5,pos=4,"A/T, non-syn, drastic AA change")

#T #NONSYN #no CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==1],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==1],pch=22,col=2,bg=brewer.pal(11, "Spectral")[5])

#A #NONSYN #CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==1&OverviewDF$num>297&OverviewDF$bigAAChange==1],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="a"&OverviewDF$makesCpG==1&OverviewDF$num>297&OverviewDF$bigAAChange==1],pch=22,col=2,bg=brewer.pal(11, "Spectral")[5])

#T #NONSYN #CPG
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==1&OverviewDF$num>297&OverviewDF$bigAAChange==1],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="t"&OverviewDF$makesCpG==1&OverviewDF$num>297&OverviewDF$bigAAChange==1],pch=22,col=2,bg=brewer.pal(11, "Spectral")[5])


#no big AA change

#C #NONSYN no big AA change
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="c"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==0],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="c"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==0],pch=21,col=2,bg=brewer.pal(11, "Spectral")[3])

points(298+180,1*10^-5,pch=21,bg=brewer.pal(11, "Spectral")[3],col=2)
text(298+183,1*10^-5,pos=4,"C, non-syn, no drastic AA change")

#G #NONSYN no big AA change
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="g"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==0],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="g"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==0],pch=21,col=2,bg=brewer.pal(11, "Spectral")[1])

points(298+400,2*10^-5,pch=21,bg=brewer.pal(11, "Spectral")[1],col=2)
text(298+403,2*10^-5,pos=4,"G, non-syn, no drastic AA change")

#C #NONSYN bigAAChange
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="c"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==1],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="c"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==1],pch=22,col=2,bg=brewer.pal(11, "Spectral")[3])

points(298+400,1.4*10^-5,pch=22,bg=brewer.pal(11, "Spectral")[3],col=2)
text(298+403,1.4*10^-5,pos=4,"C, non-syn, drastic AA change")

#G #NONSYN bigAAChange
points(OverviewDF$num[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="g"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==1],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="nonsyn"&OverviewDF$WTnt=="g"&OverviewDF$makesCpG==0&OverviewDF$num>297&OverviewDF$bigAAChange==1],pch=22,col=2,bg=brewer.pal(11, "Spectral")[1])

points(298+400,1*10^-5,pch=22,bg=brewer.pal(11, "Spectral")[1],col=2)
text(298+403,1*10^-5,pos=4,"G, non-syn, drastic AA change")

#nonsense
points(OverviewDF$num[OverviewDF$TypeOfSite=="stop"&OverviewDF$num>297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="stop"&OverviewDF$num>297],pch=22,col=2,bg=1)

points(298+620,2*10^-5,pch=22,bg=1,col=1)
text(298+623,2*10^-5,pos=4,"nonsense")

#resistance
points(OverviewDF$num[OverviewDF$TypeOfSite=="res"&OverviewDF$num>297],OverviewDF$EstSelCoeff[OverviewDF$TypeOfSite=="res"&OverviewDF$num>297],pch=8,col=1,bg=1)

points(298+620,1.4*10^-5,pch=8,bg=1,col=1)
text(298+623,1.4*10^-5,pos=4,"resistance")

dev.off()
}

#Make a figure with single site frequency spectra for Protease AA 58
if (TRUE){
pdf("../Output/SingleSiteFrequencySpectraPRO_58_July2017.pdf",width=8,height=4)
zerobar=50; h2=22; x1=0.25
cols <- c(0,brewer.pal(6, "Set2")[c(2, 1)])
par(mfrow=c(2,3))
par(mar = c(1,3,4,2))
for (i in 172:174){
    #first create empty plot with title
    if (i == 172){
        #par(fig=c(0,2/3,0,1))
        t=paste("nonsense mutation",sep="")
        hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),yaxt="n",
             col=cols[1],border=0,
             #    main= bquote(paste(.(t),(C %->% T ))), cex=1.3,
             main="",cex=1.2,
             xlab="Frequency", ylab="Count",cex.lab=1.4)
        title(t,cex=1.2,line=0)
        text(x1,30,"observed data",cex=1.3)
        text(x1,h2,"(C172T)",cex=1.3)
    }
    if (i == 173){
        t=paste("         non-synonymous mutation",sep="")
        #    t=paste("Protease: site ", i,"\n non-synonymous mutation",sep="")
        hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),yaxt="n",
             col=cols[2],border=0,
             #    main = bquote(paste(.(t),(A %->% G ))), cex=1.3,
             main= "", cex=1.2,
             xlab="Frequency", ylab="Count",cex.lab=1.4)
        title(t,cex=1.2,line=0)
        text(x1,30,"observed data",cex=1.3)
        text(x1,h2,"(A173G)",cex=1.3)
        
    }
    if (i == 174){
        t=paste("   synonymous mutation",sep="")
        #    t=paste("Protease: site ", i,"\n synonymous mutation",sep="")
        hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),yaxt="n",
             col=cols[3],border=0,
             #    main = bquote(paste(.(t),(G %->% A ))), cex=1.3,
             main= "", cex=1.2,
             xlab="Frequency", ylab="Count",cex.lab=1.4)
        title(t,cex=1.2,line=0)
        text(x1,30,"observed data",cex=1.3)
        text(x1,h2,"(G174A)",cex=1.3)
        
    }
    #Next, show  0 bar
    if (i == 172){
        hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),
             yaxt="n",col=OverviewDF$color[which(OverviewDF$num==i)],add=T)}
    if (i == 173){
        hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),
             yaxt="n",col=cols[2],add=T)}
    if (i == 174){
        hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),
             yaxt="n",col=cols[3],add=T)}
    
    #next show all data (unfiltered), but only until 50 for 0 cat
    if (i == 172){
        hist(c(rep(0,min(zerobar-10,length(which(freqPatTs0[,i]<0.02)))),freqPatTs0[,i][which(freqPatTs0[,i]>0)]),
             breaks=seq(0,1,by=0.02),add=T,
             col=OverviewDF$color[which(OverviewDF$num==i)])}
    if (i == 173){
        hist(c(rep(0,min(zerobar-10,length(which(freqPatTs0[,i]<0.02)))),freqPatTs0[,i][which(freqPatTs0[,i]>0)]),
             breaks=seq(0,1,by=0.02),add=T,
             col=cols[2])}
    if (i == 174){
        hist(c(rep(0,min(zerobar-10,length(which(freqPatTs0[,i]<0.02)))),freqPatTs0[,i][which(freqPatTs0[,i]>0)]),
             breaks=seq(0,1,by=0.02),add=T,
             col=cols[3])}
    
    axis(2,labels = c(10,20,30,max(zerobar,length(which(freqPatTs0[,i]<0.02)))), 
         at = c(10,20,30,zerobar), las=1)
    if (length(which(freqPatTs0[,i]<0.02))>=zerobar){
        axis.break(axis=2,breakpos=zerobar-10,bgcol="white",breakcol="black",style="slash",brw=0.02)
        points(c(0.01,0.02),c(zerobar-10,zerobar-10),pch=15,cex=2.5,col="white")
    }else{axis(2,labels = zerobar-10,at=zerobar-10,las=1)}
    
}
#next, show simulated frequencies
#    par(fig=c(2/3,1,0,1), new = TRUE)
par(mar = c(4.5,3,1,2))
for (i in 172:174){
    if (i ==172)Freqs<-read.csv("../Output/SimFreqs172.csv",row.names=1)[1][,1]
    if (i ==173)Freqs<-read.csv("../Output/SimFreqs173.csv",row.names=1)[1][,1]
    if (i ==174)Freqs<-read.csv("../Output/SimFreqs174.csv",row.names=1)[1][,1]
    t=paste("simulated data",sep="")
    hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),yaxt="n",
         col=cols[1],border=0,
         main="",cex=1.2,
         xlab="Frequency", ylab="Count",cex.lab=1.4)
    #title(t,cex=1.2,line=0)
    text(x1,30,"simulated data",cex=1.3)
    if (i ==172)text(x1,h2,"(s=1)",cex=1.3)
    if (i ==173)text(x1,h2,paste("(s=",round(OverviewDF$EstSelCoeff[173],3),")",sep=""),cex=1.3)
    if (i ==174)text(x1,h2,paste("(s=",round(OverviewDF$EstSelCoeff[174],3),")",sep=""),cex=1.3)
    
    hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),
         yaxt="n",col=brewer.pal(4, "Set2")[4],add=T)
    hist(c(rep(0,
               min(zerobar,length(which(Freqs<0.02)))
               ),
           Freqs[which(Freqs>=0.02)]),
         breaks=seq(0,1,by=0.02),add=T,
        col=brewer.pal(4, "Set2")[4])
    axis(2,labels = c(10,20,30,max(zerobar,length(which(Freqs<0.02)))), 
     at = c(10,20,30,zerobar), las=1)
    if (length(which(Freqs<0.02))>=zerobar){
    axis.break(axis=2,breakpos=zerobar-10,bgcol="white",breakcol="black",style="slash",brw=0.02)
    points(c(0.01,0.02),c(zerobar-10,zerobar-10),pch=15,cex=2.5,col="white")
    }else{axis(2,labels = zerobar-10,at=zerobar-10,las=1)}
}

dev.off()
}
