#Script to analyse the frequency data and associate with features. 
#Using the Bacheler et al data


setwd("~/Documents/Git/bachelerProject/Rscripts")
source('./baseRscript.R')
library(RColorBrewer)

read.table("../Output/freqPatTs_Bacheler.csv",sep=",",header=TRUE,row.names=1)->freqPatTs0
read.csv("../Output/OverviewSelCoeff_Bacheler.csv")->OverviewDF

ytxt<-"Estimated Selection Coefficient (cost)"
ylab<-c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0)

#Make a figure with the selection coefficients across Protease

OverviewDF$LC<-OverviewDF$lowerConf
OverviewDF$LC[OverviewDF$LC==0]<-OverviewDF$TSmutrate[OverviewDF$LC==0] #Use LC to make conf interval include 1 (= lethal)

for (ProRT in c("Pro","RT"))
    {
    xrange<-40:297; if (ProRT =="RT") xrange<-298:984
    xtxt<-"Position in Protease"; if (ProRT =="RT")  xtxt<-"Position in RT"
    xlab<-c(14,35,55,75,95); if (ProRT =="RT")  xlab<-seq(20,228,by=20)
    xlabat<-3*xlab-1; if (ProRT =="RT")  xlabat<-3*xlab-1+298
    RToffset = 0 ; if (ProRT =="RT") RToffset = 99
    
pdf(paste("../Output/EstSelCoeffConfIntervals",ProRT,".pdf",sep=""),width=12,height=8)

plottitle<-"#not G #SYN #no CPG" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)

mtext(side = 3,plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt!="g"&OverviewDF$makesCpG==0&OverviewDF$ProRT==ProRT,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])
a=1:length(data$num) #PSP 2017 June use this to add conf intervals
arrows(x0=data$num[a],y0=data$TSmutrate[a]/data$LC[a],x1=data$num[a],y1=data$TSmutrate[a]/data$upperConf[a],code=3,length = 0.05,angle=90)
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3,labels=ceiling(data$num[a]/3)-RToffset,cex=0.5)

plottitle<-"#not G #SYN #CPG" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt!="g"&OverviewDF$makesCpG==1&OverviewDF$ProRT==ProRT,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])
a=1:length(data$num) #PSP 2017 June use this to add conf intervals
arrows(x0=data$num[a],y0=data$TSmutrate[a]/data$LC[a],x1=data$num[a],y1=data$TSmutrate[a]/data$upperConf[a],code=3,length = 0.05,angle=90)
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3,labels=ceiling(data$num[a]/3)-RToffset,cex=0.5)

plottitle<-"#G #SYN" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="g"&OverviewDF$ProRT==ProRT,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])
a=1:length(data$num) #PSP 2017 June use this to add conf intervals
arrows(x0=data$num[a],y0=data$TSmutrate[a]/data$LC[a],x1=data$num[a],y1=data$TSmutrate[a]/data$upperConf[a],code=3,length = 0.05,angle=90)
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3,labels=ceiling(data$num[a]/3)-RToffset,cex=0.5)

plottitle<-"#A or T #NONSYN #no CPG #not drastic" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="nonsyn"&(OverviewDF$WTnt=="a"|OverviewDF$WTnt=="t")&OverviewDF$makesCpG==0&OverviewDF$bigAAChange==0&OverviewDF$ProRT==ProRT,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])
a=1:length(data$num) #PSP 2017 June use this to add conf intervals
arrows(x0=data$num[a],y0=data$TSmutrate[a]/data$LC[a],x1=data$num[a],y1=data$TSmutrate[a]/data$upperConf[a],code=3,length = 0.05,angle=90)
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3,labels=ceiling(data$num[a]/3)-RToffset,cex=0.5)

plottitle<-"#A or T #NONSYN #CPG #not drastic" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="nonsyn"&(OverviewDF$WTnt=="a"|OverviewDF$WTnt=="t")&OverviewDF$makesCpG==1&OverviewDF$bigAAChange==0&OverviewDF$ProRT==ProRT,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])
a=1:length(data$num) #PSP 2017 June use this to add conf intervals
arrows(x0=data$num[a],y0=data$TSmutrate[a]/data$LC[a],x1=data$num[a],y1=data$TSmutrate[a]/data$upperConf[a],code=3,length = 0.05,angle=90)
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3,labels=ceiling(data$num[a]/3)-RToffset,cex=0.5)

plottitle<-"#A or T #NONSYN #no CPG #drastic" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="nonsyn"&(OverviewDF$WTnt=="a"|OverviewDF$WTnt=="t")&OverviewDF$makesCpG==0&OverviewDF$bigAAChange==1&OverviewDF$ProRT==ProRT,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])
a=1:length(data$num) #PSP 2017 June use this to add conf intervals
arrows(x0=data$num[a],y0=data$TSmutrate[a]/data$LC[a],x1=data$num[a],y1=data$TSmutrate[a]/data$upperConf[a],code=3,length = 0.05,angle=90)
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3,labels=ceiling(data$num[a]/3)-RToffset,cex=0.5)

plottitle<-"#A or T #NONSYN #CPG #drastic" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="nonsyn"&(OverviewDF$WTnt=="a"|OverviewDF$WTnt=="t")&OverviewDF$makesCpG==1&OverviewDF$bigAAChange==1&OverviewDF$ProRT==ProRT,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])
a=1:length(data$num) #PSP 2017 June use this to add conf intervals
arrows(x0=data$num[a],y0=data$TSmutrate[a]/data$LC[a],x1=data$num[a],y1=data$TSmutrate[a]/data$upperConf[a],code=3,length = 0.05,angle=90)
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3,labels=ceiling(data$num[a]/3)-RToffset,cex=0.5)


plottitle<-"#C or G #NONSYN #not drastic" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="nonsyn"&(OverviewDF$WTnt=="c"|OverviewDF$WTnt=="g")&OverviewDF$bigAAChange==0&OverviewDF$ProRT==ProRT,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])
a=1:length(data$num) #PSP 2017 June use this to add conf intervals
arrows(x0=data$num[a],y0=data$TSmutrate[a]/data$LC[a],x1=data$num[a],y1=data$TSmutrate[a]/data$upperConf[a],code=3,length = 0.05,angle=90)
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3,labels=ceiling(data$num[a]/3)-RToffset,cex=0.5)

plottitle<-"#C or G #NONSYN #drastic" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="nonsyn"&(OverviewDF$WTnt=="c"|OverviewDF$WTnt=="g")&OverviewDF$bigAAChange==1&OverviewDF$ProRT==ProRT,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])
a=1:length(data$num) #PSP 2017 June use this to add conf intervals
arrows(x0=data$num[a],y0=data$TSmutrate[a]/data$LC[a],x1=data$num[a],y1=data$TSmutrate[a]/data$upperConf[a],code=3,length = 0.05,angle=90)
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3,labels=ceiling(data$num[a]/3)-RToffset,cex=0.5)

dev.off()

}


