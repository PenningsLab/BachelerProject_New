#Script to analyse the frequency data with conf intervals.
#Using the Bacheler et al data


setwd("~/Documents/Git/bachelerProject/Rscripts")
source('./baseRscript.R')
library(RColorBrewer)

read.table("../Output/freqPatTs_Bacheler.csv",sep=",",header=TRUE,row.names=1)->freqPatTs0
read.csv("../Output/OverviewSelCoeff_Bacheler.csv")->OverviewDF

ytxt<-"Estimated Selection Coefficient (cost)"
ylab<-c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0)

perc = 0.05 #(show 5% outliers)

#Make a figure with the selection coefficients across Protease

OverviewDF$LC<-OverviewDF$lowerConf
OverviewDF$LC[OverviewDF$LC==0]<-OverviewDF$TSmutrate[OverviewDF$LC==0] #Use LC to make conf interval include 1 (= lethal)
ListOutliers<-c()

xrange <-40:984
xtxt<-"AA Position in Protease / RT"
xlab<-c(c(15,35,55,75,95),seq(20,228,by=20));
xlabat<-3*xlab-1+c(rep(0,5),rep(298,11)); 
RToffset = 0 

pdf(paste("../Output/EstSelCoeffConfIntervals_005",".pdf",sep=""),width=12,height=8)

plottitle<-"#not G #SYN #no CPG" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1.2))
axis(1,at=xlabat,labels=xlab)
axis(2,at=ylab,labels=ylab)
abline(v = 297.5,lty=2); text(150,10^-5,"<--- Protease --->");text(600,10^-5,"<--- RT --->")

mtext(side = 3,plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt!="g"&OverviewDF$makesCpG==0,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11]) 
#show outliers
cutoff<-data$EstSelCoeff[order(data$EstSelCoeff)][ceiling(length(data$num)*(1-perc))]
points(data$num[data$EstSelCoeff>=cutoff],data$EstSelCoeff[data$EstSelCoeff>=cutoff],pch=21,bg="red")
data$num[data$EstSelCoeff<=cutoff]
#show outliers 
cutoff<-data$EstSelCoeff[order(data$EstSelCoeff)][ceiling(length(data$num)*(1-perc))] 
points(data$num[data$EstSelCoeff>=cutoff],data$EstSelCoeff[data$EstSelCoeff>=cutoff],pch=21,bg="red")
arrows(x0=data$num,y0=data$TSmutrate/data$LC,x1=data$num,y1=data$TSmutrate/data$upperConf,code=3,length = 0.05,angle=90)
#add locations for outliers
a=data$EstSelCoeff>=cutoff #PSP 2017 June use this to add conf intervals
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3^runif(length(data$num[a]),0.7,1.3),labels=data$AAnum[a],cex=0.3)
ListOutliers<-c(ListOutliers,data$num[a])

plottitle<-"#not G #SYN #CPG" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1.2))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)
abline(v = 297.5,lty=2); text(150,10^-5,"<--- Protease --->");text(600,10^-5,"<--- RT --->")

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt!="g"&OverviewDF$makesCpG==1,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])  
#show outliers 
cutoff<-data$EstSelCoeff[order(data$EstSelCoeff)][ceiling(length(data$num)*(1-perc))] 
points(data$num[data$EstSelCoeff>=cutoff],data$EstSelCoeff[data$EstSelCoeff>=cutoff],pch=21,bg="red")
arrows(x0=data$num,y0=data$TSmutrate/data$LC,x1=data$num,y1=data$TSmutrate/data$upperConf,code=3,length = 0.05,angle=90)
#add locations for outliers
a=data$EstSelCoeff>=cutoff #PSP 2017 June use this to add conf intervals
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3^runif(length(data$num[a]),0.9,1.1),labels=data$AAnum[a],cex=0.5)
ListOutliers<-c(ListOutliers,data$num[a])

plottitle<-"#G #SYN" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1.2))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)
abline(v = 297.5,lty=2); text(150,10^-5,"<--- Protease --->");text(600,10^-5,"<--- RT --->")

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="syn"&OverviewDF$WTnt=="g",]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])  
#show outliers 
cutoff<-data$EstSelCoeff[order(data$EstSelCoeff)][ceiling(length(data$num)*(1-perc))] 
points(data$num[data$EstSelCoeff>=cutoff],data$EstSelCoeff[data$EstSelCoeff>=cutoff],pch=21,bg="red")
arrows(x0=data$num,y0=data$TSmutrate/data$LC,x1=data$num,y1=data$TSmutrate/data$upperConf,code=3,length = 0.05,angle=90)
#add locations for outliers
a=data$EstSelCoeff>=cutoff #PSP 2017 June use this to add conf intervals
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3^runif(length(data$num[a]),0.9,1.1),labels=data$AAnum[a],cex=0.5)
ListOutliers<-c(ListOutliers,data$num[a])

plottitle<-"#A or T #NONSYN #no CPG #not drastic" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1.2))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)
abline(v = 297.5,lty=2); text(150,10^-5,"<--- Protease --->");text(600,10^-5,"<--- RT --->")

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="nonsyn"&(OverviewDF$WTnt=="a"|OverviewDF$WTnt=="t")&OverviewDF$makesCpG==0&OverviewDF$bigAAChange==0,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])  
#show outliers 
cutoff<-data$EstSelCoeff[order(data$EstSelCoeff)][ceiling(length(data$num)*(1-perc))] 
points(data$num[data$EstSelCoeff>=cutoff],data$EstSelCoeff[data$EstSelCoeff>=cutoff],pch=21,bg="red")
arrows(x0=data$num,y0=data$TSmutrate/data$LC,x1=data$num,y1=data$TSmutrate/data$upperConf,code=3,length = 0.05,angle=90)
#add locations for outliers
a=data$EstSelCoeff>=cutoff #PSP 2017 June use this to add conf intervals
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3^runif(length(data$num[a]),0.9,1.1),labels=data$AAnum[a],cex=0.5)
ListOutliers<-c(ListOutliers,data$num[a])

plottitle<-"#A or T #NONSYN #CPG #not drastic" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1.2))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)
abline(v = 297.5,lty=2); text(150,10^-5,"<--- Protease --->");text(600,10^-5,"<--- RT --->")

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="nonsyn"&(OverviewDF$WTnt=="a"|OverviewDF$WTnt=="t")&OverviewDF$makesCpG==1&OverviewDF$bigAAChange==0,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])  
#show outliers 
cutoff<-data$EstSelCoeff[order(data$EstSelCoeff)][ceiling(length(data$num)*(1-perc))] 
points(data$num[data$EstSelCoeff>=cutoff],data$EstSelCoeff[data$EstSelCoeff>=cutoff],pch=21,bg="red")
arrows(x0=data$num,y0=data$TSmutrate/data$LC,x1=data$num,y1=data$TSmutrate/data$upperConf,code=3,length = 0.05,angle=90)
#add locations for outliers
a=data$EstSelCoeff>=cutoff #PSP 2017 June use this to add conf intervals
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3^runif(length(data$num[a]),0.9,1.1),labels=data$AAnum[a],cex=0.5)
ListOutliers<-c(ListOutliers,data$num[a])

plottitle<-"#A or T #NONSYN #no CPG #drastic" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1.2))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)
abline(v = 297.5,lty=2); text(150,10^-5,"<--- Protease --->");text(600,10^-5,"<--- RT --->")

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="nonsyn"&(OverviewDF$WTnt=="a"|OverviewDF$WTnt=="t")&OverviewDF$makesCpG==0&OverviewDF$bigAAChange==1,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])  
#show outliers 
cutoff<-data$EstSelCoeff[order(data$EstSelCoeff)][ceiling(length(data$num)*(1-perc))] 
points(data$num[data$EstSelCoeff>=cutoff],data$EstSelCoeff[data$EstSelCoeff>=cutoff],pch=21,bg="red")
arrows(x0=data$num,y0=data$TSmutrate/data$LC,x1=data$num,y1=data$TSmutrate/data$upperConf,code=3,length = 0.05,angle=90)
#add locations for outliers
a=data$EstSelCoeff>=cutoff #PSP 2017 June use this to add conf intervals
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3^runif(length(data$num[a]),0.9,1.1),labels=data$AAnum[a],cex=0.5)
ListOutliers<-c(ListOutliers,data$num[a])

plottitle<-"#A or T #NONSYN #CPG #drastic" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1.2))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)
abline(v = 297.5,lty=2); text(150,10^-5,"<--- Protease --->");text(600,10^-5,"<--- RT --->")

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="nonsyn"&(OverviewDF$WTnt=="a"|OverviewDF$WTnt=="t")&OverviewDF$makesCpG==1&OverviewDF$bigAAChange==1,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11])  
#show outliers 
cutoff<-data$EstSelCoeff[order(data$EstSelCoeff)][ceiling(length(data$num)*(1-perc))] 
points(data$num[data$EstSelCoeff>=cutoff],data$EstSelCoeff[data$EstSelCoeff>=cutoff],pch=21,bg="red")
arrows(x0=data$num,y0=data$TSmutrate/data$LC,x1=data$num,y1=data$TSmutrate/data$upperConf,code=3,length = 0.05,angle=90)
#add locations for outliers
a=data$EstSelCoeff>=cutoff #PSP 2017 June use this to add conf intervals
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3^runif(length(data$num[a]),0.9,1.1),labels=data$AAnum[a],cex=0.5)
ListOutliers<-c(ListOutliers,data$num[a])

plottitle<-"#C or G #NONSYN #not drastic" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1.2))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)
abline(v = 297.5,lty=2); text(150,10^-5,"<--- Protease --->");text(600,10^-5,"<--- RT --->")

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="nonsyn"&(OverviewDF$WTnt=="c"|OverviewDF$WTnt=="g")&OverviewDF$bigAAChange==0,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11]) 
#show outliers 
cutoff<-data$EstSelCoeff[order(data$EstSelCoeff)][ceiling(length(data$num)*(1-perc))] 
points(data$num[data$EstSelCoeff>=cutoff],data$EstSelCoeff[data$EstSelCoeff>=cutoff],pch=21,bg="red")
arrows(x0=data$num,y0=data$TSmutrate/data$LC,x1=data$num,y1=data$TSmutrate/data$upperConf,code=3,length = 0.05,angle=90)
#add locations for outliers
a=data$EstSelCoeff>=cutoff #PSP 2017 June use this to add conf intervals
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3^runif(length(data$num[a]),0.9,1.1),labels=data$AAnum[a],cex=0.5)
ListOutliers<-c(ListOutliers,data$num[a])

plottitle<-"#C or G #NONSYN #drastic" #make plot
plot(OverviewDF$num[xrange],OverviewDF$EstSelCoeff[xrange],log="y", ylab=ytxt, xlab =xtxt, xaxt="n",yaxt="n",col="darkgrey",t="c",pch=".", ylim=c(0.9*10^-5,1.5))
axis(1,at=xlabat,labels=xlab);axis(2,at=ylab,labels=ylab)
abline(v = 297.5,lty=2); text(150,10^-5,"<--- Protease --->");text(600,10^-5,"<--- RT --->")

mtext(side = 3, plottitle)
data<-OverviewDF[OverviewDF$TypeOfSite=="nonsyn"&(OverviewDF$WTnt=="c"|OverviewDF$WTnt=="g")&OverviewDF$bigAAChange==1,]
points(data$num,data$EstSelCoeff,pch=21,bg=brewer.pal(11, "Spectral")[11]) 
#show outliers 
cutoff<-data$EstSelCoeff[order(data$EstSelCoeff)][ceiling(length(data$num)*(1-perc))] 
points(data$num[data$EstSelCoeff>=cutoff],data$EstSelCoeff[data$EstSelCoeff>=cutoff],pch=21,bg="red")
arrows(x0=data$num,y0=data$TSmutrate/data$LC,x1=data$num,y1=data$TSmutrate/data$upperConf,code=3,length = 0.05,angle=90)
#add locations for outliers
a=data$EstSelCoeff>=cutoff #PSP 2017 June use this to add conf intervals
text(data$num[a],y=data$TSmutrate[a]/data$LC[a]*1.3^runif(length(data$num[a]),0.9,1.1),labels=data$AAnum[a],cex=0.5)
ListOutliers<-c(ListOutliers,data$num[a])

#text(data$num[data$EstSelCoeff>=cutoff],
#    y=data$TSmutrate[data$EstSelCoeff>=cutoff]/data$LC[data$EstSelCoeff>=cutoff]*1.3^runif(length(data$num[data$EstSelCoeff>=cutoff]),0.5,2),labels=ceiling(data$num[a]/3)-RToffset,cex=0.5)

dev.off()

OutLiers<-cbind(
OverviewDF[sort(ListOutliers),c(7,2,7,10,15,11,12,13,14)],
round(OverviewDF[sort(ListOutliers),c(9)],3)
)
names(OutLiers)[10]<-"SelCoeff"
for (i in 1:length(OutLiers[,1])){
    OutLiers$WTnt.1[i]<-transition(as.character(OutLiers$WTnt[i]))
    }

print(OutLiers)

write.csv(OutLiers,file="../Output/ListOfOutliers.csv",row.names = FALSE)
