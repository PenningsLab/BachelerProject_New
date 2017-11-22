setwd("/Users/pleuni/Documents/Git/bachelerProject")
source("Rscripts/baseRscript.R")

Bach.dat<-read.csv("Output/OverviewSelCoeff_BachelerFilter.csv")

head(Bach.dat)
breaks=seq(-4,0,by=.2)

if (TRUE){
pdf("Output/DFE_log.pdf")

for (typeofsite in c("syn", "nonsyn")){
par(mfrow=c(2,2))
i=1
for (wtnt in c("a", "t", "c", "g")){
datavector<- Bach.dat$EstSelCoeff[Bach.dat$TypeOfSite==typeofsite & Bach.dat$WTnt==wtnt]
hist(log10(datavector), breaks = breaks, xlim=c(-4,-0),col=cols[i], 
     main = paste0(wtnt, " : ", typeofsite), ylab="Count", xlab="",xaxt="n")

#abline(v=median(log10(datavector)),lwd=3,lty=2)

x=seq(-4,0,by=1)
axis(1, at=x,labels=paste("10^",x,sep=""), col.axis="black", las=1, line=-.2)
mtext(text = "Selection Coefficient", side = 1, line = 2, cex=.9)

i=i+1
}
}

dev.off()
}
