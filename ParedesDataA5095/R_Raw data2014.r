#This R script compares expected with real Freq Dist. I am not sure how useful it still is (PP 2016 Feb).
#Maybe it is useful to get the data from only the random subcohort.
#MeanFreqData has also data from the patients who were picked because they failed treatment.

#to estimate the number of particles that were measured, I assume that the estimated frequencies are the result of binomial sampling.
#if the data are x=c(0.02,0.03) for a given patient (2% and 3%), then the estimated number of particles is 
#nest = mean(x)*(1-mean(x))/var(x)
#which in this case is 487. 
#now i will do this for all patients for K103N and Y181C
###2014: now we know that the two estimates were not independent. the project was put on ice. 
###2014: maybe the data can be useful for an NIH proposal

###########################
#Get list of PIDs of inds from the random subcohort
###########################
read.table("DataJonLi.csv",sep="\t",header=TRUE)->realdata
read.csv("MeanFreqData.csv",sep=",",header=TRUE)->Mutdata

dataPare<-realdata[realdata$Study==8,]
###there are 280 patients in this list. 
dataPareRan<-dataPare[which(!is.na(dataPare$VF2)),]
ListRanPID<-as.character(dataPareRan[,2])
###there are 172 inds in this list

	
###########################
##REMOVE INDS WHO ARE NOT FROM THE RANDOM COHORT
###########################
	indstokeep<-vector()
	for (i in 1:length(Mutdata[,1])){
		if (length(which(Mutdata$PATID[i]==ListRanPID))>0   ) {
			indstokeep<-c(indstokeep,i)
		}
	}
#check: there are 172 inds in indstokeep
	Mutdata<-Mutdata[indstokeep,]

#add predicted distribution

write.csv(Mutdata,file="MeanFreqDataRanSubCohort.csv")

Ne=15000; mu=2*10^-6;

plessthanepsilon<-function(eps,Ne,mu,sd){
	theta=2*Ne*mu; alpha_d = 2*Ne*sd;
	return((1/theta)*eps^theta + ((theta - 1 - alpha_d)/(1+theta))*eps^(1+theta))}

pepsilon1minepsilon<-function(eps,Ne,mu,sd){
	theta=2*Ne*mu; alpha_d = 2*Ne*sd;
	integrand <- function(x) {x^(theta-1)*(1-x)^(theta-1)*exp(-alpha_d*x)}
	return(integrate(integrand,eps,1-eps)$value)}

constant<-function(eps,Ne,mu,sd){
	return(plessthanepsilon(eps,Ne,mu,sd)+pepsilon1minepsilon(eps,Ne,mu,sd))}

pfreqfmax<-function(fmax,eps,Ne,mu,sd){
	theta=2*Ne*mu; alpha_d = 2*Ne*sd;
	integrand <- function(x){x^(theta-1)*(1-x)^(theta-1)*exp(-alpha_d*x)}
	return((integrate(integrand,eps,fmax)$value+plessthanepsilon(eps,Ne,mu,sd))/constant(eps,Ne,mu,sd))}
	

###############
#PLOT THE DATA
###############
#		pdf(paste(mut,"Hist_Rawdata_1.pdf",sep=""))
pdf("K103NAACFreqs2014.pdf")

for (breaksize in c(0.0002,0.00002,0.000002)){
Freqs<-Mutdata$K103AACFreq
xvals<-seq(breaksize,20*breaksize,by=breaksize)
yvalsdata<-vector()
for (i in 1:length(xvals)){
	if (i==1)yvalcumprev=0; 
	yvalcum=length(which(Freqs<xvals[i]))/length(Freqs)
	yvalsdata[i]=yvalcum-yvalcumprev
	yvalcumprev=yvalcum	
}
yvals<-vector()
for (i in 1:length(xvals)){
	if (i==1)yvalcumprev=0; 
	yvalcum=pfreqfmax(xvals[i],breaksize,Ne,mu,mu/mean(Freqs))
	yvals[i]=yvalcum-yvalcumprev
	yvalcumprev=yvalcum	
}

ylim=min(1,1.2*max(yvalsdata))
plot(xvals[1:20],yvalsdata[1:20],type="h",lwd=10,ylim=c(0,ylim),col=2,main="K103NAAC",ylab="frequency of mutation",xlab="count / total pats")
points(xvals[1:20],yvals[1:20],t="b",pch=16,cex=0.5)
abline(h=0)
text(xvals[10],ylim*0.8,paste("undetected in",length(which(Freqs==0)),"of",length(Freqs),"patients"))	
text(xvals[10],ylim*0.6,paste("mean freq",round(mean(Freqs),5)))
abline(v=mean(Freqs),col=2,lty=2)
}
dev.off()


pdf("K103NAATFreqs2014.pdf")

for (breaksize in c(0.0002,0.00002,0.000002)){
	Freqs<-Mutdata$K103AATFreq
	xvals<-seq(breaksize,20*breaksize,by=breaksize)
	yvalsdata<-vector()
	for (i in 1:length(xvals)){
		if (i==1)yvalcumprev=0; 
		yvalcum=length(which(Freqs<xvals[i]))/length(Freqs)
		yvalsdata[i]=yvalcum-yvalcumprev
		yvalcumprev=yvalcum	
	}
	yvals<-vector()
	for (i in 1:length(xvals)){
		if (i==1)yvalcumprev=0; 
		yvalcum=pfreqfmax(xvals[i],breaksize,Ne,mu,mu/mean(Freqs))
		yvals[i]=yvalcum-yvalcumprev
		yvalcumprev=yvalcum	
	}
	
	ylim=min(1,1.2*max(yvalsdata))
	plot(xvals[1:20],yvalsdata[1:20],type="h",lwd=10,ylim=c(0,ylim),col=2,main="K103NAAC",ylab="frequency of mutation",xlab="count / total pats")
	points(xvals[1:20],yvals[1:20],t="b",pch=16,cex=0.5)
	abline(h=0)
	text(xvals[10],ylim*0.8,paste("undetected in",length(which(Freqs==0)),"of",length(Freqs),"patients"))	
	text(xvals[10],ylim*0.6,paste("mean freq",round(mean(Freqs),5)))
	abline(v=mean(Freqs),col=2,lty=2)
}
dev.off()

pdf("Y181CFreqs2014.pdf")

for (breaksize in c(0.001,0.0002,0.00002,0.000002)){
	Freqs<-Mutdata$Y181CFreq
	xvals<-seq(breaksize,20*breaksize,by=breaksize)
	yvalsdata<-vector()
	for (i in 1:length(xvals)){
		if (i==1)yvalcumprev=0; 
		yvalcum=length(which(Freqs<xvals[i]))/length(Freqs)
		yvalsdata[i]=yvalcum-yvalcumprev
		yvalcumprev=yvalcum	
	}
	yvals<-vector()
	for (i in 1:length(xvals)){
		if (i==1)yvalcumprev=0; 
		yvalcum=pfreqfmax(xvals[i],breaksize,Ne,mu,mu/mean(Freqs))
		yvals[i]=yvalcum-yvalcumprev
		yvalcumprev=yvalcum	
	}
	
	ylim=min(1,1.2*max(yvalsdata))
	plot(xvals[1:20],yvalsdata[1:20],type="h",lwd=10,ylim=c(0,ylim),col=2,main="K103NAAC",ylab="frequency of mutation",xlab="count / total pats")
	points(xvals[1:20],yvals[1:20],t="b",pch=16,cex=0.5)
	abline(h=0)
	text(xvals[10],ylim*0.8,paste("undetected in",length(which(Freqs==0)),"of",length(Freqs),"patients"))	
	text(xvals[10],ylim*0.6,paste("mean freq",round(mean(Freqs),5)))
	abline(v=mean(Freqs),col=2,lty=2)
}
dev.off()


