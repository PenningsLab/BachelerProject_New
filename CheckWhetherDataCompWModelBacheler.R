#July 2017 , July 11th. 

#My goal is to see if the observed Single Site Freq Spec (Bacheler data) fits with predicted for fitness cost paper. 

setwd("~/Documents/Git/bachelerProject/SimulationsEstimatingSelCoeffSims/")
system("./Code_and_shellscript/make_HIV1site")    #compile the code
#system("rm ./SimData/Data*") #remove old sim files

#Read the data
#We read a file with WT 0 as threshold, meaning no threshold July 12, now with filter)
read.table("../Output/freqPatTs_Bacheler_filter.csv",sep=",",header=TRUE,row.names=1)->freqPatTsFilter
#calculate mean frequencies
MeanFreq<-apply(freqPatTsFilter, 2 , mean, na.rm=TRUE)
read.csv("../Output/OverviewSelCoeff_BachelerFilter.csv")->OverviewDF

######
#First, let's get the sample sizes for Bacheler data and store them in Numseqs
#Read the correct fastafiles.
library(seqinr)
listfastafiles<-list.files("../Data/BachelerFiles//FASTAfiles")
NumseqsBach<-c()
for (i in 1:length(listfastafiles)){ #for each fastafile
    patfasta<-read.fasta(paste("../Data/BachelerFiles/FASTAfiles/",substr(listfastafiles[i],1,6),".fasta",sep="")) #read the file
    if (length(patfasta)>=5) NumseqsBach<-c(NumseqsBach,length(patfasta))#remove patients with fewer than 5 sequences
}

pdf ("comparisonDataSims.pdf")

#Choose a site 
for (site in 15:110){
realFreqs<-freqPatTsFilter[,site]

#Mut rates and sel coefficients
mu=OverviewDF$TSmutrate[site]
cost=OverviewDF$EstSelCoeff[site]

#Use code from 
#SimulationsEstimatingSelCoeffSims
#To create simulated data. 

Ne = 100000 #currently Ne cannot be changed in the sims, has to be 10000 # 07/11/17 I changed it to 100000
numoutputs = length(NumseqsBach) #(max num of patients, should be around 200)
theta = mu*Ne

#create data through simulations
seed =100
#make script 
x<-"#!/bin/bash"
x<-c(x,paste("mu=",mu,sep=""))
outputfrequency=min(c(2*Ne,ceiling(5/cost)))
x<-c(x,paste("output_every_Xgen=",outputfrequency,sep=""))
x<-c(x,paste("numgen_inN=",(numoutputs+2)*outputfrequency/Ne,sep=""))
x<-c(x,paste("start_output=",2*outputfrequency/Ne,sep=""))
x<-c(x,paste("cost=",cost,sep=""))
x<-c(x,paste("for seed in",seed))
x<-c(x,"do",
     "echo \"", "$seed", "$mu", "$cost",
     "$output_every_Xgen", "$numgen_inN", "$start_output",
     paste("\" | ./Code_and_shellscript/HIVevolution_HIV1site >./SimData/Data_T_", theta, "_cost_", cost,".txt",sep=""), 
     "done")
write(x,file="./Code_and_shellscript/tempscript.sh")
system("chmod 775 ./Code_and_shellscript/tempscript.sh")
system("./Code_and_shellscript/tempscript.sh")     #Run tempscript.sh

#Read the data
list.files("SimData/")->listfiles
filename=paste("SimData/",listfiles[1],sep="")
read.csv(filename,sep="\t")$freq->Freqs
system("rm ./SimData/Data*")

#Get sample Freqs
#use the real n from the Bacheler data. 
numoutputs<-length(which(!is.na(realFreqs)))
Freqs<-rbinom(numoutputs,NumseqsBach[!is.na(realFreqs)],Freqs)/NumseqsBach[!is.na(realFreqs)]

#Plot
br<-seq(-10^-7,1.01,by=0.01)
br2<-seq(0,1,by=0.05)
xlims=c(0,0.01+max(c(Freqs,realFreqs),na.rm=TRUE))
ylims = c(0,max(length(which(realFreqs==0)),length(which(Freqs==0))))
#ylims = c(0,5)

   
par(mfrow=c(2,1))

hist(realFreqs,breaks=br,xlim=xlims, ylim = ylims,
     ylab="", xlab = "", main=paste("Single-SFS, site",site),
     xaxt="n", col="blue"
)
axis(1, at=br2,labels=br2, col.axis="black", las=2)
hist(rep(0,length(which(realFreqs==0))),breaks=br, add = T, col=1)
mtext("Observed frequency", side=1, line=4, cex.lab=1,las=1, col="black")
mtext("Number of patients", side=2, line=2.5, cex.lab=1,las=3, col="black")
abline(v=mean(realFreqs,na.rm=TRUE),col=2)

hist(Freqs,breaks=br,xlim=xlims,ylim = ylims,
     ylab="", xlab = "", main=paste("Simulated Data, cost",round(cost,5)),
     xaxt="n", col="blue"
)
axis(1, at=br2,labels=br2, col.axis="black", las=2)
hist(rep(0,length(which(Freqs==0))),breaks=br, add = T, col=1)
mtext("Observed frequency", side=1, line=4, cex.lab=1,las=1, col="black")
mtext("Number of patients", side=2, line=2.5, cex.lab=1,las=3, col="black")
abline(v=mean(Freqs),col=2)
}

dev.off()
