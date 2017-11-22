#July 2017 

#My goal is to see if the observed Single Site Freq Spec fits with predicted for fitness cost paper. 

setwd("~/Documents/Git/bachelerProject/ParedesDataA5095")
#system("./Code_and_shellscript/make_HIV1site")    #compile the code
system("rm ./SimData/Data*")

Data<-read.csv("MeanFreqDataRanSubCohort.csv")

Freqs103AAC<-Data$K103AACFreq[!is.na(Data$K103AACFreq)]
Freqs103AAT<-Data$K103AATFreq[!is.na(Data$K103AATFreq)]
#Freqs103AAT<-sort(Freqs103AAT)[1:(length(Freqs103AAT)-2)]

#Mut rates and sel coefficients
read.csv("../Data/HIVMutRates/HIVMutRates.csv")->mutrates
#AC mut rate
mutrates$Probability[mutrates$Nucleotide.substitution=="AC"]

mu_AC=mutrates$Probability[mutrates$Nucleotide.substitution=="AC"]
mu_AT=mutrates$Probability[mutrates$Nucleotide.substitution=="AU"] #divide by 5?
#mu_AC=9*10^-7
mu_AC_Zan=mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="AC"]
mu_AT_Zan=mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="AU"]/5


#mu_AT=(7*10^-7)
#mu_AT=(7*10^-7)/5
s_AAC<-mu_AC/mean(Freqs103AAC)
s_AAT<-mu_AT/mean(Freqs103AAT)
#(cost of AAC mutation 0.7 %)

#Use code from 
#SimulationsEstimatingSelCoeffSims
#To create simulated data. 

Ne = 300000 #currently Ne cannot be changed in the sims, currently 300000 Nov 2017
nmax = 200 #(max num of patients, should be around 200)
theta_AC = mu_AC*Ne
theta_AT = mu_AT*Ne
NUMRUNS=1; numoutputs=nmax

#for AAC
if (TRUE){ #create data through simulations
    numsites = 1
    cost = s_AAC; theta = theta_AC
    seed =100
    #make script 
    x<-"#!/bin/bash"
    x<-c(x,paste("mu=",mu_AC,sep=""))
    outputfrequency=min(c(2*Ne,ceiling(5/cost)))
    x<-c(x,paste("output_every_Xgen=",outputfrequency,sep=""))
    x<-c(x,paste("numgen_inN=",(numoutputs+2)*outputfrequency/Ne,sep=""))
    x<-c(x,paste("start_output=",2*outputfrequency/Ne,sep=""))
    x<-c(x,paste("cost=",cost,sep=""))
    #x<-c(x,paste("for cost in", paste(costs,collapse = " ")))
    #x<-c(x,"do")
    x<-c(x,paste("for seed in",seed))
    x<-c(x,"do",
         "echo \"", "$seed", "$mu", "$cost",
         "$output_every_Xgen", "$numgen_inN", "$start_output",
         paste("\" | ./Code_and_shellscript/HIVevolution_HIV1site >./SimData/Data_T_", theta, "_cost_", cost,".txt",sep=""), 
         "done")
    #x<-c(x,"done")
    write(x,file="./Code_and_shellscript/tempscript.sh")
    system("chmod 775 ./Code_and_shellscript/tempscript.sh")
    #Run tempscript.sh
    system("./Code_and_shellscript/tempscript.sh")
}
#Read the data
list.files("SimData/")->listfiles
filename=paste("SimData/",listfiles[1],sep="")
read.csv(filename,sep="\t")$freq->PopFreqs103AAC_sim
system("rm ./SimData/Data*")

#for AAT
if (TRUE){ #create data through simulations
    numsites = 1
    cost = s_AAT; theta = theta_AT
    seed =100
    #make script 
    x<-"#!/bin/bash"
    x<-c(x,paste("mu=",mu_AT,sep=""))
    outputfrequency=min(c(2*Ne,ceiling(5/cost)))
    x<-c(x,paste("output_every_Xgen=",outputfrequency,sep=""))
    x<-c(x,paste("numgen_inN=",(numoutputs+2)*outputfrequency/Ne,sep=""))
    x<-c(x,paste("start_output=",2*outputfrequency/Ne,sep=""))
    x<-c(x,paste("cost=",cost,sep=""))
    #x<-c(x,paste("for cost in", paste(costs,collapse = " ")))
    #x<-c(x,"do")
    x<-c(x,paste("for seed in",seed))
    x<-c(x,"do",
         "echo \"", "$seed", "$mu", "$cost",
         "$output_every_Xgen", "$numgen_inN", "$start_output",
         paste("\" | ./Code_and_shellscript/HIVevolution_HIV1site >./SimData/Data_T_", theta, "_cost_", cost,".txt",sep=""), 
         "done")
    #x<-c(x,"done")
    write(x,file="./Code_and_shellscript/tempscript.sh")
    system("chmod 775 ./Code_and_shellscript/tempscript.sh")
    #Run tempscript.sh
    system("./Code_and_shellscript/tempscript.sh")
}
#Read the data
list.files("SimData/")->listfiles
filename=paste("SimData/",listfiles[1],sep="")
read.csv(filename,sep="\t")$freq->PopFreqs103AAT_sim


#Get sample Freqs
numparticles=50000
Freqs103AAC_sim<-rbinom(n=200,size=numparticles,prob=PopFreqs103AAC_sim)/numparticles
Freqs103AAT_sim<-rbinom(n=200,size=numparticles,PopFreqs103AAT_sim)/numparticles
#Plot
br<-seq(-0.00001+10^-7,0.1,by=0.00001)

if (TRUE){
xlimset=c(-0.00001+10^-7,0.0008)    
    
pdf("ParedesData_Sims.pdf",width = 10, height = 10)
par(mfrow=c(2,2))
#png("Freqs103AAC.png")
hist(Freqs103AAC,xlim=xlimset,breaks=br, 
     ylab="", xlab = "", main="SFS RT103AAC",
     xaxt="n", col="blue"
)
axis(1, at=c(round(br[2],5)-0.000005,round(br[3:16],5)),labels=round(br[2:16],5), col.axis="black", las=2)
hist(rep(0,length(which(Freqs103AAC==0))),breaks=br, add = T, col=1)
mtext("Observed frequency", side=1, line=4, cex.lab=1,las=1, col="black")
mtext("Number of patients", side=2, line=2.5, cex.lab=1,las=3, col="black")
abline(v=mean(Freqs103AAC),col=2)

#dev.off()

hist(Freqs103AAC_sim,xlim=xlimset,breaks=br, 
     ylab="", xlab = "", main="Simul.SFS RT103AAC",
     xaxt="n", col="blue"
)
axis(1, at=c(round(br[2],5)-0.000005,round(br[3:16],5)),labels=round(br[2:16],5), col.axis="black", las=2)
hist(rep(0,length(which(PopFreqs103AAC_sim==0))),breaks=br, add = T, col=1)
mtext("Observed frequency", side=1, line=4, cex.lab=1,las=1, col="black")
mtext("Number of patients", side=2, line=2.5, cex.lab=1,las=3, col="black")
abline(v=mean(Freqs103AAC_sim),col=2)

hist(Freqs103AAT,xlim=xlimset,breaks=br, 
     ylab="", xlab = "", main="SFS RT103AAT",
     xaxt="n", col="blue"
)
axis(1, at=c(round(br[2],5)-0.000005,round(br[3:16],5)),labels=round(br[2:16],5), col.axis="black", las=2)
hist(rep(0,length(which(Freqs103AAT==0))),breaks=br, add = T, col=1)
mtext("Observed frequency", side=1, line=4, cex.lab=1,las=1, col="black")
mtext("Number of patients", side=2, line=2.5, cex.lab=1,las=3, col="black")
abline(v=mean(Freqs103AAT),col=2)

hist(Freqs103AAT_sim,xlim=xlimset,breaks=br, 
     ylab="", xlab = "", main="Simul.SFS RT103AAT",
     xaxt="n", col="blue"
)
axis(1, at=c(round(br[2],5)-0.000005,round(br[3:16],5)),labels=round(br[2:16],5), col.axis="black", las=2)
hist(rep(0,length(which(PopFreqs103AAT_sim==0))),breaks=br, add = T, col=1)
mtext("Observed frequency", side=1, line=4, cex.lab=1,las=1, col="black")
mtext("Number of patients", side=2, line=2.5, cex.lab=1,las=3, col="black")
abline(v=mean(Freqs103AAT_sim),col=2)


dev.off()
}
