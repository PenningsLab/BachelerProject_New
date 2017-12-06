#cvbcxc
#pdf("Sims_estimating_mu_over_s_DIFF_s_values.pdf")
system("./Code_and_shellscript/make_HIV1site")    #compile the code 
library(ggplot2)

#DataOverview<-data.frame("N"=Ne,"mu"=0,"cost"=0,"num_runs"=0,"datapointsperrun"=0,"equiPi"=0,"VarPi"=0,"expectedPi"=0,"t_half"=0)

Ne = 10000 #currently Ne cannot be changed in the sims 
nmax = 200 #(max num of patients, should be around 200)
NUMRUNS=1; numoutputs=nmax
if (FALSE){ #create data through simulations
    thetas = c(0.1,1, 10)
    for (theta in thetas){
        mu = theta / Ne
        numsites = 300
        costs = sort(round(runif(numsites),2))
        seed =1
        for (site in 1:numsites){
            print(paste("site",site))
            cost = costs[site]
            #make script 
            x<-"#!/bin/bash"
            x<-c(x,paste("mu=",mu,sep=""))
            outputfrequency=min(c(2*Ne,ceiling(5/cost)))
            #outputfrequency=2*Ne
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
                 paste("\" | ./Code_and_shellscript/HIVevolution_HIV1site >./Data/Data_s_", site, "_T_", theta, "_cost_", cost,".txt",sep=""), 
                 "done")
            #x<-c(x,"done")
            write(x,file="./Code_and_shellscript/tempscript.sh")
            system("chmod 775 ./Code_and_shellscript/tempscript.sh")
            #Run tempscript.sh
            system("./Code_and_shellscript/tempscript.sh")
        }
    }
}

##########Read data

#ListFiles<-list.files("../Data/")
#ListCosts<-rep(0,length(ListFiles))
#ListMeanFreqs<-rep(0,length(ListFiles))
#for (i in 1:length(ListFiles)){
#    ListCosts[i]<-as.numeric(substr(ListFiles[i], regexpr("cost",ListFiles[i])[1]+5, regexpr(".txt",ListFiles[i])[1]+-1))
    #read the file
#    ListMeanFreqs[i]<-mean(read.csv(paste("../Data/",ListFiles[i],sep=""),sep = "\t")$freq)
#}
#plot(ListCosts,mu/ListMeanFreqs)
#print(cor.test(ListCosts,mu/ListMeanFreqs))
#> plot(ListCosts,mu/ListMeanFreqs,ylim=c(0,1))
#> plot(ListCosts,mu/ListMeanFreqs,ylim=c(0,2))
#> plot(ListCosts,ListMeanFreqs,ylim=c(0,1))
#> plot(ListCosts,ListMeanFreqs,log="y")

###Hm, it looks like the very deleterious mutations have way too low frequencies. Fixed. 
###I need to check how i did the simulations SWAPPED ORDER OF MUTATION AND SELECTION. Now I sample after mutation, but before selection. 
###Also create a github repository. DONE

#####function to go from mu/meanfreq to estimated s (including case where meanfreq = 0)
EstimatedS <- function(mu, meanfreq){
    if (meanfreq == 0) return (1)
    else return (min(c(mu/meanfreq,1)))
}

#####Looking at smaller number of patients. 
thetas<-c(0.1,1,10)
ListFiles<-list.files("../Data/")
NumPatsList<-c(1,2,3,5,7,seq(10,200,by=20))   

CorData<-data.frame("Numpats"=as.numeric(),"Theta"=as.numeric(),"CorCoeff"=as.numeric())

#for (theta in thetas){
theta = 1
    mu = theta / Ne
    thetapattern=paste("T_",theta,"_",sep="")
    ListFilesTheta<-ListFiles[grep(thetapattern,ListFiles)]
    for (Numpats in NumPatsList){
        ListCosts<-rep(0,length(ListFilesTheta))
        ListMeanFreqs<-rep(0,length(ListFilesTheta))
        Pats<-sample(200,Numpats)
        for (i in 1:length(ListFilesTheta)){
            ListCosts[i]<-as.numeric(substr(ListFilesTheta[i], regexpr("cost",ListFilesTheta[i])[1]+5, regexpr(".txt",ListFilesTheta[i])[1]+-1))
            #read the file
            ListMeanFreqs[i]<-mean(read.csv(paste("../Data/",ListFilesTheta[i],sep=""),sep = "\t")$freq[Pats])
        }
        ListEstimatedCosts<-sapply(ListMeanFreqs,function(x) EstimatedS(mu,x))
        plot(ListCosts,ListEstimatedCosts,main=paste("Numpats",Numpats),ylim=c(0,1),pch=16, col=2) #alpha=0.5,col=2)
        CorData[nrow(CorData)+1,]<-c(Numpats,theta,round(cor.test(ListCosts,sapply(ListMeanFreqs,function(x) EstimatedS(mu,x)))$estimate,3))
        #text(0.6,0.2,paste("cor = ", Cor))
        print(paste("theta",theta,"numpats",Numpats,"cor",CorData$CorCoeff[nrow(CorData)]))
        #CorList<-c(CorList,Cor)
    }
#}

CorData$Theta<-as.factor(CorData$Theta)

p <- ggplot(CorData, aes(Numpats, CorCoeff, colour = Theta))
p + geom_line() +geom_point() +ylim(0,1)

#Ran with diff thetas
###What if freq = 0 , then m/0 give NaN. Should be cost = 1. DONE

###Next: what if sequencing less precise! 
thetas<-c(0.1,1,10)
ListFiles<-list.files("../Data/")
NumPatsList<-c(1,2,3,5,7,seq(10,200,by=20))   
SampleDepths = c(20,100,1000,10000)

CorDataSampleDepth<-data.frame("Numpats"=as.numeric(),"Theta"=as.numeric(),"CorCoeff"=as.numeric(), "SampleDepth" = as.numeric())

#for (theta in thetas){
theta=1
    mu = theta / Ne
    thetapattern=paste("T_",theta,"_",sep="")
    ListFilesTheta<-ListFiles[grep(thetapattern,ListFiles)]
    ListCosts<-rep(0,length(ListFilesTheta))    
    for (Numpats in NumPatsList){
        Pats<-sample(200,Numpats)
        ListMeanSampleFreqs<-rep(0,length(ListFilesTheta)) # a different list of 300 meansamplefreqs for each unique set of pats
        for (SampleDepth in SampleDepths){
            print(paste("SampleDepth",SampleDepth))
            for (i in 1:length(ListFilesTheta)){
                #for each file get the cost of the mutation from the filename
                ListCosts[i]<-as.numeric(substr(ListFilesTheta[i], regexpr("cost",ListFilesTheta[i])[1]+5, regexpr(".txt",ListFilesTheta[i])[1]+-1))
                #read the file and get the population frequencies for the relevant Pats
                PopFreqs<-read.csv(paste("../Data/",ListFilesTheta[i],sep=""),sep = "\t")$freq[Pats]
                #get the samplefreqs using the PopFreqs and a sampledepth
                SampleFreqs<-rbinom(length(PopFreqs),SampleDepth,PopFreqs)/SampleDepth
                ListMeanSampleFreqs[i]<-mean(SampleFreqs)
                }
            if (length( which(ListMeanSampleFreqs>0) )==0) Cor = NA
            else Cor = round(cor.test(ListCosts,sapply(ListMeanSampleFreqs,function(x) EstimatedS(mu,x)))$estimate,3)
            CorDataSampleDepth[nrow(CorDataSampleDepth)+1,]<-c(Numpats,theta,Cor,SampleDepth)
            print(paste("theta",theta,"numpats",Numpats,"cor",Cor,"SampleDepth",SampleDepth))
        }
    }
#}

CorDataSampleDepth$Theta<-as.factor(CorDataSampleDepth$Theta)
CorDataSampleDepth$SampleDepth<-as.factor(CorDataSampleDepth$SampleDepth)
i=0

i=i+1
SampleDepth= SampleDepths[i]
Data<-CorDataSampleDepth[CorDataSampleDepth$SampleDepth == SampleDepth,]
p <- ggplot(Data, aes(Numpats, CorCoeff, colour = Theta))
p + geom_point() + geom_line() + ggtitle(paste("Sample Depth = ",SampleDepth) )+ ylim(0, 1)

Data<-CorDataSampleDepth[CorDataSampleDepth$Theta == 1,]
p <- ggplot(Data, aes(Numpats, CorCoeff, colour = SampleDepth))
p + geom_point(aes(size=2)) + geom_line() + ggtitle(paste("theta = ",1) )+ ylim(0, 1)
p + geom_point() + geom_line() + ggtitle(paste("theta = ",1) )+ ylim(0, 1)



#####Old code from 2014, I don't think I need this. 
if (FALSE){
#READ FREQS FILE
read.csv(paste("../Data/Link",cost,"_",seed,".txt",sep=""),sep="\t",header=TRUE)->simdata
#in stead of the real frequencies, lets assume we have a sample from each patient of, say, 100, seqs. 
plot(c(0,0),col=0,xlim=c(-0.5,3.5),ylim=c(-1*mu/cost,8*mu/cost),xlab="Num Patients",xaxt="n",ylab="95% of observed average freq of allele",main=paste("Ne",Ne,", mu",mu,", Theta",2*Ne*mu ,", cost",cost))
axis(1, at=log10(samplesizes), labels=samplesizes)
abline(h=mu/cost,lty=2)

diffnumseqs <- c(50,100,200)
for (num_seqs_per_patient in diffnumseqs){
    for (i in 1:length(simdata$freq)){
        simdata$est_freq[i]<-rbinom(1,num_seqs_per_patient,simdata$freq[i])/num_seqs_per_patient}
    
    #system(paste("rm ","../Data/Link",cost,"_",seed,".txt",sep=""))
    samplesizes<-c(1,3,10,30,100,300,1000)
    for (num_patients in samplesizes){
        list_averages<-vector()
        for (i in 1:1000){ 
            list_averages<-c(list_averages,mean(sample(simdata$est_freq,num_patients)))}
        co=which(diffnumseqs==num_seqs_per_patient)
        X=(which(diffnumseqs==num_seqs_per_patient)-2)*0.1
        print(paste(num_seqs_per_patient,num_patients,log10(num_patients)-0.02+X))
        rect(log10(num_patients)-0.02+X,sort(list_averages)[25],log10(num_patients)+0.02 +X,sort(list_averages)[975],col=co)
        #	text(log10(num_patients)+X,sort(list_averages)[975]+0.15*mu/cost+(X+0.2)/200,paste(round(sort(list_averages)[25]/(mu/cost),2),"-",round(sort(list_averages)[975]/(mu/cost),2)),cex=0.6)
        
        text(log10(num_patients)+X,(-1+co/8)*mu/cost,paste(round(sort(list_averages)[25]/(mu/cost),2),"-",round(sort(list_averages)[975]/(mu/cost),2)),cex=0.5)	
    }
    text(log10(num_patients),6.5*mu/cost,"bars and",cex=0.8)
    text(log10(num_patients),6*mu/cost,"numbers indicate range",cex=0.8)
    text(log10(num_patients),5.5*mu/cost,"of 95% of estimates",cex=0.8)
    text(log10(num_patients),5.*mu/cost,"black:50, red: 100, ",cex=0.8)
    text(log10(num_patients),4.5*mu/cost,"green: 200 sequences/pat",cex=0.8)
    text(-0.3,1.1*mu/cost,"mu/s",cex=0.8)
    
}

dev.off()
}
#rbinom(n, size, prob)

