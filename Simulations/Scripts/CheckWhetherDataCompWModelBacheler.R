#July 2017 , July 11th. 
#Dec 2017

#My goal is to see if the observed Single Site Freq Spec (Bacheler data) fits with predicted for fitness cost paper. 

setwd("~/Documents/Git/bachelerProject/SimulationsEstimatingSelCoeffSims/")
#system("./Code_and_shellscript/make_HIV1site")    #compile the code
#system("rm ./SimData/Data*") #remove old sim files

#Read the data
#We read a file with WT 1 as threshold
read.table("../Output/freqPatTs_Bacheler_Threshold1.csv",sep=",",header=TRUE,row.names=1)->freqPatTsFilter
#calculate mean frequencies
MeanFreq<-apply(freqPatTsFilter, 2 , mean, na.rm=TRUE)
read.csv("../Output/OverviewSelCoeff_BachelerFilter.csv")->OverviewDFilter

######
#First, let's get the sample sizes for Bacheler data and store them in NumseqsBach
#Read the correct fastafiles.
library(seqinr)
listfastafiles<-list.files("../Data/BachelerFiles//FASTAfiles")
NumseqsBach<-c()
for (i in 1:length(listfastafiles)){ #for each fastafile
    patfasta<-read.fasta(paste("../Data/BachelerFiles/FASTAfiles/",substr(listfastafiles[i],1,6),".fasta",sep="")) #read the file
    if (length(patfasta)>=5) NumseqsBach<-c(NumseqsBach,length(patfasta))#remove patients with fewer than 5 sequences
}

pdf ("comparisonDataSims_New_Dec2017.pdf") #create a pdf file with freq distr for real and sim data
for (Ne in c(5000)){

    listPvalues<-c()
    #Choose a site 
    for (site in OverviewDFilter$num[OverviewDFilter$TypeOfSite%in%c("nonsyn","syn")]){
    #for (site in 172:174){
        system("rm ./SimData/Data*")
        
        print(site)
        realFreqs<-freqPatTsFilter[,site]
        
        #Mut rates and sel coefficients
        mu=round(OverviewDFilter$TSmutrate[site],5)
        print(mu)
        cost=round(OverviewDFilter$EstSelCoeff[site],5)
        print(cost)
        
        #Use code from 
        #SimulationsEstimatingSelCoeffSims
        #To create simulated data. 
        numoutputs = length(NumseqsBach) #(max num of patients, should be around 200)
        theta = mu*Ne
        
        if (TRUE){
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
            
            if (Ne == 1000)sentence = paste("\" | ./Code_and_shellscript/HIVevolution_HIV1site1000 >./SimData/Data_T_", theta, "_cost_", cost,".txt",sep="")
            if (Ne == 5000)sentence = paste("\" | ./Code_and_shellscript/HIVevolution_HIV1site5000 >./SimData/Data_T_", theta, "_cost_", cost,".txt",sep="")
            if (Ne == 10000)sentence = paste("\" | ./Code_and_shellscript/HIVevolution_HIV1site10000 >./SimData/Data_T_", theta, "_cost_", cost,".txt",sep="")
            if (Ne == 50000)sentence = paste("\" | ./Code_and_shellscript/HIVevolution_HIV1site50000 >./SimData/Data_T_", theta, "_cost_", cost,".txt",sep="")
            if (Ne == 100000)sentence = paste("\" | ./Code_and_shellscript/HIVevolution_HIV1site100000 >./SimData/Data_T_", theta, "_cost_", cost,".txt",sep="")
            
            x<-c(x,"do",
                 "echo \"", "$seed", "$mu", "$cost",
                 "$output_every_Xgen", "$numgen_inN", "$start_output", sentence, 
                 "done")
            
            write(x,file="./Code_and_shellscript/tempscript.sh")
            system("chmod 775 ./Code_and_shellscript/tempscript.sh")
            system("./Code_and_shellscript/tempscript.sh")     #Run tempscript.sh
        }
        #Read the data
        filename=paste("./SimData/Data_T_", theta, "_cost_", cost,".txt",sep="")
        read.csv(filename,sep="\t")$freq->Freqs
        system("rm ./SimData/Data*")
        
        #Get sample Freqs (sample from pop freqs)
        #use the real n from the Bacheler data. 
        numoutputs<-length(which(!is.na(realFreqs)))
        simFreqs<-rbinom(numoutputs,NumseqsBach[!is.na(realFreqs)],Freqs)/NumseqsBach[!is.na(realFreqs)]
        realFreqs<-realFreqs[!is.na(realFreqs)]
        
        #Plot
        br<-seq(-10^-7,1.01,by=0.01)
        br2<-seq(0,1,by=0.02)
        xlims=c(0,0.01+max(c(simFreqs,realFreqs),na.rm=TRUE))
        ylims = c(0,max(length(which(realFreqs==0)),length(which(simFreqs==0))))
        #ylims = c(0,5)
        
        par(mfrow=c(2,1))
        maxyheight=50
        
        hist(rep(0,maxyheight),breaks=br,xlim=xlims, #ylim = ylims, 
             ylab="", xlab = "", main=paste("Single-SFS, site",site),
             xaxt="n", yaxt="n",col=OverviewDFilter$color[site])
        hist(c(rep(0,min(maxyheight-10,length(which(realFreqs<0.01)))),realFreqs[which(realFreqs>=0.01)]), breaks=br, col=OverviewDFilter$color[site], add=TRUE)
        axis(2,labels = c(10,20,30,40,max(maxyheight,length(which(realFreqs<0.01)))), 
             at = c(10,20,30,40,maxyheight), las=1)
        if (length(which(realFreqs<0.01))>=maxyheight){
            axis.break(axis=2,breakpos=maxyheight-20,bgcol="white",breakcol="black",style="slash",brw=0.02)
            rect(-0.001, maxyheight-12, 0.0101, maxyheight-8,col="white",border="white")
        }else{axis(2,labels = maxyheight-10,at=maxyheight-10,las=1)}
        
        axis(1, at=br2,labels=br2, col.axis="black", las=2)
        #hist(rep(0,length(which(realFreqs==0))),breaks=br, add = T, col=1)
        mtext("Observed frequency", side=1, line=4, cex.lab=1,las=1, col="black")
        mtext("Number of patients", side=2, line=2.5, cex.lab=1,las=3, col="black")
        abline(v=mean(realFreqs,na.rm=TRUE),col=2)
        
        hist(rep(0,maxyheight),breaks=br,xlim=xlims, #ylim = ylims, 
             ylab="", xlab = "", main=paste("Simulated data, cost",round(cost,4)),
             xaxt="n", yaxt="n",col="pink")
        hist(c(rep(0,min(maxyheight-10,length(which(simFreqs<0.01)))),simFreqs[which(simFreqs>=0.01)]), breaks=br, col="pink", add=TRUE)
        axis(2,labels = c(10,20,30,40,max(maxyheight,length(which(simFreqs<0.01)))), 
             at = c(10,20,30,40,maxyheight), las=1)
        if (length(which(simFreqs<0.01))>=maxyheight){
            axis.break(axis=2,breakpos=maxyheight-20,bgcol="white",breakcol="black",style="slash",brw=0.02)
            rect(-0.001, maxyheight-12, 0.0101, maxyheight-8,col="white",border="white")
        }else{axis(2,labels = maxyheight-10,at=maxyheight-10,las=1)}
        
        text(x=mean(xlims),y=maxyheight/2,paste("pvalue=",round(wilcox.test(simFreqs,realFreqs)[[3]],3)))
        axis(1, at=br2,labels=br2, col.axis="black", las=2)
        #hist(rep(0,length(which(simFreqs==0))),breaks=br, add = T, col=2)
        mtext("Observed frequency", side=1, line=4, cex.lab=1,las=1, col="black")
        mtext("Number of patients", side=2, line=2.5, cex.lab=1,las=3, col="black")
        abline(v=mean(simFreqs),col=2)
        
        listPvalues<-c(listPvalues,wilcox.test(simFreqs,realFreqs)[[3]])
        #listPvalues<-c(listPvalues,ks.test(simFreqs, realFreqs)[[2]])
        #ks.test(simFreqs, realFreqs)[[2]]
        
        #This was done to keep the data for positions 172, 173 and 174 for the first figure in the paper. 
        if (site ==172) {
            Freqs172<-simFreqs; #write.csv(Freqs172,"../Output/SimFreqs172.csv")
            }
        if (site ==173) {
            Freqs173<-simFreqs; #write.csv(Freqs173,"../Output/SimFreqs173.csv")
            }
        if (site ==174) {
            Freqs174<-simFreqs; #write.csv(Freqs174,"../Output/SimFreqs174.csv")
            }
    }
    
    qqnorm(listPvalues)
    qqline(listPvalues)
    print(paste("N",Ne,"num p values",round(length(which(listPvalues<.05))/length(listPvalues),3)))
    
    if (Ne == 1000) listPvalues1000<-listPvalues 
    if (Ne == 5000) listPvalues5000<-listPvalues 
    if (Ne == 10000) listPvalues10000<-listPvalues 
    if (Ne == 50000) listPvalues50000<-listPvalues 
    if (Ne == 100000) listPvalues100000<-listPvalues 
    
}    

dev.off()

for (Ne in c(5000)){
    
    if (Ne == 1000) listPvalues1000->listPvalues 
    if (Ne == 5000) listPvalues5000->listPvalues 
    if (Ne == 10000) listPvalues10000->listPvalues 
    if (Ne == 50000) listPvalues50000->listPvalues 
    if (Ne == 100000) listPvalues100000->listPvalues 
    
    Freqs<-OverviewDFilter$MeanFreq[OverviewDFilter$num[OverviewDFilter$TypeOfSite%in%c("nonsyn","syn")]]
    listPvalues[which(Freqs<median(Freqs))]
    
    #how many sites have sign different freq distributions?
    
    print(paste("N",Ne,"fraction p values <5% all frequencies",
                round(length(which(listPvalues<.05))/length(which(!is.na(listPvalues))),3)
    ))
}
