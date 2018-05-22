source("Rscripts/baseRscript.R")
source("Rscripts/RResistanceMutations.r")
read.csv("Data/HIVMutRates/HIVMutRates.csv")->mutrates

AllMuts<-rbind(RTImuts,PImuts)
AllMuts$wtcodon<-""
AllMuts$rescodon<-""
AllMuts$easiestMut<-""
AllMuts$mutrate<-0
AllMuts$type<-c(rep("NN",9),rep("NR",11),"NA","184", rep("PI",12))
#WTCodonList<-list()

#Get WT AA at all positions
for (i in 1:nrow(AllMuts)){
    AApos = AllMuts$pos[i]+99
    if (AllMuts$type[i]=="PI"){AApos = AllMuts$pos[i]}
    AllMuts$wt[i]<-translate(consensusB)[AApos]
    AllMuts$wtcodon[i]<-paste(consensusB[(3*AApos-2):(3*AApos)],collapse = "")
}

#get all possible resistance codons
allcodons<-c(); allAA<-c()
for (i in c("a","c","g","t")){
    for (j in c("a","c","g","t")){
        for (k in c("a","c","g","t")){
            codon=paste(i,j,k,collapse ="",sep="")
            allcodons<-c(allcodons,codon)
            allAA<-c(allAA,translate(c(i,j,k)))
        }
    }
}

#Which are the resistant codons?
for (i in 1:nrow(AllMuts)){
resAA<-as.vector(strsplit(AllMuts$mut[i],split=""))
rescodons<-allcodons[which(allAA%in%resAA[[1]])]

#from the WT, all 1-step mutants (i include the wt because  easier to program)
onestepmuts<-c();wtcodon=AllMuts$wtcodon[i]
resmutrates=c()
for (j in 1:3){
    if (substr(wtcodon, j, j) =="a"){muts<-c("g","t","c")}
    if (substr(wtcodon, j, j) =="c"){muts<-c("t","a","g")}
    if (substr(wtcodon, j, j) =="g"){muts<-c("a","c","t")}
    if (substr(wtcodon, j, j) =="t"){muts<-c("c","a","g")}
    for (l in muts){
        mutcodon = wtcodon
        substr(mutcodon, j, j) <- l
        onestepmuts<-c(onestepmuts, mutcodon)
        x=toupper(substr(wtcodon, j, j));if (x=="T") x="U"
        y=toupper(l); if (y=="T") y="U"
        mutationtype=paste(x,y,collapse = "",sep="")
        resmutrates=c(resmutrates,mutrates$Probability[which(mutrates$Nucleotide.substitution==mutationtype)])
}}

#get the highest mut rate for a one step mut that leads to a res codon!
whichonesteps<-which(onestepmuts%in%rescodons)
onestepmuts<-onestepmuts[whichonesteps]
resmutrates<-resmutrates[whichonesteps]

if (length(whichonesteps)>0){
m=resmutrates[which.max(resmutrates)]
o=onestepmuts[which.max(resmutrates)]
AllMuts$mutrate[i]<-m
AllMuts$rescodon[i]<-o
}
}


CummMutRateNNRTI = 1- prod(1-AllMuts$mutrate[which(AllMuts$type=="NN")])
CummMutRateNRTI = 1- prod(1-AllMuts$mutrate[which(AllMuts$type=="NR")])
CummMutRate184 = 1- prod(1-AllMuts$mutrate[which(AllMuts$type=="184")])
CummMutRatePI = 1- prod(1-AllMuts$mutrate[which(AllMuts$type=="PI")])

#For NNRTI treatments
ProbNNRTIFirst<-CummMutRateNNRTI/(CummMutRateNNRTI+CummMutRateNRTI+CummMutRate184)
ProbNRTIFirst<-CummMutRateNRTI/(CummMutRateNNRTI+CummMutRateNRTI+CummMutRate184)
Prob184First<-CummMutRate184/(CummMutRateNNRTI+CummMutRateNRTI+CummMutRate184)

Pred<-data.frame(NNRTI_based = c(0,ProbNNRTIFirst,ProbNRTIFirst,Prob184First), PI_based = rep(0,4),
                 row.names=c("PI","NNRTI","NRTI","184"))

#For PI treatments
ProbPIFirst<-CummMutRatePI/(CummMutRatePI+CummMutRateNRTI+CummMutRate184)
ProbNRTIFirst<-CummMutRateNRTI/(CummMutRatePI+CummMutRateNRTI+CummMutRate184)
Prob184First<-CummMutRate184/(CummMutRatePI+CummMutRateNRTI+CummMutRate184)
Pred[c(3,4,1),2]<-c(ProbNRTIFirst,Prob184First,ProbPIFirst)

pdf("whichMutationFirst.pdf",width=10,height=8)
barplot(as.matrix(Pred),col=c("yellow","red","lightblue","blue"),
        horiz=TRUE, cex.axis=2, main = "Which resistance comes first? Predictions",cex.main=2,cex.lab=1.6, cex.names=1.6)
#legend(x=0.2,y=0.2,legend=row.names(Pred))

text(0.14,0.65,"NNRTI first",cex=1.5)
text(0.58,0.65,"NRTI first",cex=1.5)
text(0.92,0.65,"184 first",,cex=1.5)

text(0.2,1.9,"PI first",cex=1.5)
text(0.65,1.9,"NRTI first",cex=1.5)
text(0.935,1.9,"184 first",,cex=1.5)
dev.off()

