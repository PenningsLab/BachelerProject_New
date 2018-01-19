#Script to analyse the frequency data and associate with features for Zanini data. 
#Read the csv files 

source('Rscripts/baseRscript.R')

#Read the stored frequencies rather than calculating frequencies again
read.table("Output/freqPatTs_Zanini.csv",sep=",",header=TRUE,row.names=1)->freqPatTsZanini
colMeansTsZanini<-apply(freqPatTsZanini, 2 , mean, na.rm=TRUE)

## Create overview dataframe and plot site frequency spectra 
#Only synonymous, non-synomous and stop codons are considered
#- for each mutation, determine whether it is synonymous, non-synonymous or creates a stop
#- add information on resistance  positions

#PSP Nov 11 2015 I renamed x OverviewDFZanini and newdata OverviewDFZaniniOrderedByFreq
numsitesZanini<-length(colMeansTsZanini)
OverviewDFZanini<-data.frame(num=1:numsitesZanini,colMeansTsZanini)
OverviewDFZanini$TypeOfSite<-TypeOfSite[1:numsitesZanini]
OverviewDFZanini$TypeOfSite[1:39]<-"overlap"
OverviewDFZanini$TypeOfSite[which(consensusB != consensusC | consensusC != consensus01AE)]<-"exclude"
OverviewDFZanini$WTnt<-consensusB[1:numsitesZanini]

#Mut rates and sel coefficients from Abrams paper
read.csv("Data/HIVMutRates/HIVMutRates.csv")->mutrates
OverviewDFZanini$TSmutrate<-0
OverviewDFZanini$TSmutrate[OverviewDFZanini$WTnt=="a"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="AG"]
OverviewDFZanini$TSmutrate[OverviewDFZanini$WTnt=="c"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="CU"]
OverviewDFZanini$TSmutrate[OverviewDFZanini$WTnt=="g"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="GA"]
OverviewDFZanini$TSmutrate[OverviewDFZanini$WTnt=="t"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="UC"]

OverviewDFZanini$TSmutZan<-0
OverviewDFZanini$TSmutZan[OverviewDFZanini$WTnt=="a"]<-mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="AG"]
OverviewDFZanini$TSmutZan[OverviewDFZanini$WTnt=="c"]<-mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="CU"]
OverviewDFZanini$TSmutZan[OverviewDFZanini$WTnt=="g"]<-mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="GA"]
OverviewDFZanini$TSmutZan[OverviewDFZanini$WTnt=="t"]<-mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="UC"]


for (i in 1:984){
    OverviewDFZanini$EstSelCoeff[i] = EstimatedS(OverviewDFZanini$TSmutrate[i],OverviewDFZanini$colMeansTsZanini[i])
    OverviewDFZanini$EstSelCoeffZan[i] = EstimatedS(OverviewDFZanini$TSmutZan[i],OverviewDFZanini$colMeansTsZanini[i])}

#OverviewDFZanini$EstSelCoeff= OverviewDFZanini$TSmutrate/OverviewDFZanini$colMeansTsZanini
#OverviewDFZanini$EstSelCoeff[OverviewDFZanini$EstSelCoeff>1]<-1

#WT AAs 
OverviewDFZanini$WTAA<-""
for (i in 1:numsitesZanini){
    if (i%%3==1) OverviewDFZanini$WTAA[i] = seqinr::translate(OverviewDFZanini$WTnt[c(i,i+1,i+2)])
    if (i%%3==2) OverviewDFZanini$WTAA[i] = seqinr::translate(OverviewDFZanini$WTnt[c(i-1,i,i+1)])
    if (i%%3==0) OverviewDFZanini$WTAA[i] = seqinr::translate(OverviewDFZanini$WTnt[c(i-2,i-1,i)])
}

OverviewDFZanini$MUTAA<-""
#MUT AAs
for (i in 1:numsitesZanini){
    if (i%%3==1) OverviewDFZanini$MUTAA[i] = seqinr::translate(c(transition(OverviewDFZanini$WTnt[i]),OverviewDFZanini$WTnt[c(i+1,i+2)]))
    if (i%%3==2) OverviewDFZanini$MUTAA[i] = seqinr::translate(c(OverviewDFZanini$WTnt[c(i-1)],transition(OverviewDFZanini$WTnt[i]),OverviewDFZanini$WTnt[c(i+1)]))
    if (i%%3==0) OverviewDFZanini$MUTAA[i] = seqinr::translate(c(OverviewDFZanini$WTnt[c(i-2,i-1)],transition(OverviewDFZanini$WTnt[i])))
}

write.csv(OverviewDFZanini,"Output/OverviewSelCoeffZanini.csv")