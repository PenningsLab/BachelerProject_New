#Read the stored frequencies rather than calculating frequencies again
#read the stored data
#We read a file with WT 0 as threshold, meaning no threshold freqPatTs_Bacheler_Threshold05.csv)

read.table("./Output/freqPatTs_Bacheler_Threshold1.csv",sep=",",header=TRUE,row.names=1)->freqPatTsFilter
#calculate mean frequencies
MeanFreq<-apply(freqPatTsFilter, 2 , mean, na.rm=TRUE)

#Calculate mean freq for each site and bootstrap / no longer done

#Get conf intervals by bootsstrapping
#PSP Nov 11 2015: I don't think we need te bootstrapping, so I put it in an if (FALSE) statement.
#PSP June 2017: it looks like we may need it after all. 
if (FALSE){
    btmeansTs<-data.frame(row.names=names(freqPatTsFilter)[1:984]) #each row is a site, each column is a bootstrapped mean
    numbootstraps=1000
    for (j in 1:984){# start with the first site
        print(j)
        for (k in 1:numbootstraps){ # first iteration
            btmeansTs[j,k]=mean(sample(freqPatTsFilter[,j],length(freqPatTsFilter[,j]),replace = TRUE),na.rm=TRUE)
        }
    }
    btmeansSorted<-t(apply(btmeansTs,1,sort))
    lowerConf<-btmeansSorted[,floor(0.025*numbootstraps)]
    upperConf<-btmeansSorted[,floor((1-0.025)*numbootstraps)]
}

## Create overview dataframe and plot site frequency spectra 
#Only synonymous, non-synomous and stop codons are considered
#- for each mutation, determine whether it is synonymous, non-synonymous or creates a stop
#- add information on resistance  positions

#PSP Nov 11 2015 I removed lowerConf and upperConf here because we no longer calculate them
#PSP June 2017 added them again (lowerConf and upperConf)
#PSP Nov 11 2015 I renamed x OverviewDFilter and newdata OverviewDFilterOrderedByFreq
if (exists("lowerConf")){
    OverviewDFilter<-data.frame(num=1:984,MeanFreq,TypeOfSite,lowerConf,upperConf)
}else OverviewDFilter<-data.frame(num=1:984,MeanFreq,TypeOfSite)

OverviewDFilter$WTnt<-consensusB[1:984]

#Mut rates and sel coefficients
read.csv("./Data/HIVMutRates/HIVMutRates.csv")->mutrates
OverviewDFilter$TSmutrate<-0
OverviewDFilter$TSmutrate[OverviewDFilter$WTnt=="a"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="AG"]
OverviewDFilter$TSmutrate[OverviewDFilter$WTnt=="c"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="CU"]
OverviewDFilter$TSmutrate[OverviewDFilter$WTnt=="g"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="GA"]
OverviewDFilter$TSmutrate[OverviewDFilter$WTnt=="t"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="UC"]

#Add mutrates from Zanini paper
OverviewDFilter$TSmutZan<-0
OverviewDFilter$TSmutZan[OverviewDFilter$WTnt=="a"]<-mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="AG"]
OverviewDFilter$TSmutZan[OverviewDFilter$WTnt=="c"]<-mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="CU"]
OverviewDFilter$TSmutZan[OverviewDFilter$WTnt=="g"]<-mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="GA"]
OverviewDFilter$TSmutZan[OverviewDFilter$WTnt=="t"]<-mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="UC"]


for (i in 1:984){
    OverviewDFilter$EstSelCoeff[i] = EstimatedS(OverviewDFilter$TSmutrate[i],OverviewDFilter$MeanFreq[i])
    OverviewDFilter$EstSelCoeffZan[i] = EstimatedS(OverviewDFilter$TSmutZan[i],OverviewDFilter$MeanFreq[i])
}

#WT AAs 
OverviewDFilter$WTAA<-""
for (i in 1:984){
    if (i%%3==1) OverviewDFilter$WTAA[i] = seqinr::translate(OverviewDFilter$WTnt[c(i,i+1,i+2)])
    if (i%%3==2) OverviewDFilter$WTAA[i] = seqinr::translate(OverviewDFilter$WTnt[c(i-1,i,i+1)])
    if (i%%3==0) OverviewDFilter$WTAA[i] = seqinr::translate(OverviewDFilter$WTnt[c(i-2,i-1,i)])
}

OverviewDFilter$MUTAA<-""
#MUT AAs
for (i in 1:984){
    if (i%%3==1) OverviewDFilter$MUTAA[i] = seqinr::translate(c(transition(OverviewDFilter$WTnt[i]),OverviewDFilter$WTnt[c(i+1,i+2)]))
    if (i%%3==2) OverviewDFilter$MUTAA[i] = seqinr::translate(c(OverviewDFilter$WTnt[c(i-1)],transition(OverviewDFilter$WTnt[i]),OverviewDFilter$WTnt[c(i+1)]))
    if (i%%3==0) OverviewDFilter$MUTAA[i] = seqinr::translate(c(OverviewDFilter$WTnt[c(i-2,i-1)],transition(OverviewDFilter$WTnt[i])))
}

#Add whether AA change is drastic 
OverviewDFilter$bigAAChange<-0

for(i in 1:nrow(OverviewDFilter)){
    WT <- amCat(OverviewDFilter[i,'WTAA'])
    MUT <- amCat(OverviewDFilter[i,'MUTAA'])
    if (WT == MUT){ OverviewDFilter$bigAAChange[i] <- 0 
    }else{
        OverviewDFilter$bigAAChange[i] <- 1
    }
}

#Add whetehr makes CpG 
OverviewDFilter$makesCpG <- 0
for(i in 1:nrow(OverviewDFilter)){
    trip <- OverviewDFilter$WTnt[c(i-1, i, i + 1)]
    if(trip[1] == "c" & trip[2] == "a" ){
        OverviewDFilter$makesCpG[i] <- 1
    }
    if(trip[2] == "t" & trip[3] == "g"){
        OverviewDFilter$makesCpG[i] <- 1
    }
}

OverviewDFilter$ProRT<-c(rep("Pro",297),rep("RT",687))
OverviewDFilter$AAnum<-c(sort(rep(1:99,3)),sort(rep(1:229,3)))

#colors
OverviewDFilter$color<-""
for (i in 1:984){
    if (OverviewDFilter$TypeOfSite[i]=="syn") OverviewDFilter$color[i] = "darkolivegreen3"
    if (OverviewDFilter$TypeOfSite[i]=="nonsyn") OverviewDFilter$color[i] = "red"
    if (OverviewDFilter$TypeOfSite[i]=="stop") OverviewDFilter$color[i] = "black"
    if (OverviewDFilter$TypeOfSite[i]=="res") OverviewDFilter$color[i] = "purple"
    if (OverviewDFilter$TypeOfSite[i]=="overlap") OverviewDFilter$color[i] = "grey"
}

write.csv(OverviewDFilter,"./Output/OverviewSelCoeff_BachelerFilter.csv")

