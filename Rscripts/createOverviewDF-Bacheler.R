#Read the stored frequencies rather than calculating frequencies again
#read the stored data
#We read a file with WT 0 as threshold, meaning no threshold)
read.table("../Output/freqPatTs_Bacheler.csv",sep=",",header=TRUE,row.names=1)->freqPatTs0
#calculate mean frequencies
MeanFreq<-apply(freqPatTs0, 2 , mean, na.rm=TRUE)

#Calculate mean freq for each site and bootstrap / no longer done

#Get conf intervals by bootsstrapping
#PSP Nov 11 2015: I don't think we need te bootstrapping, so I put it in an if (FALSE) statement.
#PSP June 2017: it looks like we may need it after all. 
if (TRUE){
    btmeansTs<-data.frame(row.names=names(freqPatTs0)[1:984]) #each row is a site, each column is a bootstrapped mean
    numbootstraps=1000
    for (j in 1:984){# start with the first site
        print(j)
        for (k in 1:numbootstraps){ # first iteration
            btmeansTs[j,k]=mean(sample(freqPatTs0[,j],length(freqPatTs0[,j]),replace = TRUE),na.rm=TRUE)
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
#PSP Nov 11 2015 I renamed x OverviewDF and newdata OverviewDFOrderedByFreq
OverviewDF<-data.frame(num=1:984,MeanFreq,TypeOfSite,lowerConf,upperConf)

OverviewDF$WTnt<-consensusB[1:984]

#Mut rates and sel coefficients
read.csv("../Data/HIVMutRates/HIVMutRates.csv")->mutrates
OverviewDF$TSmutrate<-0
OverviewDF$TSmutrate[OverviewDF$WTnt=="a"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="AG"]
OverviewDF$TSmutrate[OverviewDF$WTnt=="c"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="CU"]
OverviewDF$TSmutrate[OverviewDF$WTnt=="g"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="GA"]
OverviewDF$TSmutrate[OverviewDF$WTnt=="t"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="UC"]

for (i in 1:984){
    OverviewDF$EstSelCoeff[i] = EstimatedS(OverviewDF$TSmutrate[i],OverviewDF$MeanFreq[i])
}

#WT AAs 
OverviewDF$WTAA<-""
for (i in 1:984){
    if (i%%3==1) OverviewDF$WTAA[i] = seqinr::translate(OverviewDF$WTnt[c(i,i+1,i+2)])
    if (i%%3==2) OverviewDF$WTAA[i] = seqinr::translate(OverviewDF$WTnt[c(i-1,i,i+1)])
    if (i%%3==0) OverviewDF$WTAA[i] = seqinr::translate(OverviewDF$WTnt[c(i-2,i-1,i)])
}

OverviewDF$MUTAA<-""
#MUT AAs
for (i in 1:984){
    if (i%%3==1) OverviewDF$MUTAA[i] = seqinr::translate(c(transition(OverviewDF$WTnt[i]),OverviewDF$WTnt[c(i+1,i+2)]))
    if (i%%3==2) OverviewDF$MUTAA[i] = seqinr::translate(c(OverviewDF$WTnt[c(i-1)],transition(OverviewDF$WTnt[i]),OverviewDF$WTnt[c(i+1)]))
    if (i%%3==0) OverviewDF$MUTAA[i] = seqinr::translate(c(OverviewDF$WTnt[c(i-2,i-1)],transition(OverviewDF$WTnt[i])))
}

#Add whether AA change is drastic 
OverviewDF$bigAAChange<-0

for(i in 1:nrow(OverviewDF)){
    WT <- amCat(OverviewDF[i,'WTAA'])
    MUT <- amCat(OverviewDF[i,'MUTAA'])
    if (WT == MUT){ OverviewDF$bigAAChange[i] <- 0 
    }else{
        OverviewDF$bigAAChange[i] <- 1
    }
}

#Add whetehr makes CpG 
OverviewDF$makesCpG <- 0
for(i in 1:nrow(OverviewDF)){
    trip <- OverviewDF$WTnt[c(i-1, i, i + 1)]
    if(trip[1] == "c" & trip[2] == "a" ){
        OverviewDF$makesCpG[i] <- 1
    }
    if(trip[2] == "t" & trip[3] == "g"){
        OverviewDF$makesCpG[i] <- 1
    }
}

OverviewDF$ProRT<-c(rep("Pro",297),rep("RT",687))
OverviewDF$AAnum<-c(sort(rep(1:99,3)),sort(rep(1:229,3)))

#colors
OverviewDF$color<-""
for (i in 1:984){
    if (OverviewDF$TypeOfSite[i]=="syn") OverviewDF$color[i] = "darkolivegreen3"
    if (OverviewDF$TypeOfSite[i]=="nonsyn") OverviewDF$color[i] = "red"
    if (OverviewDF$TypeOfSite[i]=="stop") OverviewDF$color[i] = "black"
    if (OverviewDF$TypeOfSite[i]=="res") OverviewDF$color[i] = "purple"
    if (OverviewDF$TypeOfSite[i]=="overlap") OverviewDF$color[i] = "grey"
}

write.csv(OverviewDF,"../Output/OverviewSelCoeff_Bacheler.csv")
