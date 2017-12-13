#This R script reads the Lehman SRA fasta files and creates a csv file that has frequencies for each patient and each site for the transition mutations. 
# We filter the data before writing the file. 

listLehmanfiles<-list.files("Data/LehmanData/PleuniAlignments")

SraRunInfo<-read.csv("Data//LehmanData/SraRunInfo.csv")
SraRunInfo<-SraRunInfo[,which(names(SraRunInfo)%in%c("Run","SampleName"))]
#    data.frame(Run=SraRunInfo$Run, SampleName=SraRunInfo$SampleName)

#get unique part of filename
#substr(listLehmanfiles[1],1,gregexpr(".fasta",listLehmanfiles[1])[[1]]-1)

#which fasta files correspond to +1M data?
listtokeep<-c(); shortfilenames<-c()
for (i in 1:length(listLehmanfiles)){
filename=listLehmanfiles[i]
shortfilename<-substr(listLehmanfiles[i],1,gregexpr(".fasta",listLehmanfiles[1])[[1]]-1)
shortfilenames<-c(shortfilenames,shortfilename)
#find line in SraRunInfo
line<-which(SraRunInfo$Run==shortfilename)
samplename<-SraRunInfo$SampleName[line]
if (regexpr("\\+",samplename)[[1]]!=-1) {
    #print(samplename)
    listtokeep<-c(listtokeep,i) 
}
}

listLehmanfiles<-listLehmanfiles[listtokeep]
shortfilenames<-shortfilenames[listtokeep]

#get all the names of the fasta files for each patient in Lehman dataset
#determine all the frequencies and store in freqPatSite data.frame
#only needed if the stored data are not O

#Load libraries and necessary files from the baseRscript.Rmd
source('Rscripts/baseRscript.R')

###############
#make dataframe with frequencies for all non-muts for all patients for all sites filtered with the WT threshold.
freqPatTs_Lehman<-data.frame(row.names=shortfilenames)
CountDataLehman<-data.frame(pat=character(), pos=integer(),WTnt=character(),MutCount=integer(),WTcount=integer(),stringsAsFactors = F)

################
#OK, now we go through the fasta files again to get the frequencies

for (i in 1:length(listLehmanfiles)){ #for each fastafile
    filename=paste("../Data/LehmanData/PleuniAlignments/",listLehmanfiles[i],sep="")
    print(filename)
    patfasta<-read.dna(filename, format = "fasta",as.character=TRUE) #read the file
    
    #   MUST BE SHORTER maybe? 
    for (j in 1:687){#for each site in the sequence
        #print(paste('j',j))
        WT =  consensusB[j+297] #what is WT at site j? start at RT
        transitionnuc = transition(WT) #which nuc is the transition mutation?
        #check wether the neigboring sequences are the same
        # in order to check whether j is the first, second or third position of the codon, you can do
        if(j %in% seq(1,982,by=3)) {# first position
            #use only those sequences that are WT at pos 2 and 3 of the triplet
            goodsequences<-which(paste(patfasta[,j+1],patfasta[,j+2]) == paste(consensusB[c(j+297+1)],consensusB[(j+297+2)]) & patfasta[,j]%in% c(WT,transitionnuc))
            freqPatTs_Lehman[i,j]=length(which(patfasta[goodsequences,j]== transitionnuc))/length(goodsequences)
            CountDataLehman[length(CountDataLehman[,1])+1,]<-c(shortfilenames[i],j,WT,length(which(patfasta[goodsequences,j]== transitionnuc)),length(goodsequences))
        }
        if((j %in% seq(2,983,by=3))) {# second position
            #use only those sequences that are WT at pos 1 and 3 of the triplet
            goodsequences<-which(paste(patfasta[,j-1],patfasta[,j+1]) == paste(consensusB[c(j+297-1)],consensusB[(j+297+1)])& patfasta[,j]%in% c(WT,transitionnuc))
            freqPatTs_Lehman[i,j]=length(which(patfasta[goodsequences,j]==transitionnuc))/length(goodsequences)
            CountDataLehman[length(CountDataLehman[,1])+1,]<-c(shortfilenames[i],j,WT,length(which(patfasta[goodsequences,j]== transitionnuc)),length(goodsequences))
        }
        if((j %in% seq(3,984,by=3))) {# third position
            #use only those sequences that are WT at pos 1 and 2 of the triplet
            goodsequences<-which(paste(patfasta[,j-2],patfasta[,j-1]) == paste(consensusB[c(j+297-2)],consensusB[(j+297-1)])& patfasta[,j]%in% c(WT,transitionnuc))
            freqPatTs_Lehman[i,j]=length(which(patfasta[goodsequences,j]==transitionnuc))/length(goodsequences)
            CountDataLehman[length(CountDataLehman[,1])+1,]<-c(shortfilenames[i],j,WT,length(which(patfasta[goodsequences,j]== transitionnuc)),length(goodsequences))
        }
        if (length(goodsequences)==0)  freqPatTs_Lehman[i,j]<-NA
    }
}
write.csv(freqPatTs_Lehman,file="../Output/freqPatTs_LehmanRT.csv")    
write.csv(CountDataLehman,"../Output/CountDataLehmanRT.csv")

    