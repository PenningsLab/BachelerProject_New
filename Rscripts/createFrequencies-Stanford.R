#----
#title: "Calculating mutation frequencies using the Stanford fasta files"
#version: "May 2016"
#output: html_document
----

# Output as end result 
###Output of the file are 5 files with frequencies
###### one for PR naive: subtypeB-frequencies_stanford_pr_naive.csv
###### one for PR treated: subtypeB-frequencies_stanford_pr_treated.csv
###### one for RT naive: subtypeB-frequencies_stanford_rt_naive.csv
###### one for RT treated: subtypeB-frequencies_stanford_rt_treated.csv
###### a combined file with PR naive and RT naive together : subtypeB-frequencies_stanford_pr_rt_naive.csv

# Creating files with mutation frequencies of PR and RT 
# we want one file with frequencies for 984 sites (PR and part of RT)


#Load libraries and necessary files from the baseRscript.R
source('baseRscript.R')

## Read in the fasta files 
# PROTEASE
pr_B_naive<-read.dna('../Data/StanfordData/subtypeB-pr_naive_aligned.fasta', format = "fasta",as.character=TRUE)
nrow(pr_B_naive)

# REVERSE TRANSCRIPTASE
rt_B_naive<-read.dna('../Data/StanfordData/subtypeB-rt_naive_aligned.fasta', format = "fasta",as.character=TRUE)
nrow(rt_B_naive)

## Calculate frequencies 
# We have a separate fasta file for PR and one for RT, so we need separate loops

### PROTEASE
#CALCULATE FREQUENCIES OF TRANSITIONS
## for naive sequences 
patfasta<-pr_B_naive  #naive
i<-1
freqPatTs_threshold_pr_naive<-data.frame()

for (j in 1:297){#for each site in the sequence
    WT=	consensusB[j] #what is WT at site j?
    #print(WT)
    if(j %in% seq(1,295,by=3)) {# first position
        goodsequences<-which(paste(patfasta[,j+1],patfasta[,j+2]) == paste(consensusB[c(j+1)],consensusB[(j+2)]))
        if (WT=="c"){freqPatTs_threshold_pr_naive[i,j]=length(which(patfasta[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_pr_naive[i,j]=length(which(patfasta[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_pr_naive[i,j]=length(which(patfasta[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_pr_naive[i,j]=length(which(patfasta[goodsequences,j]=="a"))/length(goodsequences)
        }}
    
    if((j %in% seq(2,296,by=3))) {# second position
        goodsequences<-which(paste(patfasta[,j-1],patfasta[,j+1]) == paste(consensusB[c(j-1)],consensusB[(j+1)]))
        if (WT=="c"){freqPatTs_threshold_pr_naive[i,j]=length(which(patfasta[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_pr_naive[i,j]=length(which(patfasta[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_pr_naive[i,j]=length(which(patfasta[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_pr_naive[i,j]=length(which(patfasta[goodsequences,j]=="a"))/length(goodsequences)
        }}
    
    if((j %in% seq(3,297,by=3))) {# third position
        goodsequences<-which(paste(patfasta[,j-2],patfasta[,j-1]) == paste(consensusB[c(j-2)],consensusB[(j-1)]))
        if (WT=="c"){freqPatTs_threshold_pr_naive[i,j]=length(which(patfasta[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_pr_naive[i,j]=length(which(patfasta[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_pr_naive[i,j]=length(which(patfasta[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_pr_naive[i,j]=length(which(patfasta[goodsequences,j]=="a"))/length(goodsequences)
        }}
}

## write the frequencies  of PR into a file
#write.table(as.data.frame(freqPatTs_threshold_pr_naive),file='subtypeB-frequencies_stanford_pr_naive.csv',row.names=F,quote=F,sep=',')


### REVERSE TRANSCRIPTASE
#Limit the length of RT to that one of the bacheler dataset, which is 984 nucleotide position, of whih 687 are from RT>


#length of the RT region
length_rt<-687

#actually, an alignment against both PR and RT is the best for the numbering. For now only against RT so make a consensusB of only RT 
consensusB_rt<-consensusB[298:984]


#CALCULATE FREQUENCIES OF TRANSITIONS  
## for naive sequences 
patfasta_rt<-rt_B_naive  #naive
i<-1
freqPatTs_threshold_rt_naive<-data.frame()

for (j in 1:length_rt){#for each site in the sequence
    #print(j)
    WT=	consensusB_rt[j] #what is WT at site j?
    
    if(j %in% seq(1,length_rt-2,by=3)) {# first position
        goodsequences<-which(paste(patfasta_rt[,j+1],patfasta_rt[,j+2]) == paste(consensusB_rt[c(j+1)],consensusB_rt[(j+2)]))
        
        if (WT=="c"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="a"))/length(goodsequences)
        }}
    
    if((j %in% seq(2,length_rt-1,by=3))) {# second position
        goodsequences<-which(paste(patfasta_rt[,j-1],patfasta_rt[,j+1]) == paste(consensusB_rt[c(j-1)],consensusB_rt[(j+1)]))
        if (WT=="c"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="a"))/length(goodsequences)
        }}
    
    if((j %in% seq(3,length_rt,by=3))) {# third position
        goodsequences<-which(paste(patfasta_rt[,j-2],patfasta_rt[,j-1]) == paste(consensusB_rt[c(j-2)],consensusB_rt[(j-1)]))
        if (WT=="c"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="a"))/length(goodsequences)
        }}
}

# export the result for RT
# write.table(as.data.frame(freqPatTs_threshold_rt_naive),file='subtypeB-frequencies_stanford_rt_naive.csv',row.names=F,quote=F,sep=',')


## Combine frequencies of PR and of RT into a single file 
# add PR and RT together.
colnames(freqPatTs_threshold_pr_naive)<-1:297
colnames(freqPatTs_threshold_rt_naive)<-298:984

freqPatTs_threshold_pr_rt_naive<-as.data.frame(c(freqPatTs_threshold_pr_naive,freqPatTs_threshold_rt_naive))
freqPatTs_threshold_pr_rt_naive2<-data.frame(num=1:984, colMeansTs0_stanford=t(freqPatTs_threshold_pr_rt_naive))

write.table(as.data.frame(freqPatTs_threshold_pr_rt_naive2),file='../Output/freqPatTs_Stanford.csv',row.names=F,quote=F,sep=',')

