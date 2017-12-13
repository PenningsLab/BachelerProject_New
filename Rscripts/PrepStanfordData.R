# Analysis using Stanford website data

## Part 1 :  Download and preprocessing

### Protease download and fasta file generation 
#The dataset that you can download from Stanford (http://hivdb.stanford.edu/pages/geno-rx-datasets.html) is called. PR.txt and RT.txt (see below)

#We will process this file to remove non-B subtypes, empty sequences, treated sequences, multiple sequences per patient and sequences containing bad characters. 


#### retrieve sequences from PR.xt and clean up 
data<-read.csv('PR.txt',sep=',')
# only subtype B
dataB<-subset(data,data$Subtype=='B')
# remove empty sequences 
dataBn<-subset(dataB,dataB$NASeq!="")
# only naive sequences 
dataBnaive<-subset(dataBn,dataBn$PIList=='None')
# only 1 sequence per patient 
dataBnaive_unique<-dataBnaive[!duplicated(dataBnaive$PtID),]
# remove sequences with bad characters
dataBnaive_unique2<-dataBnaive_unique[-grep('~',dataBnaive_unique$NASeq),]

# combine  identifier en sequence
##export naive sequences but still not aligned
dataBnaive_unique3<-cbind(as.character(dataBnaive_unique2$PtID),as.character(dataBnaive_unique2$NASeq))
write.table(as.data.frame(dataBnaive_unique3),file='subtypeB-pr_naive.csv',row.names=F,quote=F,sep=',')


#### generate a aligned fasta file from the sequences 
#Several steps to convert the csv with sequences into a fasta file

#** to fasta format **
#    vi subtypeB-pr_naive.csv # add >
#vi subtypeB-pr_treated.csv
#tr ',' '\n' < subtypeB-pr_naive.csv > subtypeB-pr_naive.fasta
#tr ',' '\n' < subtypeB-pr_treated.csv > subtypeB-pr_treated.fasta

#* transfer to server **
#    scp *fasta cai:/home/ktheys0
#grep '\.' subtypeB-rt_naive.fasta  |wc

#** Align sequences using C++ library  **
#    This tool for pairwise alignment against HXB2 is available at http://regatools.med.kuleuven.be/sequencetool/sequencetool.wt
#
#SequenceTool PR subtypeB-pr_naive.fasta --exportKind GlobalAlignment --exportAlphabet Nucleotides --exportWithInsertions no  > subtypeB-pr_naive_aligned.fasta
#SequenceTool PR subtypeB-pr_treated.fasta --exportKind GlobalAlignment --exportAlphabet Nucleotides --exportWithInsertions no  > subtypeB-pr_treated_aligned.fasta

#** remove HXB2 from alignement ** 
#    vi subtypeB-pr_naive_aligned.fasta
#vi subtypeB-pr_treated_aligned.fasta

#** Convert nucleotide sequences to amino acids:   **
#    SequenceTool PR subtypeB-pr_naive_aligned.fasta --exportKind PositionTable  > subtypeB-pr_naive_mut.csv
#SequenceTool PR subtypeB-pr_treated_aligned.fasta --exportKind PositionTable   > subtypeB-pr_treated_mut.csv

#** get from server  ** 
#    scp cai:/home/ktheys0/subtypeB-pr* .
#scp cai:/home/ktheys0/subtypeB-pr*csv .


### Reverse transcriptase  download and fasta file generation 
#Similar steps as for the PR.txt file
#### retrieve sequences from RT.txt
#same steps as for PR
data<-read.csv('RT.txt',sep='\t')
dataB<-subset(data,data$Subtype=='B')
dataBn<-subset(dataB,dataB$NASeq!="")
dataBnaive<-subset(dataBn,dataBn$RTIList=='None')

dataBnaive_unique<-dataBnaive[!duplicated(dataBnaive$PtID),]

# remove sequences with bad character
dataBnaive_unique2<-dataBnaive_unique[-grep('~',dataBnaive_unique$NASeq),]

# combine  identifier en sequence
##export naive sequences but still not aligned
dataBnaive_unique3<-cbind(as.character(dataBnaive_unique2$PtID),as.character(dataBnaive_unique2$NASeq))
write.table(as.data.frame(dataBnaive_unique3),file='subtypeB-rt_naive.csv',row.names=F,quote=F,sep=',')

#### generate a aligned fasta file from the sequences 
#vi to add '>' 
#tr ',' '\n' < subtypeB-rt_naive.csv > subtypeB-rt_naive.fasta
#tr ',' '\n' < subtypeB-rt_treated.csv > subtypeB-rt_treated.fasta

#** transfer to server ** 
#    scp *fasta cai:/home/ktheys0
#grep '\.' subtypeB-rt_naive.fasta  |wc

#** pairwise alignment against HXB2 ** 
#    SequenceTool RT subtypeB-rt_naive.fasta --exportKind GlobalAlignment --exportAlphabet Nucleotides --exportWithInsertions no  > subtypeB-rt_naive_aligned.fasta
#SequenceTool RT subtypeB-rt_treated.fasta --exportKind GlobalAlignment --exportAlphabet Nucleotides --exportWithInsertions no  > subtypeB-rt_treated_aligned.fasta

#** remove hxb2 sequence**
#    vi subtypeB-rt_naive_aligned.fasta
#vi subtypeB-rt_treated_aligned.fasta

#** Convert nucleotide sequences to amino acids  **
#    SequenceTool PR subtypeB-rt_naive_aligned.fasta --exportKind PositionTable  > subtypeB-rt_naive_mut.csv
#SequenceTool PR subtypeB-rt_treated_aligned.fasta --exportKind PositionTable   > subtypeB-rt_treated_mut.csv

#** download from server ** 
#    scp cai:/home/ktheys0/subtypeB-rt*fasta .
#scp cai:/home/ktheys0/subtypeB-rt*csv .

## Step 2:  R analysis for calculating frequencies  

### At the amino acid level  

#### Protease frequencies AA 
# NAIVE SEQUENCES
datamutnaive<-read.csv('subtypeB-pr_naive_mut.csv')
mutnaive<-c()
for(i in 2:ncol(datamutnaive))
{
    ref<-datamutnaive[1,i]
    compare<-datamutnaive[2:nrow(datamutnaive),i]
    mutnaive<-c(mutnaive,sum(compare[compare!=""]==ref)/length(compare[compare!=""])*100)
}
plot(mutnaive,xaxt='n')
axis(1,at=seq(1,100,by=10),seq(1,100,by=10))

# TREATED SEQUENCES 
datamuttreated<-read.csv('subtypeB-pr_treated_mut.csv')
muttreat<-c()
for(i in 2:ncol(datamuttreated))
{
    ref<-datamuttreated[1,i]
    compare<-datamuttreated[2:nrow(datamuttreated),i]
    muttreat<-c(muttreat,sum(compare[compare!=""]==ref)/length(compare[compare!=""])*100)
}
plot(muttreat,xaxt='n')
axis(1,at=seq(1,100,by=10),seq(1,100,by=10))
points(mutnaive,col='green')

plot(muttreat,mutnaive)
```

#### Reverse transcriptase AA frequencies 
```{r,eval=F}
# NAIVE
datamutnaive<-read.csv('subtypeB-rt_naive_mut.csv')
mutnaive<-c()
for(i in 2:ncol(datamutnaive))
{
    ref<-datamutnaive[1,i]
    compare<-datamutnaive[2:nrow(datamutnaive),i]
    mutnaive<-c(mutnaive,sum(compare[compare!=""]==ref)/length(compare[compare!=""])*100)
}
plot(mutnaive,xaxt='n')
axis(1,at=seq(1,100,by=10),seq(1,100,by=10))

# TREATED
datamuttreated<-read.csv('subtypeB-pr_treated_mut.csv')
muttreat<-c()
for(i in 2:ncol(datamuttreated))
{
    ref<-datamuttreated[1,i]
    compare<-datamuttreated[2:nrow(datamuttreated),i]
    muttreat<-c(muttreat,sum(compare[compare!=""]==ref)/length(compare[compare!=""])*100)
}
plot(muttreat,xaxt='n')
axis(1,at=seq(1,100,by=10),seq(1,100,by=10))
points(mutnaive,col='green')

plot(muttreat,mutnaive)
```

### At the nucleotide level 

#### Protease frequencies nt frequencies  
```{r}
#READ IN FILES  
library(ape);library(seqinr);library(scales)
consensusB<-read.dna('consprrt.fasta',format = "fasta",as.character=TRUE)
pr_B_naive<-read.dna('subtypeB-pr_naive_aligned.fasta', format = "fasta",as.character=TRUE)
pr_B_treated<-read.dna('subtypeB-pr_treated_aligned.fasta', format = "fasta",as.character=TRUE)
nrow(pr_B_naive);nrow(pr_B_treated)
```

```{r}
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

## for treated sequences 
patfasta<-pr_B_treated
i<-1
freqPatTs_threshold_pr_treated<-data.frame()

for (j in 1:297){#for each site in the sequence
    #print(j)
    WT=	consensusB[j] #what is WT at site j?
    
    if(j %in% seq(1,295,by=3)) {# first position
        goodsequences<-which(paste(patfasta[,j+1],patfasta[,j+2]) == paste(consensusB[c(j+1)],consensusB[(j+2)]))
        if (WT=="c"){freqPatTs_threshold_pr_treated[i,j]=length(which(patfasta[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_pr_treated[i,j]=length(which(patfasta[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_pr_treated[i,j]=length(which(patfasta[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_pr_treated[i,j]=length(which(patfasta[goodsequences,j]=="a"))/length(goodsequences)
        }}
    
    if((j %in% seq(2,296,by=3))) {# second position
        goodsequences<-which(paste(patfasta[,j-1],patfasta[,j+1]) == paste(consensusB[c(j-1)],consensusB[(j+1)]))
        if (WT=="c"){freqPatTs_threshold_pr_treated[i,j]=length(which(patfasta[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_pr_treated[i,j]=length(which(patfasta[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_pr_treated[i,j]=length(which(patfasta[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_pr_treated[i,j]=length(which(patfasta[goodsequences,j]=="a"))/length(goodsequences)
        }}
    
    if((j %in% seq(3,297,by=3))) {# third position
        goodsequences<-which(paste(patfasta[,j-2],patfasta[,j-1]) == paste(consensusB[c(j-2)],consensusB[(j-1)]))
        if (WT=="c"){freqPatTs_threshold_pr_treated[i,j]=length(which(patfasta[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_pr_treated[i,j]=length(which(patfasta[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_pr_treated[i,j]=length(which(patfasta[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_pr_treated[i,j]=length(which(patfasta[goodsequences,j]=="a"))/length(goodsequences)
        }}
}

## write 
write.table(as.data.frame(freqPatTs_threshold_pr_naive),file='subtypeB-frequencies_stanford_pr_naive.csv',row.names=F,quote=F,sep=',')
write.table(as.data.frame(freqPatTs_threshold_pr_treated),file='subtypeB-frequencies_stanford_pr_treated.csv',row.names=F,quote=F,sep=',')
```

#### Reverse transcriptase  frequencies nt frequencies  
```{r}
#READ IN FILES  
rt_B_naive<-read.dna('subtypeB-rt_naive_aligned.fasta', format = "fasta",as.character=TRUE)
rt_B_treated<-read.dna('subtypeB-rt_treated_aligned.fasta', format = "fasta",as.character=TRUE)
nrow(rt_B_naive);nrow(rt_B_treated)

#actually, an alignment against both PR and RT is the best for the numbering. For now only against RT so make a consensusB of only RT 
consensusB_rt<-consensusB[298:1977] 
```

```{r}
#CALCULATE FREQUENCIES OF TRANSITIONS  
## for naive sequences 
patfasta_rt<-rt_B_naive  #naive
i<-1
freqPatTs_threshold_rt_naive<-data.frame()

for (j in 1:1680){#for each site in the sequence
    #print(j)
    WT=	consensusB_rt[j] #what is WT at site j?
    
    if(j %in% seq(1,1678,by=3)) {# first position
        goodsequences<-which(paste(patfasta_rt[,j+1],patfasta_rt[,j+2]) == paste(consensusB_rt[c(j+1)],consensusB_rt[(j+2)]))
        
        if (WT=="c"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="a"))/length(goodsequences)
        }}
    
    if((j %in% seq(2,1679,by=3))) {# second position
        goodsequences<-which(paste(patfasta_rt[,j-1],patfasta_rt[,j+1]) == paste(consensusB_rt[c(j-1)],consensusB_rt[(j+1)]))
        if (WT=="c"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="a"))/length(goodsequences)
        }}
    
    if((j %in% seq(3,1680,by=3))) {# third position
        goodsequences<-which(paste(patfasta_rt[,j-2],patfasta_rt[,j-1]) == paste(consensusB_rt[c(j-2)],consensusB_rt[(j-1)]))
        if (WT=="c"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_rt_naive[i,j]=length(which(patfasta_rt[goodsequences,j]=="a"))/length(goodsequences)
        }}
}

## for treated sequences 
patfasta<-rt_B_treated  #treated
i<-1
freqPatTs_threshold_rt_treated<-data.frame()
for (j in 1:1680){#for each site in the sequence
    #print(j)
    WT=	consensusB_rt[j] #what is WT at site j?
    
    if(j %in% seq(1,1678,by=3)) {# first position
        goodsequences<-which(paste(patfasta_rt[,j+1],patfasta_rt[,j+2]) == paste(consensusB_rt[c(j+1)],consensusB_rt[(j+2)]))
        
        if (WT=="c"){freqPatTs_threshold_rt_treated[i,j]=length(which(patfasta_rt[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_rt_treated[i,j]=length(which(patfasta_rt[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_rt_treated[i,j]=length(which(patfasta_rt[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_rt_treated[i,j]=length(which(patfasta_rt[goodsequences,j]=="a"))/length(goodsequences)
        }}
    
    if((j %in% seq(2,1679,by=3))) {# second position
        goodsequences<-which(paste(patfasta_rt[,j-1],patfasta_rt[,j+1]) == paste(consensusB_rt[c(j-1)],consensusB_rt[(j+1)]))
        if (WT=="c"){freqPatTs_threshold_rt_treated[i,j]=length(which(patfasta_rt[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_rt_treated[i,j]=length(which(patfasta_rt[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_rt_treated[i,j]=length(which(patfasta_rt[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_rt_treated[i,j]=length(which(patfasta_rt[goodsequences,j]=="a"))/length(goodsequences)
        }}
    
    if((j %in% seq(3,1680,by=3))) {# third position
        goodsequences<-which(paste(patfasta_rt[,j-2],patfasta_rt[,j-1]) == paste(consensusB_rt[c(j-2)],consensusB_rt[(j-1)]))
        if (WT=="c"){freqPatTs_threshold_rt_treated[i,j]=length(which(patfasta_rt[goodsequences,j]=="t"))/length(goodsequences)}
        if (WT=="t"){freqPatTs_threshold_rt_treated[i,j]=length(which(patfasta_rt[goodsequences,j]=="c"))/length(goodsequences)}
        if (WT=="a"){freqPatTs_threshold_rt_treated[i,j]=length(which(patfasta_rt[goodsequences,j]=="g"))/length(goodsequences)}
        if (WT=="g"){freqPatTs_threshold_rt_treated[i,j]=length(which(patfasta_rt[goodsequences,j]=="a"))/length(goodsequences)
        }}
}``

# export 
write.table(as.data.frame(freqPatTs_threshold_rt_naive),file='subtypeB-frequencies_stanford_rt_naive.csv',row.names=F,quote=F,sep=',')
write.table(as.data.frame(freqPatTs_threshold_rt_treated),file='subtypeB-frequencies_stanford_rt_treated.csv',row.names=F,quote=F,sep=',')
```