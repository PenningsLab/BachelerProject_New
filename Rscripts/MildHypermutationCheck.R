#Check for mild G to A hypermutation 

#This R script reads the Bacheler et al fasta files and creates a csv file that has frequencies for each patient and each site for the transition mutations. 
# We filter the data before writing the file. 

#new addition
#1. only take patients with a consensus sequence at Day0
#2. only consider changed nucleotides when the two other nucleotides in a codon are not mutated

#get all the names of the fasta files for each patient in Bacheler2000 dataset
#determine all the frequencies and store in freqPatSite data.frame
#only needed if the stored data are not O

#Load libraries and necessary files from the baseRscript.Rmd
setwd("~/Documents/Git/bachelerProject")
source('Rscripts/baseRscript.R')

#Read the correct fastafiles.
listfastafiles<-list.files("Data/BachelerFiles//FASTAfiles")
lengthlatersequenceslist<-c();lengthallsequenceslist<-c()

###############   
#First, let's remove patients with too few sequences (less than 5). 
Numseqs<-c()
for (i in 1:length(listfastafiles)){ #for each fastafile
    filename=paste("Data/BachelerFiles/FASTAfiles/",substr(listfastafiles[i],1,6),".fasta",sep="")
    patfasta<-read.dna(filename, format = "fasta",as.character=TRUE) #read the file
    Numseqs<-c(Numseqs,nrow(patfasta))
}
#Remove from listfastafiles the files / patients with too few sequences
if (length(which(Numseqs<5))>0){listfastafiles<-listfastafiles[-which(Numseqs<5)]}

WTthreshold = 1
#0 means no threshold

#make empty vectors to count the num of muts
NumGA<-c(); NumAG<-c(); NumCT<-c(); NumTC<-c()

for (i in 1:length(listfastafiles)){ #for each fastafile
    #for (i in 1:2){ #for each fastafile
    filename=paste("Data/BachelerFiles/FASTAfiles/",substr(listfastafiles[i],1,6),".fasta",sep="")
    #print(filename)
    patfasta<-read.dna(filename, format = "fasta",as.character=TRUE) #read the file
    
    # find which seqs are from first day of sampling
    days=sort(unique(as.numeric(substr(row.names(patfasta),5,7))))
    day0sequences<-which(as.numeric(substr(row.names(patfasta),5,7))==days[1])
    latersequences<-which(as.numeric(substr(row.names(patfasta),5,7))>days[1])
    allsequences=c(day0sequences,latersequences)
    
    for ( s in allsequences){
        NumGA<-c(NumGA,length(which(consensusB[1:984]=="g"&patfasta[s,1:984]=="a")))
        NumAG<-c(NumAG,length(which(consensusB[1:984]=="a"&patfasta[s,1:984]=="g")))
        NumCT<-c(NumCT,length(which(consensusB[1:984]=="c"&patfasta[s,1:984]=="t")))
        NumTC<-c(NumTC,length(which(consensusB[1:984]=="t"&patfasta[s,1:984]=="c")))
    }         
}

print(paste("Num A to G muts per genome: ", "mean: ", round(mean(NumAG),2),"var: ",round(var(NumAG),2), "var/mean ratio: ",round(var(NumAG)/mean(NumAG),3)))
print(paste("Num G to A muts per genome: ", "mean: ",round(mean(NumGA),2),"var: ",round(var(NumGA),2), "var/mean ratio: ",round(var(NumGA)/mean(NumGA),3)))
print(paste("Num C to T muts per genome: ", "mean: ",round(mean(NumCT),2),"var: ",round(var(NumCT),2), "var/mean ratio: ",round(var(NumCT)/mean(NumCT),3)))
print(paste("Num T to C muts per genome: ", "mean: ",round(mean(NumTC),2),"var: ",round(var(NumTC),2), "var/mean ratio: ",round(var(NumTC)/mean(NumTC),3)))
