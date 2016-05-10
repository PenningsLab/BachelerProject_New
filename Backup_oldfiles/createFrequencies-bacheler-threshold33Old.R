#new addition
#1. only take patienys with a consensus sequence at Day0
#2. only consider changed nucleotides when the two other nucleotides in a codon are not mutated

#get all the names of the fasta files for each patient in Bacheler2000 dataset
#determine all the frequencies and store in freqPatSite data.frame
#only needed if the stored data are not O


#Load libraries and necessary files from the baseRscript.Rmd
source('baseRscript.R')

#Read the correct fastafiles. 
if (TRUE){

list.files(path="../Data/BachelerFiles/FASTAfiles/")->listfastafiles

lengthlatersequenceslist<-c();lengthallsequenceslist<-c()

#make dataframe with frequencies for all non-muts for all patients for all sites. 
freqPatTs_threshold<-data.frame(row.names=substr(listfastafiles,1,6))
#make dataframe that keeps track, for each position and patients, what happened according to criteria and which data points are included
filterSummary<-data.frame(row.names=substr(listfastafiles,1,6)) 

Patientswith1sequence<-0
Nonconsensusday0_pat_pos<-c()
ListPatientsWithoutData<-c()

for (i in 1:length(listfastafiles)){ #for each fastafile 
filename=paste("../Data/BachelerFiles/FASTAfiles/",substr(listfastafiles[i],1,6),".fasta",sep="")   
	print(filename)
	patfasta<-read.dna(filename, format = "fasta",as.character=TRUE) #read the file 
	
	if(nrow(patfasta)==1){
	freqPatTs_threshold[i,1:984]<-NA
	filterSummary[i,1:984]<-'1seq'
	Patientswith1sequence<-Patientswith1sequence+1  
	}
	else{
	#which seqs are from first day of sampling
	days=sort(unique(as.numeric(substr(row.names(patfasta),5,7))))
	day0sequences<-which(as.numeric(substr(row.names(patfasta),5,7))==days[1])
	latersequences<-which(as.numeric(substr(row.names(patfasta),5,7))>days[1])
	allsequences=c(day0sequences,latersequences)


	for (j in 1:984){#for each site in the sequence
 #print(paste('j',j))
	  		  WT=	consensusB[j] #what is WT at site j?
	  		  
	  		  if(length(names(prop.table(table(patfasta[day0sequences,j] ==WT))))==1 & names(prop.table(table(patfasta[day0sequences,j] ==WT)))==FALSE){
	  		    freqPatTs_threshold[i,j]<-NA
	  		    #keep track of pat and pos that are excluded
	  		    Nonconsensusday0_pat_pos<-rbind(Nonconsensusday0_pat_pos,c(i,j))
	  		    filterSummary[i,j]<-'nonconsensus'
	  		  }
	  		  else if(prop.table(table(patfasta[day0sequences,j]==WT))[which(names(prop.table(table(patfasta[day0sequences,j] ==WT)))==TRUE)] < 0.66){	  
		      freqPatTs_threshold[i,j]<-NA
		      #keep track of pat and pos that are excluded
	  		  Nonconsensusday0_pat_pos<-rbind(Nonconsensusday0_pat_pos,c(i,j))
	  		  filterSummary[i,j]<-'nonconsensus'
	  		    }
		  else { #check wether the neigboring sequences are the same 
		    # in order to check whether j is the first, second or third position of the codon, you can do 
		     if(j %in% seq(1,982,by=3)) {# first position
		        goodsequences<-which(paste(patfasta[,j+1],patfasta[,j+2]) == paste(consensusB[c(j+1)],consensusB[(j+2)]))
		        filterSummary[i,j]<-length(goodsequences)
		        if (WT=="c"){freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]=="t"))/length(goodsequences)}
		        if (WT=="t"){freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]=="c"))/length(goodsequences)}
		        if (WT=="a"){freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]=="g"))/length(goodsequences)}
		        if (WT=="g"){freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]=="a"))/length(goodsequences)
		        }}

		     if((j %in% seq(2,983,by=3))) {# second position
		        goodsequences<-which(paste(patfasta[,j-1],patfasta[,j+1]) == paste(consensusB[c(j-1)],consensusB[(j+1)]))
		        filterSummary[i,j]<-length(goodsequences)
		        if (WT=="c"){freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]=="t"))/length(goodsequences)}
		        if (WT=="t"){freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]=="c"))/length(goodsequences)}
		        if (WT=="a"){freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]=="g"))/length(goodsequences)}
		        if (WT=="g"){freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]=="a"))/length(goodsequences)
		        }}
		      
		     if((j %in% seq(3,984,by=3))) {# third position
		         goodsequences<-which(paste(patfasta[,j-2],patfasta[,j-1]) == paste(consensusB[c(j-2)],consensusB[(j-1)]))
		         filterSummary[i,j]<-length(goodsequences)
		         if (WT=="c"){freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]=="t"))/length(goodsequences)}
		        if (WT=="t"){freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]=="c"))/length(goodsequences)}
		        if (WT=="a"){freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]=="g"))/length(goodsequences)}
		        if (WT=="g"){freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]=="a"))/length(goodsequences)}}
		        }}



	}
#remove Patients without data 
for (pat in ListPatientsWithoutData){
  freqPatSite<-freqPatSite[-which(row.names(freqPatSite)==pat),]
  freqPatTs_threshold<-freqPatTs_threshold[-which(row.names(freqPatTs)==pat),]
}
write.csv(freqPatTs_threshold,file="../Output/freqPatTsInclDay0-threshold.csv")
}	  



}












