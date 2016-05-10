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
freqPatTs_verion2<-data.frame(row.names=substr(listfastafiles,1,6))
freqPatTs_verion3<-data.frame(row.names=substr(listfastafiles,1,6))


```{r}
ListPatientsWithoutData<-c()
for (i in 1:length(listfastafiles)){ #for each fastafile 
filename=paste("../Data/BachelerFiles/FASTAfiles/",substr(listfastafiles[i],1,6),".fasta",sep="")   
	print(filename)
	patfasta<-read.dna(filename, format = "fasta",as.character=TRUE) #read the file 
	
	if(nrow(patfasta)==1){
	freqPatTs_verion2[i,1:984]<-NA
	  }
	else{
	#which seqs are from first day of sampling
	days=sort(unique(as.numeric(substr(row.names(patfasta),5,7))))
	day0sequences<-which(as.numeric(substr(row.names(patfasta),5,7))==days[1])
	latersequences<-which(as.numeric(substr(row.names(patfasta),5,7))>days[1])
	allsequences=c(day0sequences,latersequences)


#option 1: considere positions seperately  (traditional manner)
	
	for (j in 1:984){#for each site in the sequence
 #print(paste('j',j))
	  		  WT=	consensusB[j] #what is WT at site j?
		  if(length(grep(FALSE,patfasta[day0sequences,j]==WT))>0){  # if the nucleotide of the first sequenceS is different from the consensus, then skip that patient 
	  	    #freqPatSite[i,j]='NA'
		      freqPatTs_verion2[i,j]='NA'
		    }
		  else { #check wether the neigboring sequences are the same 
		    # in order to check whether j is the first, second or third position of the codon, you can do 
		     if(j %in% seq(1,982,by=3)) {# first position
		            goodsequences<-which(paste(patfasta[,j+1],patfasta[,j+2]) == paste(consensusB[c(j+1)],consensusB[(j+2)]))
		            if (WT=="c"){freqPatTs_verion2[i,j]=length(which(patfasta[goodsequences,j]=="t"))/length(goodsequences)}  
		            if (WT=="t"){freqPatTs_verion2[i,j]=length(which(patfasta[goodsequences,j]=="c"))/length(goodsequences)}
		            if (WT=="a"){freqPatTs_verion2[i,j]=length(which(patfasta[goodsequences,j]=="g"))/length(goodsequences)}
		            if (WT=="g"){freqPatTs_verion2[i,j]=length(which(patfasta[goodsequences,j]=="a"))/length(goodsequences)}}

		     if((j %in% seq(2,983,by=3))) {# second position
		            goodsequences<-which(paste(patfasta[,j-1],patfasta[,j+1]) == paste(consensusB[c(j-1)],consensusB[(j+1)]))
		            if (WT=="c"){freqPatTs_verion2[i,j]=length(which(patfasta[goodsequences,j]=="t"))/length(goodsequences)}  
		            if (WT=="t"){freqPatTs_verion2[i,j]=length(which(patfasta[goodsequences,j]=="c"))/length(goodsequences)}
		            if (WT=="a"){freqPatTs_verion2[i,j]=length(which(patfasta[goodsequences,j]=="g"))/length(goodsequences)}
		            if (WT=="g"){freqPatTs_verion2[i,j]=length(which(patfasta[goodsequences,j]=="a"))/length(goodsequences)}}
		      
		     if((j %in% seq(3,984,by=3))) {# third position
		            goodsequences<-which(paste(patfasta[,j-2],patfasta[,j-1]) == paste(consensusB[c(j-2)],consensusB[(j-1)]))
		            if (WT=="c"){freqPatTs_verion2[i,j]=length(which(patfasta[goodsequences,j]=="t"))/length(goodsequences)}  
		            if (WT=="t"){freqPatTs_verion2[i,j]=length(which(patfasta[goodsequences,j]=="c"))/length(goodsequences)}
		            if (WT=="a"){freqPatTs_verion2[i,j]=length(which(patfasta[goodsequences,j]=="g"))/length(goodsequences)}
		            if (WT=="g"){freqPatTs_verion2[i,j]=length(which(patfasta[goodsequences,j]=="a"))/length(goodsequences)}}
		  }
}


#option 2: consider the codon as such: (Marion idea)
	for (h in 1:328){  # for each codon 

	  codonposition<-((h-1)*3+1):(h*3)
	  if(FALSE %in% names(table(consensusB[codonposition]==t(patfasta[day0sequences,codonposition])))){
	      freqPatTs_verion3[i,(((h-1)*3+1):(h*3))[1]]<-'NA'
	      freqPatTs_verion3[i,(((h-1)*3+1):(h*3))[2]]<-'NA'
	      freqPatTs_verion3[i,(((h-1)*3+1):(h*3))[3]]<-'NA'
	      }
	
	 
	  else {   #if first sequence -patient codon is equal to consensus codon, then continue  #this is ok
	        good2sequences<-which(apply(consensusB[codonposition]==(apply(patfasta[,codonposition], 1, paste)),2,sum) >1)   ## number of changes in the codon is 0 or 1 at most
	         for (cpos in codonposition){
	        		  if (consensusB[cpos]=="c"){freqPatTs_verion3[i,cpos]=length(which(patfasta[good2sequences,cpos]=="t"))/length(good2sequences)}  
		            if (consensusB[cpos]=="t"){freqPatTs_verion3[i,cpos]=length(which(patfasta[good2sequences,cpos]=="c"))/length(good2sequences)}
		            if (consensusB[cpos]=="a"){freqPatTs_verion3[i,cpos]=length(which(patfasta[good2sequences,cpos]=="g"))/length(good2sequences)}
		            if (consensusB[cpos]=="g"){freqPatTs_verion3[i,cpos]=length(which(patfasta[good2sequences,cpos]=="a"))/length(good2sequences)}
	              }
	        }    

	  

	  
	}
#remove Patients without data 
for (pat in ListPatientsWithoutData){
  freqPatSite<-freqPatSite[-which(row.names(freqPatSite)==pat),]
  freqPatTs<-freqPatTs[-which(row.names(freqPatTs)==pat),]
}
write.csv(freqPatTs_verion2,file="../Output/freqPatTsInclDay0-version2.csv")
write.csv(freqPatTs_verion3,file="../Output/freqPatTsInclDay0-version3.csv")

}	  }

#}}}   	     
```











