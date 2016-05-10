#get all the names of the fasta files for each patient in Bacheler2000 dataset
#determine all the frequencies and store in freqPatSite data.frame
#only needed if the stored data are not O


#Load libraries and necessary files from the baseRscript.Rmd
source('baseRscript.R')

#Read the correct fastafiles. 
if (FALSE){

list.files(path="../Data/BachelerFiles/FASTAfiles/")->listfastafiles

lengthlatersequenceslist<-c();lengthallsequenceslist<-c()

#make dataframe with frequencies for all non-muts for all patients for all sites. 
freqPatSite<-data.frame(row.names=substr(listfastafiles,1,6))
freqPatTs<-data.frame(row.names=substr(listfastafiles,1,6))
freqPatTv<-data.frame(row.names=substr(listfastafiles,1,6))

ListPatientsWithoutData<-c()
for (i in 1:length(listfastafiles)){ #for each fastafile 
filename=paste("../Data/BachelerFiles/FASTAfiles/",substr(listfastafiles[i],1,6),".fasta",sep="")   ## kristof: waarom?   ## kristof changer
	print(filename)
	patfasta<-read.dna(filename, format = "fasta",as.character=TRUE) #read the file 
	#which seqs are from first day of sampling
	days=sort(unique(as.numeric(substr(row.names(patfasta),5,7))))
	day0sequences<-which(as.numeric(substr(row.names(patfasta),5,7))==days[1])
	latersequences<-which(as.numeric(substr(row.names(patfasta),5,7))>days[1])
	allsequences=c(day0sequences,latersequences)
	for (j in 1:984){#for each site in the sequence
		WT=	consensusB[j] #what is WT at site j?
		freqPatSite[i,j]=length(which(patfasta[allsequences,j]!=WT))/length(allsequences)#if WT, what is freq of non WT at later time points?  ## kristof: any change, ts and tv
        if (WT=="c"){freqPatTs[i,j]=length(which(patfasta[allsequences,j]=="t"))/length(allsequences)}   ## kristof only ts
		if (WT=="t"){freqPatTs[i,j]=length(which(patfasta[allsequences,j]=="c"))/length(allsequences)}
		if (WT=="a"){freqPatTs[i,j]=length(which(patfasta[allsequences,j]=="g"))/length(allsequences)}
		if (WT=="g"){freqPatTs[i,j]=length(which(patfasta[allsequences,j]=="a"))/length(allsequences)}
        if (WT=="a" | WT=="g"){freqPatTv[i,j]=length(which(patfasta[allsequences,j]=="c" | patfasta[allsequences,j]=="t"))/length(allsequences)}
        if (WT=="c" | WT=="t"){freqPatTv[i,j]=length(which(patfasta[allsequences,j]=="a" | patfasta[allsequences,j]=="g"))/length(allsequences)}


        }
}
#remove Patients without data 
for (pat in ListPatientsWithoutData){
  freqPatSite<-freqPatSite[-which(row.names(freqPatSite)==pat),]
  freqPatTs<-freqPatTs[-which(row.names(freqPatTs)==pat),]
}
write.csv(freqPatSite,file="../Output/freqPatSiteInclDay0.csv")
write.csv(freqPatTs,file="../Output/freqPatTsInclDay0.csv")
write.csv(freqPatTv,file="../Output/freqPatTvInclDay0.csv")
}

