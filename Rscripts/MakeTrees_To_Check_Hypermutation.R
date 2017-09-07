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

pdf("Output/trees_to_check_hypermutation.pdf")
pdf("trees.pdf")
for (i in 1:length(listfastafiles)){ #for each fastafile
    #for (i in 1:2){ #for each fastafile
    filename=paste("Data/BachelerFiles/FASTAfiles/",substr(listfastafiles[i],1,6),".fasta",sep="")
    print(filename)
    patfasta_alignment<-read.alignment(file=filename, format = "fasta") #read the file
    XYZ <- dist.alignment(patfasta_alignment)
    tree<-nj(XYZ)
    plot(tree,main=filename) 
}
dev.off()
