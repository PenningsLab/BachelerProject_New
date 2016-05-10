#Will read Zanini files, determine where POL is. 

#This is just for one patient. Later add others.
ZaniniFiles<-list.files("..//Data/ZaniniNeherData/",recursive = TRUE,pattern="tsv")


#SeqData<-read.csv(paste("..//Data/ZaniniNeherData/",ZaniniFiles[1],sep=""),sep="\t",skip=1)
#SeqData$MajNt<-""
#SeqData$MajNt<-apply(X[,1:4],1,function(x) c("a","c","g","t")[which.max(x)])

#i=0
#while(i , 10000){
#    print(i)
#    polstart=regexpr("cctca",paste(X$MajNt[i:9060],collapse="")) [[1]]
#    print(polstart)
#    print(seqinr::translate(X$MajNt[(i+polstart-1):(i+polstart+100)]))
#    i=i+polstart+1
#}

#start of POL is 1719 :-)
#seqinr::translate(X$MajNt[(1719):(1821)])
#seqinr::translate(X$MajNt[(1719+297):(1821+297)])

Prostart<-data.frame(Pat=paste("act_p",c(1:3,5:6,8:11),sep=""))
Prostart$start<-c(1719,1714,1728,1706,1687, 1705,1698, 1707,1713)

    #make dataframe with frequencies for all non-muts for all patients for all sites filtered with the WT threshold.
freqPatTs_Zanini<-data.frame(row.names=ZaniniFiles)

for (i in 1:length(ZaniniFiles)){
    print(i)
    print(ZaniniFiles[i])
    SeqData<-read.csv(paste("..//Data/ZaniniNeherData/",ZaniniFiles[i],sep=""),sep="\t",skip=1)
    #Pol starts at 1719. First let's just take 100 nt. Later 984.
    pat = substr(ZaniniFiles[i],1,regexpr("/",ZaniniFiles[i])[[1]]-1)
    print(pat)
    SeqData<-SeqData[Prostart$start[which(Prostart$Pat==pat)]:(Prostart$start[which(Prostart$Pat==pat)]+984-1),]
    SeqData$MajNt<-""
    SeqData$MajNt<-apply(SeqData[,1:4],1,function(x) c("a","c","g","t")[which.max(x)])
    print(seqinr::translate(SeqData$MajNt[1:30]))
    print(seqinr::translate(SeqData$MajNt[298:(298+29)]))  
    #What is transition mut?
    SeqData$consensusB<-consensusB[1:length(SeqData[,1])]
    for (j in 1:length(SeqData$consensusB)) SeqData$transition[j]<-transition(SeqData$consensusB[j])
    
    #determine Ts freq of every site. 
    SeqData$freq<-0
    for (k in 1:length(SeqData$consensusB)){#for each site in the sequence
        MutNum<- SeqData [k,which(c("a","c","g","t")==SeqData$transition[k])]
        WTNum <- SeqData [k,which(c("a","c","g","t")==SeqData$consensusB[k])]
        #check wether the neigboring sequences are the same / WE CANT DO THIS FOR THE ZANINI DATA
        SeqData$freq[k]<-MutNum/(MutNum+WTNum) 
        freqPatTs_Zanini[i,k]<-SeqData$freq[k]
    }
}

write.csv(freqPatTs_Zanini,file="../Output/freqPatTs_Zanini.csv")    

pdf("TRY.pdf")
par(mfrow=c(3,3))
for (i in 1:100) hist(freqPatTs_Zanini[,i],main=i,xlim=c(0,1),breaks=seq(0,1,by=0.01),col=2)
dev.off()



