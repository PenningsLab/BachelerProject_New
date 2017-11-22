#This R script takes RawA5095Results.xls , which has the data I got from Roger Paredes, it writes MeanFreqData.csv which has the mean freq (from 2 quasi-independent measurements) of K103AAT K103AAC and Y181C in patients from the random subcohort from A5095.

###########################
#Get list of PIDs of inds from the random subcohort
###########################
read.table("../RScriptsToLookAtData/DataJonLi.csv",sep="\t",header=TRUE)->realdata
dataPare<-realdata[realdata$Study==8,]
###there are 280 patients in this list. 
dataPareRan<-dataPare[which(!is.na(dataPare$VF2)),]
ListRanPID<-as.character(dataPareRan[,2])
###there are 172 inds in this list

##########################
#Read data from Roger Paredes file ("RawA5095Results.xls")
#which has three sheets, for three mutations. 
##########################

mutations<-c("Y181","K103_AAC","K103_AAT")
library(xlsx);

for (m in 1:3){
	read.xlsx("RawA5095Results.xls",sheetIndex=m,header=TRUE)->Mutdata
	names(Mutdata)[1:3]<-c("VL","TestNum","PATID")
	x=length(unique(Mutdata$PATID))
	read.xlsx("RawA5095Results.xls",sheetIndex=m,header=TRUE,rowIndex=2:x)->Mutdata
	names(Mutdata)[1:3]<-c("VL","TestNum","PATID")
	
####################
#Clean up Mutdata
####################
	names(Mutdata)[which(names(Mutdata)=="X1.1.")]<-"Perc1"
	names(Mutdata)[which(names(Mutdata)=="X2.2.")]<-"Perc2"
	names(Mutdata)[which(names(Mutdata)=="Copy..1")]<-"n1"
	names(Mutdata)[which(names(Mutdata)=="Copy..2")]<-"n2"
	names(Mutdata)[which(names(Mutdata)=="Copy..1.1")]<-"b1"
	names(Mutdata)[which(names(Mutdata)=="Copy..2.1")]<-"b2"
	
#there rare 330 patients now. 	
	Mutdata$Freq1<-Mutdata$Perc1/100
	Mutdata$Freq2<-Mutdata$Perc2/100
	Mutdata$VarFreqs<-0;Mutdata$DiffInFreq<-0; Mutdata$MeanFreq<-0
	Mutdata$EstPartNum<- 0
	for (i in 1:length(Mutdata$EstPartNum)){
		Mutdata$EstPartNum[i]=mean(c(Mutdata$Freq1[i],Mutdata$Freq2[i]))*(1-mean(c(Mutdata$Freq1[i],Mutdata$Freq2[i])))/var(c(Mutdata$Freq1[i],Mutdata$Freq2[i]))
		Mutdata$VarFreqs[i]=var(c(Mutdata$Freq1[i],Mutdata$Freq2[i]))
		Mutdata$DiffInFreq[i]=abs(Mutdata$Freq1[i]-Mutdata$Freq2[i])
		Mutdata$MeanFreq[i]=(Mutdata$Freq1[i]+Mutdata$Freq2[i])/2
	}
	
	if (m==1)MeanFreqData<-data.frame(PATID=Mutdata$PATID,VL=Mutdata$VL,K103AACFreq=Mutdata$MeanFreq)
	if (m==2){
		for (pat in MeanFreqData$PATID){
			MeanFreqData$K103AATFreq[which(MeanFreqData$PATID==pat)]<-Mutdata$MeanFreq[which(Mutdata$PATID==pat)]}}
	if (m==3){
		for (pat in MeanFreqData$PATID){
			MeanFreqData$Y181CFreq[which(MeanFreqData$PATID==pat)]<-Mutdata$MeanFreq[which(Mutdata$PATID==pat)]}}
}

MeanFreqData$CtoTFreq<-MeanFreqData$K103AACFreq/(MeanFreqData$K103AACFreq+MeanFreqData$K103AATFreq)

write.csv(MeanFreqData,file="MeanFreqData.csv")