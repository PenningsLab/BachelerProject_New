#according to Johnson2010TopHIVMed.pdf

NNRTImuts<-data.frame("pos"=numeric(9),"wt"=numeric(9),"mut"=numeric(9))
NNRTImuts$pos<-c(100,101,103,106,108,181,188,190,225)
NNRTImuts$mut<-c("I","P","N","MA","I","CI","LCH","SA","H")
NRTImuts<-data.frame("pos"=numeric(11),"wt"=numeric(11),"mut"=numeric(11))
NRTImuts$pos<-c(41, 62, 67, 70, 75, 77, 116, 151, 210, 215, 219 )
NRTImuts$mut<-c("L", "V", "N", "R", "I", "L", "Y", "M", "W", "YF", "QE")
Lamimuts<-data.frame("pos"=numeric(2),"wt"=numeric(2),"mut"=numeric(2))
Lamimuts$pos<-c(65,184)#I had a typo here (67 in stead of 65, fixed Dec 2012)
Lamimuts$mut<-c("R","VI")
Indinavirmuts<-data.frame("pos"=numeric(3),"wt"=numeric(3),"mut"=numeric(3))
Indinavirmuts$pos<-c(46,82,84)
Indinavirmuts$mut<-c("IL","AFT","V")
NumImportantMuts=9+11+2+3
PositionsRT<-sort(c(NNRTImuts$pos,NRTImuts$pos,Lamimuts$pos))
PositionsPRO<-Indinavirmuts$pos

AllNtPositionsInvolvedInResistance<-sort(c(c((PositionsRT+99),PositionsPRO)*3,c((PositionsRT+99),PositionsPRO)*3-1,c((PositionsRT+99),PositionsPRO)*3-2))

RTImuts<-rbind(NNRTImuts,NRTImuts,Lamimuts)
RTImuts$mut[which(RTImuts$pos==67)]<-"NR"



