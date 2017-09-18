#PSP: Sept 2017 Changed this so that it reads the most current data at all times

setwd("~/Documents/Git/bachelerProject/Rscripts")
#Use filtered data May 2017
dat <- read.table("../Output/BachelerCountData_Threshold1.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)
#suppDat <- read.table("../GLM_Analysis/dat/OverviewSelCoeffwProteinFeatures.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
#suppDatShape <- read.table("GLM_Analysis/dat/OverviewSelCoeffwSHAPE.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDF<-read.csv("../Output/OverviewSelCoeff_BachelerFilter.csv")
PolShapeData<-read.csv("../Data/Pol_SHAPE.csv")

#function to judge amino acid categories
pos <- "R|H|K"
neg <- "D|E"
unc <- "S|T|N|Q"
spe <- "C|U|G|P"
hyd <- "A|I|L|F|M|W|Y|V"

amCat <- function(AA){
	if(regexpr(pos, AA) > 0){ return(0) }
	if(regexpr(neg, AA) > 0){ return(1) }
	if(regexpr(unc, AA) > 0){ return(2) }
	if(regexpr(hyd, AA) > 0){ return(3) }
	if(regexpr("C", AA) > 0){ return(4) }
	if(regexpr("U", AA) > 0){ return(6) }
	if(regexpr("G", AA) > 0){ return(7) }
	if(regexpr("P", AA) > 0){ return(8) }
	return(5)
}

#Which positions will have big amino acid changes if they undergo transitions
bigChange <- rep(NA, nrow(OverviewDF))
for(i in 1:nrow(OverviewDF)){
	WT <- amCat(OverviewDF[i,'WTAA'])
	MUT <- amCat(OverviewDF[i,'MUTAA'])
	if(WT == MUT){ bigChange[i] <- 0
	}else{
		bigChange[i] <- 1
	}
}


#Which positions will create CpG sites if they undergo transitions
makesCpG <- rep(0, length = nrow(OverviewDF))
for(i in 1:nrow(OverviewDF)){
	trip <- OverviewDF$WTnt[c(i-1, i, i + 1)]
	if(trip[1] == "c" & trip[2] == "a" ){
		makesCpG[i] <- 1
	}
	if(trip[2] == "t" & trip[3] == "g"){
		makesCpG[i] <- 1
	}
}

#Helper function to do some matrix replication in creating the dataset
repMe <- function(pref, majormin, toRep){
	if(toRep == 0){ return() }
	toRet <- matrix(data = NA, nrow = toRep, ncol = length(pref) + 1)
	for(i in 1:toRep){
		toRet[i,] <- c(pref, majormin)
	}
	return(toRet)
}

names = c( "a", "t", "c", "g", "inRT", "bigAAChange", "nonsyn", "stop", "res", "shape", "pairingprob", "CpG", "helix", "alpha", "beta", "coil", "minor")
nucord <- c("a", "t", "c", "g")

numRows <- sum(dat$Totalcount) + sum(dat$MutCount)
datRows <- matrix(data = NA, nrow = numRows, ncol = length(names))
colnames(datRows) <- names

dim(datRows)
#3355728

#This part of dataprep only needs to be run once. 
if (FALSE){ #THIS ONLY NEEDS TO RUN ONCE, TO CREATE "GLM_Analysis/dat/datFitModel.csv"
relInd <- 1
#Leave out the first 40 positions
#for(i in 40:100){
for(i in 40:length(OverviewDF$num)){
        
	print(i)
	pos <- OverviewDF[i,]$num
	inRT <- as.numeric(pos < 298)
	atcg <- c(0,0,0,0)
	atcg[which(nucord == OverviewDF[i,]$WTnt)] <- 1
	bigAAChange <- bigChange[i]
	nonsyn <- as.numeric(regexpr("nonsyn",OverviewDF[i,]$TypeOfSite) > 0)
	stop <- as.numeric(regexpr("stop",OverviewDF[i,]$TypeOfSite) > 0)
	res <- as.numeric(regexpr("res",OverviewDF[i,]$TypeOfSite) > 0)
	pairingprob <- PolShapeData[i,]$PAIRINGPROB
	shape <- PolShapeData[i,]$SHAPE
	CpG <- makesCpG[i]
	proteinfeats <- c(0,0,0,0)
	proteinfeats[which(PolShapeData[i,]$Protein_feature == c("3-10-helix", "alpha", "beta", "coil"))] <- 1
	pref <- c(atcg, inRT, bigAAChange, nonsyn, stop, res, shape, pairingprob, CpG, proteinfeats)
	
	#pref is the same for all patients.
	#Now, let's go through each of the patients and do this:
	if(res != 1){
		match.is <- which(dat[,3 ] == i) #Which of the count data is for position i? This is usually around 160 lins, one for each patient
		for(j in unique(dat[,2])){#for each patient
			toExp <- dat[intersect(which(dat[,2] == j), match.is),] #get the info for that position for that patient
			toReturn<-c() 
			if (length(toExp$X)){  #Only happens if the patient has data for that position
			minorNum <- toExp$MutCount
			majorNum <- toExp$Totalcount - minorNum #Pleuni 05/17 WTcount is actually the number of good sequences. So this is correct. 
			toReturn <- rbind(repMe(pref, 1, minorNum), repMe(pref, 0, majorNum))}
			if(!is.null(toReturn)){
				datRows[relInd:(relInd + nrow(toReturn) - 1),] <- toReturn
				relInd <- relInd + nrow(toReturn)
			}
		}
	}
}

write.csv(as.data.frame(datRows[which(!is.na(datRows[,1])),]),file = "../Output/datFitModel.csv")

}
read.csv("../Output/datFitModel.csv",row.names=1)->datFitModel
