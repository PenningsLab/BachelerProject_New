#PSP: Need to change this so that it reads the most current data at all times

#Use filtered data May 2017
#dat <- read.table("Output/BachelerCountData_filter.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)
dat <- read.table("Output/BachelerCountData.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)
suppDat <- read.table("GLM_Analysis/dat/OverviewSelCoeffwProteinFeatures.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
suppDatShape <- read.table("GLM_Analysis/dat/OverviewSelCoeffwSHAPE.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)


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
bigChange <- rep(NA, nrow(suppDat))
for(i in 1:nrow(suppDat)){
	WT <- amCat(suppDat[i,'WTAA'])
	MUT <- amCat(suppDat[i,'MUTAA'])
	if(WT == MUT){ bigChange[i] <- 0
	}else{
		bigChange[i] <- 1
	}
}


#Which positions will create CpG sites if they undergo transitions
makesCpG <- rep(0, length = nrow(suppDat))
for(i in 1:nrow(suppDat)){
	trip <- suppDat$WTnt[c(i-1, i, i + 1)]
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
for(i in 40:length(suppDat$num)){
	
	print(i)
	pos <- suppDat[i,]$num
	inRT <- as.numeric(pos < 298)
	atcg <- c(0,0,0,0)
	atcg[which(nucord == suppDat[i,]$WTnt)] <- 1
	bigAAChange <- bigChange[i]
	nonsyn <- as.numeric(regexpr("nonsyn",suppDat[i,]$TypeOfSite) > 0)
	stop <- as.numeric(regexpr("stop",suppDat[i,]$TypeOfSite) > 0)
	res <- as.numeric(regexpr("res",suppDat[i,]$TypeOfSite) > 0)
	pairingprob <- suppDatShape[i,]$PAIRINGPROB
	shape <- suppDatShape[i,]$SHAPE
	CpG <- makesCpG[i]
	proteinfeats <- c(0,0,0,0)
	proteinfeats[which(suppDat[i,]$Protein_feature == c("3-10-helix", "alpha", "beta", "coil"))] <- 1
	pref <- c(atcg, inRT, bigAAChange, nonsyn, stop, res, shape, pairingprob, CpG, proteinfeats)
	
	#pref is the same for all patients.
	#Now, let's go through each of the patients and do this:
	if(res != 1){
		match.is <- which(dat[,3 ] == i)
		for(j in unique(dat[,2])){
			
			toExp <- dat[intersect(which(dat[,2] == j), match.is),] #PLEUNI WHAT IF THIS DOESN"T EXIST??
			toReturn<-c()
			if (length(toExp$X)){
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

#Pleuni changed this to write to results_filtered_data/
#write.csv(as.data.frame(datRows[which(!is.na(datRows[,1])),]),file = "GLM_Analysis/results_filtered_data/datFitModel.csv")
write.csv(as.data.frame(datRows[which(!is.na(datRows[,1])),]),file = "GLM_Analysis/dat/datFitModel.csv")

}

read.csv("GLM_Analysis/dat/datFitModel.csv",row.names=1)->datFitModel
#Pleuni changed this to read from results_filtered_data/
#read.csv("GLM_Analysis/results_filtered_data/datFitModel.csv",row.names=1)->datFitModel
