#PSP: change this so that it reads the most current data at all times
dat <- read.table("GLM_Analysis/dat/BachelerCountData.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)
#dat has one row for each combination of patient and site
suppDat <- read.table("GLM_Analysis/dat/OverviewSelCoeffwProteinFeatures.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
suppDatShape <- read.table("GLM_Analysis/dat/OverviewSelCoeffwSHAPE.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

#The groups of amino acids
pos <- "R|H|K"
neg <- "D|E"
unc <- "S|T|N|Q"
spe <- "C|U|G|P"
hyd <- "A|I|L|F|M|W|Y|V"
amCat <- function(AA){
    if(regexpr(pos, AA) > 0){ return(0) }
    if(regexpr(neg, AA) > 0){ return(1) }
    if(regexpr(unc, AA) > 0){ return(2) }
    if(regexpr(spe, AA) > 0){ return(3) }
    if(regexpr(hyd, AA) > 0){ return(4) }
    return(5)
}

#Does the mutation make a CpG? 
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

#Does the mutation make a drastic AA change?
bigChange <- rep(NA, nrow(suppDat))
for(i in 1:nrow(suppDat)){
    WT <- amCat(suppDat[i,'WTAA'])
    MUT <- amCat(suppDat[i,'MUTAA'])
    if(WT == MUT){ bigChange[i] <- 0
    }else{
        bigChange[i] <- 1
    }
}

#PSP: I am not sure whath this function does. 
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

numRows <- sum(dat$WTcount) + sum(dat$MutCount)
datRows <- matrix(data = NA, nrow = numRows, ncol = length(names))
colnames(datRows) <- names

relInd <- 1
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
    match.is <- which(dat[,3 ] == i)
    for(j in unique(dat[,2])){ #for all patients
        
        toExp <- dat[intersect(which(dat[,2] == j), match.is),]
        minorNum <- toExp$MutCount
        majorNum <- toExp$WTcount - minorNum
        toReturn <- rbind(repMe(pref, 1, minorNum), repMe(pref, 0, majorNum)) #create separate row for each observation
        
        if(!is.null(toReturn)){
            datRows[relInd:(relInd + nrow(toReturn) - 1),] <- toReturn
            relInd <- relInd + nrow(toReturn)
        }
    }
}

write.csv(as.data.frame(datRows[which(!is.na(datRows[,1])),]),file = "GLM_Analysis/dat/datFitModel.csv")