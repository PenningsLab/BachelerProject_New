Maoz.dat <- read.table("../dat/kaks.edited.res", sep = "\t", stringsAsFactors = FALSE, header = FALSE)
Bach.dat <- read.table("../dat/OverviewSelCoeffwSHAPE.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
Lehman.dat <- read.table("../dat/OverviewSelCoeffLehman.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
Zanini.dat <- read.table("../dat/OverviewSelCoeffZanini.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)




Maoz.dat[,1] == Bach.dat[1:99,2]

table(Maoz.dat[,3])
cor.test(Maoz.dat[okInds,3], log(AAsels))

dim(Bach.dat)
Maoz.dat[20,]
Bach.dat[1:3,]

AAsels <- c()
for(i in 1:99){
    relinds <- (3*i -2):(3*i)
    nsinds <- which(Bach.dat[relinds, ]$TypeOfSite  == "nonsyn")
    AAsels[i]<-  mean(Bach.dat[relinds[nsinds],]$EstSelCoeff)
}


okInds <- intersect(which(Maoz.dat[,3] < 3), which(!is.na(AAsels)))
plot(Maoz.dat[okInds,3], log(AAsels)[okInds], ylim = c(-10, 0), xlim = c(0, 1), pch = 16, col = rgb(0,0,0,.5))
cor.test(Maoz.dat[okInds,3], log(AAsels)[okInds])
abline(lm(log(AAsels)[okInds] ~ Maoz.dat[okInds,3]), col = cols[2])


Maoz.All <- c()
AAs.All <- c()
for(i in 1:99){
    relinds <- (3*i -2):(3*i)
    nsinds <- which(Bach.dat[relinds, ]$TypeOfSite  == "nonsyn")
    for(j in 1:length(nsinds)){
        AAs.All <- c(AAs.All, Bach.dat[relinds[nsinds[j]],]$EstSelCoeff)
        Maoz.All <- c(Maoz.All, Maoz.dat[i,3])
    }
}

okInds <- intersect(which(Maoz.All < 3), which(!is.na(AAs.All)))
plot(Maoz.All[okInds], log(AAs.All)[okInds], ylim = c(-10, 0), xlim = c(0, 1), pch = 16, col = rgb(0,0,0,.5))
abline(lm(log(AAs.All)[okInds] ~ Maoz.All[okInds]), col = cols[2])
cor.test(AAs.All[okInds], Maoz.All[okInds])

Maoz.dat[1,3]
cbind(AmAcid = rep(1,3), Bach.dat[1:3,c("TypeOfSite", "WTAA", "MUTAA", "EstSelCoeff")])


