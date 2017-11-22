Bach.dat<-read.csv("Output/OverviewSelCoeff_BachelerFilter.csv")

par(mfrow=c(1,1))

head(Bach.dat)

listAA<-unique(Bach.dat$WTAA)
listfrac5<-c()

for (a in listAA){
    #get the fraction of non syn mutations that have cost > 5% for this AA
#    selco<-Bach.dat$EstSelCoeff[Bach.dat$WTAA==a& Bach.dat$TypeOfSite=="nonsyn"]
    selco<-Bach.dat$EstSelCoeff[Bach.dat$WTAA==a]
    print (a)
    print (Bach.dat$bigAAChange[Bach.dat$WTAA==a])
    listfrac5<-c(listfrac5,length(which(selco>0.05))/length(selco))
}

plot(listfrac5, pch=as.character(listAA))
points(listfrac5, cex=3,col=2)

