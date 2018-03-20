setwd("/Users/pleuni/Documents/Git/bachelerProject")
source("Rscripts/baseRscript.R")

OverviewDFBacheler <- read.table("Output/OverviewSelCoeff_BachelerFilter.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFLehman <- read.table("Output/OverviewSelCoeffLehman.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFZanini <- read.table("Output/OverviewSelCoeffZanini.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFZanini$makesCpG<-OverviewDFBacheler$makesCpG
OverviewDFLehman$makesCpG<-OverviewDFBacheler$makesCpG
OverviewDFLehman$EstSelCoeffZan<-0
    
    for (i in 1:984){
        if (is.na(OverviewDFLehman$colMeansTsLehman[i]))OverviewDFLehman$EstSelCoeffZan[i]<-NA
        if (!is.na(OverviewDFLehman$colMeansTsLehman[i]))
            OverviewDFLehman$EstSelCoeffZan[i] = EstimatedS(OverviewDFBacheler$TSmutZan[i],OverviewDFLehman$colMeansTsLehman[i])
    }

    
DF<-data.frame(dat.file=character(),typeofsite=character(),MutRates=character(),nuc=character(),
               effectsize=double(),pvalue=double(),stringsAsFactors=FALSE)
n=1

for(dat.file in c("Lehman", "Zanini", "Bacheler")){
    if (dat.file == "Lehman") dat = OverviewDFLehman
    if (dat.file == "Bacheler") dat = OverviewDFBacheler
    if (dat.file == "Zanini") dat = OverviewDFZanini
    
    for (typeofsite in c("syn", "nonsyn")){
            
        for (MutRates in c("Abram","Zan")){
            if (MutRates == "Abram") selcoeffcolumn<-which(names(dat)=="EstSelCoeff")
            if (MutRates == "Zan") selcoeffcolumn<-which(names(dat)=="EstSelCoeffZan")
        
            for (nuc in c("t","c","g")){
            
            comp="a"
        
            group1<-dat[dat$TypeOfSite==typeofsite&dat$WTnt==nuc&dat$makesCpG==0,selcoeffcolumn]
            groupComp<-dat[dat$TypeOfSite==typeofsite&dat$WTnt==comp&dat$makesCpG==0,selcoeffcolumn]
            
            #print(paste("DATA:",dat.file))
            #print(paste("mutrate:",MutRates))
            #print(paste(nuc,typeofsite))
            pvalue<-wilcox.test(
                group1,groupComp,alternative = "greater", paired = FALSE)$p.value
            
            effectsize<-round(median(group1)/median(groupComp),3)
            
            DF[n,]<-list(dat.file,typeofsite,MutRates,nuc,effectsize,pvalue)
            n=n+1
            }   
        }
    }
}

pdf("EffectSizes_Datasources.pdf")
par(mfrow=c(2,2))
for (MutRates in c("Abram","Zan")){
    for (typeofsite in c("syn","nonsyn")){
        plot(1,2,xlim=c(1.5,4.5),ylim=c(0.5,14),log="y",
             main=paste("Effect size",typeofsite),
             type="n",xlab="Nuc",ylab="effect",xaxt="n")
        legend(x = 3.5,y=5,legend = c("Zan", "Bach"),col = 2:3,pch = 16,title="Data source")
        mtext(1,at=c(2,3),text = c("C","G"))
        mtext(3,at=2.5,text=paste("mut rates from:",MutRates),cex=0.8)
        abline(h=1)
        for(i in c(2,3)){
            dataname=   c("Lehman", "Zanini", "Bacheler")[i]
            for (j in 2:3){
                pointsize=1
                if(DF$pvalue[DF$dat.file==dataname&DF$typeofsite==typeofsite&DF$MutRates==MutRates][j]<0.05) pointsize=3
                points((j)+(i-2)*0.07,DF$effectsize[DF$dat.file==dataname&DF$typeofsite==typeofsite&DF$MutRates==MutRates][j],
                       pch=16,col=i,cex=pointsize)
            }
        }
    }
}
dev.off()
