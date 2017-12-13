## Step 4:  R analysis for plotting frequencies: correlation with bacheler dataset 

#Get Bacheler frequencies previously calculated ** 

setwd("~/Documents/Git/bachelerProject/Comparison-stanford_vs_bacheler-april2016")
#setwd("~/Dropbox/MarionKristofBachelerProject/GitMarionKristof/bachelerProject/Comparison-stanford_vs_bacheler-april2016")
source("../Rscripts/baseRscript.R")
OverviewDFBach<-read.csv("../Output/OverviewSelCoeff_BachelerFilter.csv",row.names=1)
OverviewDFBach$StanfordFreqs<-t(read.csv("subtypeB-frequencies_stanford_pr_rt_naive.csv"))[1:length(OverviewDFBach[,1])]

OverviewDFZanini<-read.csv("../Output/OverviewSelCoeffZanini.csv",row.names=1)
OverviewDFZanini$StanfordFreqs<-t(read.csv("subtypeB-frequencies_stanford_pr_rt_naive.csv"))[1:length(OverviewDFBach[,1])]

OverviewDFLehman<-read.csv("../Output/OverviewSelCoeffLehman.csv",row.names=1)
OverviewDFLehman$StanfordFreqs<-t(read.csv("subtypeB-frequencies_stanford_pr_rt_naive.csv"))[1:length(OverviewDFBach[,1])]

# plots all comparisons

#Take only syn and nonsyn sites after 13 AA of protease
synsites<-which(OverviewDFBach$TypeOfSite=="syn"&OverviewDFBach$num>39&consensusA==consensusB&consensusC==consensusB)
nonsynsites<-which(OverviewDFBach$TypeOfSite=="nonsyn"&OverviewDFBach$num>39&consensusA==consensusB&consensusC==consensusB)
stopsites<-which(OverviewDFBach$TypeOfSite=="stop"&OverviewDFBach$num>39&consensusA==consensusB&consensusC==consensusB)

length(which(OverviewDFBach$TypeOfSite!="res"&consensusC==consensusB&consensusB==consensus01AE&OverviewDFBach$num>39))
length(which(OverviewDFBach$TypeOfSite!="res"&consensusC==consensusB&OverviewDFBach$num>39))
length(which(OverviewDFBach$TypeOfSite!="res"&OverviewDFBach$num>39))

# plots all Zanini
plot(OverviewDFZanini$colMeansTsZanini+0.0001,OverviewDFZanini$StanfordFreqs+0.0001,
     main="Within vs between-host frequencies Zanini",col=0,
     xlim=c(0.0001,1),ylim=c(0.0001,1),
     xlab="Within-host frequencies Zanini",
     ylab="Between-host frequencies Stanford",
    log="xy"
)
points(OverviewDFZanini$colMeansTsZanini[synsites]+0.0001,OverviewDFZanini$StanfordFreqs[synsites]+0.0001,
       col="lightgreen",pch=16)
points(OverviewDFZanini$colMeansTsZanini[nonsynsites]+0.0001,OverviewDFZanini$StanfordFreqs[nonsynsites]+0.0001,
       col="orange",pch=16)

CorZaniniStanford<-cor.test(OverviewDFZanini$colMeansTsZanini[c(synsites,nonsynsites)],OverviewDFZanini$StanfordFreqs[c(synsites,nonsynsites)])$estimate
text(0.05,0.3,paste("R^2 = ",round(R2ZaniniStanford,2)))

#Better, use spearman, bc not normally distributed
cor.test(OverviewDFBach$MeanFreq[c(synsites,nonsynsites)],OverviewDFBach$StanfordFreqs[c(synsites,nonsynsites)],method="s")

#Compare Zanini
cor.test(OverviewDFZanini$colMeansTsZanini[c(synsites,nonsynsites)],OverviewDFZanini$StanfordFreqs[c(synsites,nonsynsites)],method="s")

cor.test(OverviewDFZanini$colMeansTsZanini[c(synsites,nonsynsites)],OverviewDFBach$MeanFreq[c(synsites,nonsynsites)],method="s")


#Now look at f(1-f)
cor.test(OverviewDFBach$MeanFreq[c(synsites,nonsynsites)]*(1-OverviewDFBach$MeanFreq[c(synsites,nonsynsites)]),
         OverviewDFBach$StanfordFreqs[c(synsites,nonsynsites)]*(1-OverviewDFBach$StanfordFreqs[c(synsites,nonsynsites)]),method="s")
#SAME RESULT BC RANK TEST! 

#Sep for nonsyn and syn
cor.test(OverviewDFBach$MeanFreq[c(nonsynsites)],OverviewDFBach$StanfordFreqs[c(nonsynsites)],method="p")

cor.test(OverviewDFBach$MeanFreq[c(synsites)],OverviewDFBach$StanfordFreqs[c(synsites)],method="p")




# plots all Lehman
plot(OverviewDFLehman$colMeansTsLehman+0.0001,OverviewDFLehman$StanfordFreqs+0.0001,
     main="Within vs between-host frequencies Lehman",col=0,
     xlim=c(0.0001,1),ylim=c(0.0001,1),
     xlab="Within-host frequencies Lehman",
     ylab="Between-host frequencies Stanford",
     log="xy"
)
points(OverviewDFLehman$colMeansTsLehman[synsites]+0.0001,OverviewDFLehman$StanfordFreqs[synsites]+0.0001,
       col="lightgreen",pch=16)
points(OverviewDFLehman$colMeansTsLehman[nonsynsites]+0.0001,OverviewDFLehman$StanfordFreqs[nonsynsites]+0.0001,
       col="orange",pch=16)

R2LehmanStanford<-cor.test(OverviewDFLehman$colMeansTsLehman[c(synsites,nonsynsites)],OverviewDFLehman$StanfordFreqs[c(synsites,nonsynsites)])$estimate
text(0.05,0.3,paste("R^2 = ",round(R2LehmanStanford,2)))

#Plot Bacheler
plot(OverviewDFBach$MeanFreq+0.0001,OverviewDFBach$StanfordFreqs+0.0001,
     main="Within vs between-host frequencies Bacheler",col=0,
     xlim=c(0.0001,.05),ylim=c(0.0001,.2),
     xlab="Within-host frequencies Bacheler",
     ylab="Between-host frequencies Stanford",
     log="xy"
)
points(OverviewDFBach$MeanFreq[synsites]+0.0001,OverviewDFBach$StanfordFreqs[synsites]+0.0001,
       col=cols[3],pch=16)
points(OverviewDFBach$MeanFreq[nonsynsites]+0.0001,OverviewDFBach$StanfordFreqs[nonsynsites]+0.0001,
       col=cols[5],pch=16)

R2BachStanford<-cor.test(
    OverviewDFBach$MeanFreq[c(synsites,nonsynsites)],
    OverviewDFBach$StanfordFreqs[c(synsites,nonsynsites)],
    method="s")$estimate
text(0.0002,0.1,paste("R^2 = ",round(R2BachStanford,2)))

#Plot Bacheler RT only
plot(OverviewDFBach$MeanFreq+0.0001,OverviewDFBach$StanfordFreqs+0.0001,
     main="Within vs between-host frequencies Bacheler",col=0,
     xlim=c(0.0001,1),ylim=c(0.0001,1),
     xlab="Within-host frequencies Bacheler",
     ylab="Between-host frequencies Stanford",
     log="xy"
)
points(OverviewDFBach$MeanFreq[synsites[which(synsites>297)]]+0.0001,OverviewDFBach$StanfordFreqs[synsites[which(synsites>297)]]+0.0001,
       col="lightgreen",pch=16)
points(OverviewDFBach$MeanFreq[nonsynsites[which(nonsynsites>297)]]+0.0001,OverviewDFBach$StanfordFreqs[nonsynsites[which(nonsynsites>297)]]+0.0001,
       col="orange",pch=16)

R2BachStanford<-cor.test(OverviewDFBach$MeanFreq[c(synsites,nonsynsites)],OverviewDFBach$StanfordFreqs[c(synsites,nonsynsites)])$estimate
text(0.05,0.3,paste("R^2 = ",round(R2BachStanford,2)))

