## Step 4:  R analysis for plotting frequencies: correlation with bacheler dataset 

#Get Bacheler frequencies previously calculated ** 

setwd("~/Dropbox/MarionKristofBachelerProject/GitMarionKristof/bachelerProject/Comparison-stanford_vs_bacheler-april2016")
source("../Rscripts/baseRscript.R")
OverviewDFBach<-read.csv("../Output/OverviewSelCoeff_Bacheler.csv",row.names=1)
OverviewDFBach$StanfordFreqs<-t(read.csv("subtypeB-frequencies_stanford_pr_rt_naive.csv"))[1:length(OverviewDFBach[,1])]

OverviewDFZanini<-read.csv("../Output/OverviewSelCoeffZanini.csv",row.names=1)
OverviewDFZanini$StanfordFreqs<-t(read.csv("subtypeB-frequencies_stanford_pr_rt_naive.csv"))[1:length(OverviewDFBach[,1])]

OverviewDFLehman<-read.csv("../Output/OverviewSelCoeffLehman.csv",row.names=1)
OverviewDFLehman$StanfordFreqs<-t(read.csv("subtypeB-frequencies_stanford_pr_rt_naive.csv"))[1:length(OverviewDFBach[,1])]

which(consensusA==consensusB& consensusC==consensusB)

# plots all Bacheler
#Take only syn and nonsyn sites after 13 AA of protease
synsites<-which(OverviewDFBach$TypeOfSite=="syn"&OverviewDFBach$num>39&consensusA==consensusB&consensusC==consensusB)
nonsynsites<-which(OverviewDFBach$TypeOfSite=="nonsyn"&OverviewDFBach$num>39&consensusA==consensusB&consensusC==consensusB)

plot(OverviewDFBach$colMeansTs0,OverviewDFBach$StanfordFreqs,
     main="Within vs between-host frequencies",col=0,
     xlim=c(0,.4),ylim=c(0,.4),
     xlab="Within-host frequencies Bacheler",
     ylab="Between-host frequencies Stanford"
     )
points(OverviewDFBach$colMeansTs0[synsites],OverviewDFBach$StanfordFreqs[synsites],
       col="lightgreen",pch=16)
points(OverviewDFBach$colMeansTs0[nonsynsites],OverviewDFBach$StanfordFreqs[nonsynsites],
       col="orange",pch=16)

R2BachStanford<-cor.test(OverviewDFBach$colMeansTs0[c(synsites,nonsynsites)],OverviewDFBach$StanfordFreqs[c(synsites,nonsynsites)])$estimate
text(0.05,0.3,paste("R^2 = ",round(R2BachStanford,2)))

# plots all Zanini
plot(OverviewDFZanini$colMeansTsZanini,OverviewDFZanini$StanfordFreqs,
     main="Within vs between-host frequencies Zanini",col=0,
     xlim=c(0,.4),ylim=c(0,.4),
     xlab="Within-host frequencies Zanini",
     ylab="Between-host frequencies Stanford"
)
points(OverviewDFZanini$colMeansTsZanini[synsites],OverviewDFZanini$StanfordFreqs[synsites],
       col="lightgreen",pch=16)
points(OverviewDFZanini$colMeansTsZanini[nonsynsites],OverviewDFZanini$StanfordFreqs[nonsynsites],
       col="orange",pch=16)

R2ZaniniStanford<-cor.test(OverviewDFZanini$colMeansTsZanini[c(synsites,nonsynsites)],OverviewDFZanini$StanfordFreqs[c(synsites,nonsynsites)])$estimate
text(0.05,0.3,paste("R^2 = ",round(R2ZaniniStanford,2)))


# plots all Lehman
plot(OverviewDFLehman$colMeansTsLehman,OverviewDFLehman$StanfordFreqs,
     main="Within vs between-host frequencies Lehman",col=0,
     xlim=c(0,.4),ylim=c(0,.4),
     xlab="Within-host frequencies Lehman",
     ylab="Between-host frequencies Stanford"
)
points(OverviewDFLehman$colMeansTsLehman[synsites],OverviewDFLehman$StanfordFreqs[synsites],
       col="lightgreen",pch=16)
points(OverviewDFLehman$colMeansTsLehman[nonsynsites],OverviewDFLehman$StanfordFreqs[nonsynsites],
       col="orange",pch=16)

R2LehmanStanford<-cor.test(OverviewDFLehman$colMeansTsLehman[c(synsites,nonsynsites)],OverviewDFLehman$StanfordFreqs[c(synsites,nonsynsites)])$estimate
text(0.05,0.3,paste("R^2 = ",round(R2LehmanStanford,2)))


