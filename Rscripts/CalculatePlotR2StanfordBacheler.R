## Step 4:  R analysis for plotting frequencies: correlation with bacheler dataset 

#Get Bacheler frequencies previously calculated ** 
library(scales)

#setwd("~/Documents/Git/bachelerProject/Comparison-stanford_vs_bacheler-april2016")
#setwd("~/Dropbox/MarionKristofBachelerProject/GitMarionKristof/bachelerProject/Comparison-stanford_vs_bacheler-april2016")
source("../Rscripts/baseRscript.R")
OverviewDFBach<-read.csv("Output/OverviewSelCoeff_BachelerFilter.csv",row.names=1)
OverviewDFBach$StanfordFreqs<-t(read.csv("Output/subtypeB-frequencies_stanford_pr_rt_naive.csv"))[1:length(OverviewDFBach[,1])]

OverviewDFZanini<-read.csv("../Output/OverviewSelCoeffZanini.csv",row.names=1)
OverviewDFZanini$StanfordFreqs<-t(read.csv("subtypeB-frequencies_stanford_pr_rt_naive.csv"))[1:length(OverviewDFBach[,1])]

OverviewDFLehman<-read.csv("../Output/OverviewSelCoeffLehman.csv",row.names=1)
OverviewDFLehman$StanfordFreqs<-t(read.csv("subtypeB-frequencies_stanford_pr_rt_naive.csv"))[1:length(OverviewDFBach[,1])]

png("Output/StanfordVsBacheler2017Nov27.png",width=8,height=8,units="in",res=100)

#Plot Bacheler
par(mfrow=c(1,1))
par(mar=c(5.1,6.1,4.1,2.1))
plot(OverviewDFBach$MeanFreq+0.0001,OverviewDFBach$StanfordFreqs+0.0001,
     main="Within vs between-host frequencies Bacheler",col=0,
     xlim=c(0.0001,.05),ylim=c(0.0001,.2),
     xlab="Within-host frequencies Bacheler",
#     ylab="",
ylab="",
    log="xy",
     xaxt="n", yaxt="n"
)


x <- c(10^-4,5*10^-4,10^-3, 5*10^-3, 10^-2)
axis(side=1,
     at = c(10^-4,x+0.0001), 
        labels= c(0,scientific_format(1)(x)))
x <- c(5*10^-5,10^-4,5*10^-4,10^-3, 5*10^-3, 10^-2)
axis(side=2,
     at = c(10^-4,x+0.0001), 
     labels= c(0,scientific_format(1)(x)), las=2)

axis.break(2, 1.2*10^-4, style = "slash")
axis.break(1, 1.3*10^-4, style = "slash")

mtext("Between-host frequencies Stanford", side=2, line=4, cex.lab=1,las=0, col="black")
points(OverviewDFBach$MeanFreq[SynSites]+0.0001,OverviewDFBach$StanfordFreqs[SynSites]+0.0001,
       col=cols[3],pch=16)
points(OverviewDFBach$MeanFreq[NonSynSites]+0.0001,OverviewDFBach$StanfordFreqs[NonSynSites]+0.0001,
       col=cols[5],pch=16)

R2BachStanford<-cor.test(
    OverviewDFBach$MeanFreq[c(SynSites,NonSynSites)],
    OverviewDFBach$StanfordFreqs[c(SynSites,NonSynSites)],
    method="s")$estimate
text(0.0002,0.1,paste("R^2 = ",round(R2BachStanford,2)))

dev.off()


