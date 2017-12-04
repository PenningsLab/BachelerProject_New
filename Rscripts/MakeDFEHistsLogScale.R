setwd("/Users/pleuni/Documents/Git/bachelerProject")
source("Rscripts/baseRscript.R")

OverviewDFBacheler <- read.table("Output/OverviewSelCoeff_BachelerFilter.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFLehman <- read.table("Output/OverviewSelCoeffLehman.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFZanini <- read.table("Output/OverviewSelCoeffZanini.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

breaks=seq(-5,0,by=.2)

for(dat.file in c("Lehman", "Zanini", "Bacheler")){
    for (typeofsite in c("syn", "nonsyn")){
        filename<-paste("Output/DFE_log_Nov2017", dat.file, typeofsite,".pdf", sep = "_")
        pdf(filename)
        #if (typeofsite =="nonsyn") pdf("Output/DFE_log_nonsyn.pdf")
        if (dat.file == "Lehman") dat = OverviewDFLehman
        if (dat.file == "Bacheler") dat = OverviewDFBacheler
        if (dat.file == "Zanini") dat = OverviewDFZanini
        
        par(mfrow=c(2,2))
        i=1
        for (wtnt in c("a", "t", "c", "g")){
            datavector<- dat$EstSelCoeff[OverviewDFBacheler$TypeOfSite==typeofsite & OverviewDFBacheler$WTnt==wtnt]
            hist(log10(datavector), breaks = breaks, xlim=c(-5,-0),col=cols[i], 
                 main = paste0(wtnt, " : ", typeofsite), ylab="Count", xlab="",xaxt="n")
            
            abline(v=median(log10(datavector), na.rm = TRUE),col="white",lwd=3,lty=1)
            abline(v=median(log10(datavector), na.rm = TRUE),lwd=3,lty=2)
            
            x=seq(-4,0,by=1)
            axis(1, at=x,labels=paste("10^",x,sep=""), col.axis="black", las=1, line=-.2)
            mtext(text = "Selection Coefficient", side = 1, line = 2, cex=.9)
            
            i=i+1
        }
        dev.off()
    }
}
