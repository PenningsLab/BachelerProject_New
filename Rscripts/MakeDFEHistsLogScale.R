setwd("/Users/pleuni/Documents/Git/bachelerProject")
source("Rscripts/baseRscript.R")

OverviewDFBacheler <- read.table("Output/OverviewSelCoeff_BachelerFilter.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFLehman <- read.table("Output/OverviewSelCoeffLehman.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFZanini <- read.table("Output/OverviewSelCoeffZanini.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

breaks=seq(-5,0,by=.2)

for(dat.file in c("Lehman", "Zanini", "Bacheler")){
        filename<-paste("Output/DFE_log_May2018", dat.file,".pdf", sep = "_")
        pdf(filename, width = 12, height = 7)

        par(mfcol=c(2,4))
        par(mar = c(4,4.3,2.5,1))
        for (typeofsite in c("syn", "nonsyn")){
        
            if (dat.file == "Lehman") dat = OverviewDFLehman
            if (dat.file == "Bacheler") dat = OverviewDFBacheler
            if (dat.file == "Zanini") dat = OverviewDFZanini
        
        i=1
        mycol <- rgb(190, 190, 190, max = 255, alpha = 100, names = "grey50")
        for (wtnt in c("a", "t", "c", "g")){
            
                datavector<- dat$EstSelCoeff[OverviewDFBacheler$TypeOfSite==typeofsite & OverviewDFBacheler$WTnt==wtnt]
            hist(log10(datavector), breaks = breaks, xlim=c(-4,-0),col=cols[i], 
                 main = paste0(wtnt, " : ", typeofsite), ylab="Count", xlab="",
                 xaxt="n", 
                 cex.main=1.7, las=1, cex.axis =1.2,cex.lab=1.6
                 )
            datavectorNonCpG<- dat$EstSelCoeff[OverviewDFBacheler$TypeOfSite==typeofsite & OverviewDFBacheler$WTnt==wtnt & OverviewDFBacheler$makesCpG==0]
            if (i<0) hist(log10(datavectorNonCpG), breaks = breaks, xlim=c(-5,-0),col=mycol, 
                 main = paste0(wtnt, " : ", typeofsite), ylab="Count", xlab="",xaxt="n",add=TRUE)
            #This last line highlights the effect of CpG creating mutations, but I didn't use it. 
            abline(v=median(log10(datavector), na.rm = TRUE),col="white",lwd=3,lty=1)
            abline(v=median(log10(datavector), na.rm = TRUE),lwd=3,lty=2)
            
            x=seq(-4,0,by=1)
            axis(1, at=x,labels=paste("10^",x,sep=""), col.axis="black", las=1, line=-.5, cex.axis=1.2)
            mtext(text = "Selection Coefficient", side = 1, line = 1.8, cex=1.1)
            
            i=i+1
        }
    }
        dev.off()
}

