

setwd("~/Documents/Git/bachelerProject")

#Mut rates and sel coefficients
read.csv("Data/HIVMutRates/HIVMutRates.csv")->mutrates

included<-c(which(mutrates$Nucleotide.substitution=="AG"),
            which(mutrates$Nucleotide.substitution=="UC"),
            which(mutrates$Nucleotide.substitution=="CU"),
            which(mutrates$Nucleotide.substitution=="GA"))
pdf("Output/MutationRatesUsed.pdf")
par(mar = c(4.5, 6, 4, 0.5))
barplot(t(as.matrix(mutrates[included,2:3])),beside=TRUE,ylim=c(0,6*10^-5),
        xaxt="n",yaxt="n",ylab="",main="")
xvalues=c(2,5,8,11)
axis(1,at=xvalues,labels=c("A-G","T-C","C-T","G-A"),tick = FALSE)
box()
axis(2,at=seq(0,6*10^-5,by=10^-5),labels=seq(0,6*10^-5,by=10^-5),tick = TRUE,las =2)
mtext("Estimated mutation rates", side=2, line=4,cex=1.3)
mtext("Mutation rate estimates used", side=3, line=1,cex=1.4)
text(xvalues-.5,mutrates[included,2]+0.1*10^-5,"Abram")
text(xvalues+.5,mutrates[included,3]+0.1*10^-5,"Zanini")
dev.off()
