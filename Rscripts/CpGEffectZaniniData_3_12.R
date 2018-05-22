#setwd("/Users/pleuni/Documents/Git/bachelerProject")
source("Rscripts/baseRscript.R")

OverviewDFBacheler <- read.table("Output/OverviewSelCoeff_BachelerFilter.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFLehman <- read.table("Output/OverviewSelCoeffLehman.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFZanini <- read.table("Output/OverviewSelCoeffZanini.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
OverviewDFZanini$ makesCpG <- OverviewDFBacheler$makesCpG
OverviewDFLehman$ makesCpG <- OverviewDFBacheler$makesCpG
OverviewDFZanini$bigAAChange<-OverviewDFBacheler$bigAAChange

#Group B is the one that we expect to be more costly
dat.file="Zanini"
groupARows<-which(dat$TypeOfSite=="syn"&dat$WTnt%in%c("a")&dat$makesCpG==0)
groupBRows<-which(dat$TypeOfSite=="syn"&dat$WTnt%in%c("g")&dat$makesCpG==0)
nameA = "noCpG, syn, a"
nameB = "noCpG, syn, g"

MakePlot<-function(dat.file, mutrate = "Abrams", groupARows, groupBRows, nameA, nameB){
  #if (dat.file == "Lehman") dat = OverviewDFLehman
  #if (dat.file == "Bacheler") dat = OverviewDFBacheler
  #if (dat.file == "Zanini") dat = OverviewDFZanini
  if (mutrate =="Abrams") selcoeffcolumn = which(names(dat)=="EstSelCoeff")
  if (mutrate =="Zanini") selcoeffcolumn = which(names(dat)=="EstSelCoeffZan")
#  par(mfrow=c(1,1))
  maxnuc=984; q=10;
  par(mar = c(3,5,1,2))
  plot(dat$num[40:maxnuc],dat[40:maxnuc,selcoeffcolumn],
       log="y", ylab="Estimated Selection Coefficient (cost)",cex.lab=1.3,
       xaxt="n",yaxt="n", xlab="", 
       main = paste("Frequencies=",dat.file, ", Mutrate=",mutrate),
       col="darkgrey",t="n",pch=".", ylim=c((1/q)*3.5*10^-4,1),xlim=c(40,979))
  axis(1,at=c(3*seq(15,95,by=20)-1,296+30),labels=c(seq(15,95,by=20),""))
  axis(1,at=3*seq(109,349,by=20)-1,labels=seq(109-99,349-99,by=20))
  axis(2,at=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),labels=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),las=1,line=0,tick=FALSE)
  eaxis(side = 2, at = 10^((-0):(-(5))),label=rep("",6))
  #color Protease region grey
  rect(0, 0.00001, 297.5, 2, density = NULL, angle = 45,col="grey70",border = NA)
  for(i in 1:5){abline(h = 1:10 * 10^(-i), col = "gray41")}
  
  c = cols[6]; colA = rgb(red=col2rgb(c)[1]/255,green=col2rgb(c)[2]/255,blue=col2rgb(c)[3]/255,maxColorValue = 1,alpha=0.8)
  c = 2; colB = rgb(red=col2rgb(c)[1]/255,green=col2rgb(c)[2]/255,blue=col2rgb(c)[3]/255,maxColorValue = 1,alpha=0.8)
  p=21; co=1
  for (i in groupARows){
    points(dat$num[i],dat[i,selcoeffcolumn],pch=p,col=co,
           bg=colA,cex=2)}
  c = 2
  for (i in groupBRows){
    points(dat$num[i],dat[i,selcoeffcolumn],pch=p,col=co,
           bg=colB,cex=2)}
  
  abline (h = median(dat[groupARows,selcoeffcolumn]), lty=1, col=colA, lwd=3)
  abline (h = median(dat[groupBRows,selcoeffcolumn]), lty=1, col=colB, lwd=3)
  
  #Add "Protease" and "RT" words
  rect(0, (1/q)*0.000001, 1200, (1/q)*3.5*10^-4, density = NULL, angle = 45,col=1,border = NA)
  text(55*3,(1/q)*2.9*10^-4,"PROTEASE",col="white")
  rect(297.5, (1/q)*0.000001, 1200, (1/q)*3.5*10^-4, density = NULL, angle = 45,col="grey40",border = NA)
  text(220*3,(1/q)*2.9*10^-4,"REVERSE TRANSCRIPTASE",col="white")
  
  #Add legend
  legpos=246;legposV=0.4
  rect(legpos*3, 0.8*legposV, (legpos+80)*3, 1.9*legposV, density = NULL, angle = 45,col=alpha("white",1))
  points((legpos+5)*3,legposV/0.7,pch=21,bg=colA,col=1,cex=2)
  text((legpos+9)*3,legposV/0.7,nameA,adj=0)
  points((legpos+5)*3,legposV,pch=21,bg=colB,col=1,cex=2)
  text((legpos+9)*3,legposV,nameB,adj=0)
  
  #add p-value to the plot
  pvalue<-wilcox.test(dat[groupBRows,selcoeffcolumn], dat[groupARows,selcoeffcolumn],
    alternative = "greater", paired = FALSE)$p.value
  #text(200,0.15,paste("pvalue =", pvalue))
  text(200,0.15,paste("pvalue =", format(round(pvalue, 4), nsmall = 4)))
}

#Check CpG effect 
#pdf ("CpG_Effect.pdf")
par(mfrow=c(2,2))
#for (mutrate in c("Abrams","Zanini")){
  for (dat.file in c("Bacheler","Zanini")){
    if (dat.file == "Bacheler") dat = OverviewDFBacheler
    if (dat.file == "Zanini") dat = OverviewDFZanini
    groupARows<-which(dat$TypeOfSite=="syn"&dat$WTnt%in%c("t")&dat$makesCpG==0)
    groupBRows<-which(dat$TypeOfSite=="syn"&dat$WTnt%in%c("t")&dat$makesCpG==1)
    nameA = "noCpG, syn, a/t"
    nameB = "CpG, syn, a/t"
    MakePlot(dat.file, mutrate,groupARows,groupBRows, nameA, nameB)
    text(300,0.0001,"t only")
  }
#}
#dev.off()

#Check GA vs AG effect in syn sites
#pdf ("GA_vs_AG_Effect.pdf")
for (mutrate in c("Abrams","Zanini")){
  for (dat.file in c("Bacheler","Zanini")){
    if (dat.file == "Bacheler") dat = OverviewDFBacheler
    if (dat.file == "Zanini") dat = OverviewDFZanini
    groupARows<-which(dat$TypeOfSite=="syn"&dat$WTnt%in%c("a")&dat$makesCpG==0)
    groupBRows<-which(dat$TypeOfSite=="syn"&dat$WTnt%in%c("g")&dat$makesCpG==0)
    nameA = "noCpG, syn, a"
    nameB = "noCpG, syn, g"
    MakePlot(dat.file, mutrate,groupARows,groupBRows, nameA, nameB)
  }
}
#dev.off()

#Check CT vs AG effect in syn sites
#pdf ("CT_vs_AG_Effect.pdf")
for (mutrate in c("Abrams","Zanini")){
  for (dat.file in c("Bacheler","Zanini")){
    if (dat.file == "Bacheler") dat = OverviewDFBacheler
    if (dat.file == "Zanini") dat = OverviewDFZanini
    groupARows<-which(dat$TypeOfSite=="syn"&dat$WTnt%in%c("a")&dat$makesCpG==0)
    groupBRows<-which(dat$TypeOfSite=="syn"&dat$WTnt%in%c("c")&dat$makesCpG==0)
    nameA = "noCpG, syn, a"
    nameB = "noCpG, syn, c"
    MakePlot(dat.file, mutrate,groupARows,groupBRows, nameA, nameB)
  }
}
#dev.off()
  
#Check GA vs AG effect in nonsyn sites
#pdf ("GA_vs_AG__nonsyn_nondrastic_Effect.pdf")
for (mutrate in c("Abrams","Zanini")){
  for (dat.file in c("Bacheler","Zanini")){
    if (dat.file == "Bacheler") dat = OverviewDFBacheler
    if (dat.file == "Zanini") dat = OverviewDFZanini
    groupARows<-which(dat$TypeOfSite=="nonsyn"&dat$WTnt%in%c("a")&dat$makesCpG==0&dat$bigAAChange==0)
    groupBRows<-which(dat$TypeOfSite=="nonsyn"&dat$WTnt%in%c("g")&dat$makesCpG==0&dat$bigAAChange==0)
    nameA = "noCpG, nonsyn, a"
    nameB = "noCpG, nonsyn, g"
    MakePlot(dat.file, mutrate,groupARows,groupBRows, nameA, nameB)
    if (mutrate =="Abrams") selcoeffcolumn = which(names(dat)=="EstSelCoeff")
    if (mutrate =="Zanini") selcoeffcolumn = which(names(dat)=="EstSelCoeffZan")
    print(paste(dat.file, mutrate, round((median(dat[groupBRows,selcoeffcolumn])/median(dat[groupARows,selcoeffcolumn])-1),3)))
  }
}
#dev.off()


#Check GA vs AG effect in nonsyn sites
pdf ("CT_vs_AG__nonsyn_nondrastic_Effect.pdf")
for (mutrate in c("Abrams","Zanini")){
  for (dat.file in c("Bacheler","Zanini")){
    if (dat.file == "Bacheler") dat = OverviewDFBacheler
    if (dat.file == "Zanini") dat = OverviewDFZanini
    groupARows<-which(dat$TypeOfSite=="nonsyn"&dat$WTnt%in%c("a")&dat$makesCpG==0&dat$bigAAChange==0)
    groupBRows<-which(dat$TypeOfSite=="nonsyn"&dat$WTnt%in%c("c")&dat$makesCpG==0&dat$bigAAChange==0)
    nameA = "noCpG, nonsyn, a"
    nameB = "noCpG, nonsyn, c"
    MakePlot(dat.file, mutrate,groupARows,groupBRows, nameA, nameB)
    if (mutrate =="Abrams") selcoeffcolumn = which(names(dat)=="EstSelCoeff")
    if (mutrate =="Zanini") selcoeffcolumn = which(names(dat)=="EstSelCoeffZan")
    print(paste(dat.file, mutrate, round((median(dat[groupBRows,selcoeffcolumn])/median(dat[groupARows,selcoeffcolumn])-1),3)))
  }
}
dev.off()

#Check GA vs AG effect in nonsyn sites
#filter low frequencies
#pdf ("GA_vs_AG__nonsyn_Effect.pdf")
for (mutrate in c("Abrams","Zanini")){
  for (dat.file in c("Bacheler","Zanini")){
    if (dat.file == "Bacheler") dat = OverviewDFBacheler
    if (dat.file == "Zanini") {dat = OverviewDFZanini
    groupARows<-which(dat$TypeOfSite=="nonsyn"&dat$WTnt%in%c("a")&dat$makesCpG==0)#&dat$bigAAChange==0)
    groupBRows<-which(dat$TypeOfSite=="nonsyn"&dat$WTnt%in%c("g")&dat$makesCpG==0)#&dat$bigAAChange==0)
    nameA = "noCpG, nonsyn, a"
    nameB = "noCpG, nonsyn, g"
    MakePlot(dat.file, mutrate,groupARows,groupBRows, nameA, nameB)
  }
}
dev.off()