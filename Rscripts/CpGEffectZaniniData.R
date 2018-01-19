q=10

#png("Output/EstSelCoeffZaniniDataNov27.png",width=12,height=7.5,units="in",res=100)
par(mfrow=c(1,1))
maxnuc=984
par(mar = c(3,5,1,2))
selcoeffcolumn = which(names(OverviewDFZanini)=="EstSelCoeff")
plot(OverviewDFZanini$num[40:maxnuc],OverviewDFZanini[40:maxnuc,selcoeffcolumn],
     log="y", ylab="Estimated Selection Coefficient (cost)",cex.lab=1.3,
     xaxt="n",yaxt="n", xlab="",
     col="darkgrey",t="n",pch=".", ylim=c((1/q)*3.5*10^-4,1),xlim=c(40,979))
axis(1,at=c(3*seq(15,95,by=20)-1,296+30),labels=c(seq(15,95,by=20),""))
axis(1,at=3*seq(109,349,by=20)-1,labels=seq(109-99,349-99,by=20))
axis(2,at=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),labels=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),las=1,line=0,tick=FALSE)
eaxis(side = 2, at = 10^((-0):(-(5))),label=rep("",6))

#color Protease region grey
rect(0, 0.00001, 297.5, 2, density = NULL, angle = 45,col="grey70",border = NA)
for(i in 1:5){abline(h = 1:10 * 10^(-i), col = "gray41")}

for (i in 40:maxnuc){
    c=0; co = 1
#    if (OverviewDFZanini$TypeOfSite[i]=="stop"&OverviewDFZanini$WTnt[i]%in%c("g","c")) {c=1;p=21}
#    if (OverviewDFZanini$TypeOfSite[i]=="syn"&OverviewDFZanini$WTnt[i]%in%c("g","c")) {c=cols[3];p=21; c=3}
    if (OverviewDFZanini$TypeOfSite[i]=="syn"&OverviewDFZanini$WTnt[i]%in%c("a","t")&OverviewDFBach$makesCpG[i]==0) {c=cols[6];p=21}
    if (OverviewDFZanini$TypeOfSite[i]=="syn"&OverviewDFZanini$WTnt[i]%in%c("a","t")&OverviewDFBach$makesCpG[i]==1) {c=2;p=21; print(i)}
#    if (OverviewDFZanini$TypeOfSite[i]=="nonsyn"&OverviewDFZanini$WTnt[i]%in%c("c","g")) {c=cols[4];p=21; c=3}
#    if (OverviewDFZanini$TypeOfSite[i]=="nonsyn"&OverviewDFZanini$WTnt[i]%in%c("a","t")) {c=cols[5];p=21; c=3}
    if (c!=0) points(OverviewDFZanini$num[i],OverviewDFZanini[i,selcoeffcolumn],pch=p,col=co,
                     bg=rgb(red=col2rgb(c)[1]/255,
                            green=col2rgb(c)[2]/255,
                            blue=col2rgb(c)[3]/255,
                            maxColorValue = 1,alpha=0.8),cex=2)
}
abline (h = median(OverviewDFZanini$EstSelCoeff[
    OverviewDFZanini$TypeOfSite=="syn"&OverviewDFZanini$WTnt%in%c("a","t")&OverviewDFBach$makesCpG==1]), lty=1, col=2, lwd=3)
abline (h = median(OverviewDFZanini$EstSelCoeff[
    OverviewDFZanini$TypeOfSite=="syn"&OverviewDFZanini$WTnt%in%c("a","t")&OverviewDFBach$makesCpG==0]), lty=2, col=cols[6],lwd=3)

#Add "Protease" and "RT" words
rect(0, (1/q)*0.000001, 1200, (1/q)*3.5*10^-4, density = NULL, angle = 45,col=1,border = NA)
text(55*3,(1/q)*2.9*10^-4,"PROTEASE",col="white")
rect(297.5, (1/q)*0.000001, 1200, (1/q)*3.5*10^-4, density = NULL, angle = 45,col="grey40",border = NA)
text(220*3,(1/q)*2.9*10^-4,"REVERSE TRANSCRIPTASE",col="white")

#Add legend
legpos=296;legposV=0.4
rect(legpos*3, 0.8*legposV, (legpos+42.5)*3, 1.7*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+5)*3,legposV/0.7,pch=21,bg=2,col=1,cex=2)
text((legpos+9)*3,legposV/0.7,"CpG forming",adj=0)
points((legpos+5)*3,legposV,pch=21,bg=cols[6],col=1,cex=2)
text((legpos+9)*3,legposV,"Non-CpG",adj=0)

#dev.off()

wilcox.test(
OverviewDFZanini$EstSelCoeff[
    OverviewDFZanini$TypeOfSite=="syn"&OverviewDFZanini$WTnt%in%c("a","t")&OverviewDFBach$makesCpG==1], 
OverviewDFZanini$EstSelCoeff[
    OverviewDFZanini$TypeOfSite=="syn"&OverviewDFZanini$WTnt%in%c("a","t")&OverviewDFBach$makesCpG==0], 
alternative = "greater", paired = FALSE)
#####Same for Lehman

q=10

#png("Output/EstSelCoeffLehmanDataNov27.png",width=12,height=7.5,units="in",res=100)
par(mfrow=c(1,1))
maxnuc=984
par(mar = c(3,5,1,2))
selcoeffcolumn = which(names(OverviewDFLehman)=="EstSelCoeff")
plot(OverviewDFLehman$num[40:maxnuc],OverviewDFLehman[40:maxnuc,selcoeffcolumn],
     log="y", ylab="Estimated Selection Coefficient (cost)",cex.lab=1.3,
     xaxt="n",yaxt="n", xlab="",
     col="darkgrey",t="n",pch=".", ylim=c((1/q)*3.5*10^-4,1),xlim=c(40,979))
axis(1,at=c(3*seq(15,95,by=20)-1,296+30),labels=c(seq(15,95,by=20),""))
axis(1,at=3*seq(109,349,by=20)-1,labels=seq(109-99,349-99,by=20))
axis(2,at=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),labels=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),las=1,line=0,tick=FALSE)
eaxis(side = 2, at = 10^((-0):(-(5))),label=rep("",6))

#color Protease region grey
rect(0, 0.00001, 297.5, 2, density = NULL, angle = 45,col="grey70",border = NA)
for(i in 1:5){abline(h = 1:10 * 10^(-i), col = "gray41")}

for (i in 40:maxnuc){
    c=0; co = 1
    #    if (OverviewDFLehman$TypeOfSite[i]=="stop"&OverviewDFLehman$WTnt[i]%in%c("g","c")) {c=1;p=21}
    #    if (OverviewDFLehman$TypeOfSite[i]=="syn"&OverviewDFLehman$WTnt[i]%in%c("g","c")) {c=cols[3];p=21; c=3}
    if (OverviewDFLehman$TypeOfSite[i]=="syn"&OverviewDFLehman$WTnt[i]%in%c("a","t")&OverviewDFBach$makesCpG[i]==0) {c=cols[6];p=21}
    if (OverviewDFLehman$TypeOfSite[i]=="syn"&OverviewDFLehman$WTnt[i]%in%c("a","t")&OverviewDFBach$makesCpG[i]==1) {c=2;p=21; print(i)}
    #    if (OverviewDFLehman$TypeOfSite[i]=="nonsyn"&OverviewDFLehman$WTnt[i]%in%c("c","g")) {c=cols[4];p=21; c=3}
    #    if (OverviewDFLehman$TypeOfSite[i]=="nonsyn"&OverviewDFLehman$WTnt[i]%in%c("a","t")) {c=cols[5];p=21; c=3}
    if (c!=0) points(OverviewDFLehman$num[i],OverviewDFLehman[i,selcoeffcolumn],pch=p,col=co,
                     bg=rgb(red=col2rgb(c)[1]/255,
                            green=col2rgb(c)[2]/255,
                            blue=col2rgb(c)[3]/255,
                            maxColorValue = 1,alpha=0.8),cex=2)
}
abline (h = median(OverviewDFLehman$EstSelCoeff[
    OverviewDFLehman$TypeOfSite=="syn"&OverviewDFLehman$WTnt%in%c("a","t")&OverviewDFBach$makesCpG==1]), lty=1, col=2, lwd=3)
abline (h = median(OverviewDFLehman$EstSelCoeff[
    OverviewDFLehman$TypeOfSite=="syn"&OverviewDFLehman$WTnt%in%c("a","t")&OverviewDFBach$makesCpG==0]), lty=2, col=cols[6],lwd=3)

#Add "Protease" and "RT" words
rect(0, (1/q)*0.000001, 1200, (1/q)*3.5*10^-4, density = NULL, angle = 45,col=1,border = NA)
text(55*3,(1/q)*2.9*10^-4,"PROTEASE",col="white")
rect(297.5, (1/q)*0.000001, 1200, (1/q)*3.5*10^-4, density = NULL, angle = 45,col="grey40",border = NA)
text(220*3,(1/q)*2.9*10^-4,"REVERSE TRANSCRIPTASE",col="white")

#Add legend
legpos=296;legposV=0.4
rect(legpos*3, 0.8*legposV, (legpos+42.5)*3, 1.7*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+5)*3,legposV/0.7,pch=21,bg=2,col=1,cex=2)
text((legpos+9)*3,legposV/0.7,"CpG forming",adj=0)
points((legpos+5)*3,legposV,pch=21,bg=cols[6],col=1,cex=2)
text((legpos+9)*3,legposV,"Non-CpG",adj=0)

#dev.off()

wilcox.test(
    OverviewDFLehman$EstSelCoeff[
        OverviewDFLehman$TypeOfSite=="syn"&OverviewDFLehman$WTnt%in%c("a","t")&OverviewDFBach$makesCpG==1], 
    OverviewDFLehman$EstSelCoeff[
        OverviewDFLehman$TypeOfSite=="syn"&OverviewDFLehman$WTnt%in%c("a","t")&OverviewDFBach$makesCpG==0], 
    alternative = "greater", paired = FALSE)

####Look at G to A vs A to G effect
png("../Output/EstSelCoeffZaniniDataNov27_GtoA.png",width=12,height=7.5,units="in",res=100)
par(mfrow=c(1,1))
maxnuc=984
par(mar = c(3,5,1,2))
selcoeffcolumn = which(names(OverviewDFZanini)=="EstSelCoeff")
plot(OverviewDFZanini$num[40:maxnuc],OverviewDFZanini[40:maxnuc,selcoeffcolumn],
     log="y", ylab="Estimated Selection Coefficient (cost)",cex.lab=1.3,
     xaxt="n",yaxt="n", xlab="",
     col="darkgrey",t="n",pch=".", ylim=c((1/q)*3.5*10^-4,1),xlim=c(40,979))
axis(1,at=c(3*seq(15,95,by=20)-1,296+30),labels=c(seq(15,95,by=20),""))
axis(1,at=3*seq(109,349,by=20)-1,labels=seq(109-99,349-99,by=20))
axis(2,at=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),labels=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),las=1,line=0,tick=FALSE)
eaxis(side = 2, at = 10^((-0):(-(5))),label=rep("",6))

#color Protease region grey
rect(0, 0.00001, 297.5, 2, density = NULL, angle = 45,col="grey70",border = NA)
for(i in 1:5){abline(h = 1:10 * 10^(-i), col = "gray41")}

for (i in 40:maxnuc){
    c=0; co = 1
    #    if (OverviewDFZanini$TypeOfSite[i]=="stop"&OverviewDFZanini$WTnt[i]%in%c("g","c")) {c=1;p=21}
    #    if (OverviewDFZanini$TypeOfSite[i]=="syn"&OverviewDFZanini$WTnt[i]%in%c("g","c")) {c=cols[3];p=21; c=3}
    if (OverviewDFZanini$TypeOfSite[i]=="syn"&OverviewDFZanini$WTnt[i]%in%c("g")&OverviewDFBach$makesCpG[i]==0) {c=cols[5];p=21}
    if (OverviewDFZanini$TypeOfSite[i]=="syn"&OverviewDFZanini$WTnt[i]%in%c("a")&OverviewDFBach$makesCpG[i]==0) {c=cols[6];p=21; print(i)}
    #    if (OverviewDFZanini$TypeOfSite[i]=="nonsyn"&OverviewDFZanini$WTnt[i]%in%c("c","g")) {c=cols[4];p=21; c=3}
    #    if (OverviewDFZanini$TypeOfSite[i]=="nonsyn"&OverviewDFZanini$WTnt[i]%in%c("a","t")) {c=cols[5];p=21; c=3}
    if (c!=0) points(OverviewDFZanini$num[i],OverviewDFZanini[i,selcoeffcolumn],pch=p,col=co,
                     bg=rgb(red=col2rgb(c)[1]/255,
                            green=col2rgb(c)[2]/255,
                            blue=col2rgb(c)[3]/255,
                            maxColorValue = 1,alpha=0.8),cex=2)
}
abline (h = median(OverviewDFZanini$EstSelCoeff[
    OverviewDFZanini$TypeOfSite=="syn"&OverviewDFZanini$WTnt%in%c("g")&OverviewDFBach$makesCpG==0]), lty=1, col=cols[5], lwd=3)
abline (h = median(OverviewDFZanini$EstSelCoeff[
    OverviewDFZanini$TypeOfSite=="syn"&OverviewDFZanini$WTnt%in%c("a")&OverviewDFBach$makesCpG==0]), lty=2, col=cols[6],lwd=3)

#Add "Protease" and "RT" words
rect(0, (1/q)*0.000001, 1200, (1/q)*3.5*10^-4, density = NULL, angle = 45,col=1,border = NA)
text(55*3,(1/q)*2.9*10^-4,"PROTEASE",col="white")
rect(297.5, (1/q)*0.000001, 1200, (1/q)*3.5*10^-4, density = NULL, angle = 45,col="grey40",border = NA)
text(220*3,(1/q)*2.9*10^-4,"REVERSE TRANSCRIPTASE",col="white")

#Add legend
legpos=296;legposV=0.4
rect(legpos*3, 0.8*legposV, (legpos+42.5)*3, 1.7*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+5)*3,legposV/0.7,pch=21,bg=cols[5],col=1,cex=2)
text((legpos+9)*3,legposV/0.7,"GA",adj=0)
points((legpos+5)*3,legposV,pch=21,bg=cols[6],col=1,cex=2)
text((legpos+9)*3,legposV,"AG",adj=0)

dev.off()

wilcox.test(
    OverviewDFZanini$EstSelCoeff[
        OverviewDFZanini$TypeOfSite=="syn"&OverviewDFZanini$WTnt%in%c("g")&OverviewDFBach$makesCpG==0], 
    OverviewDFZanini$EstSelCoeff[
        OverviewDFZanini$TypeOfSite=="syn"&OverviewDFZanini$WTnt%in%c("a")&OverviewDFBach$makesCpG==0], 
    alternative = "greater", paired = FALSE)


SelCoeffZanZanSynGA<-OverviewDFZanini$EstSelCoeffZan[
    OverviewDFZanini$TypeOfSite=="syn"&OverviewDFZanini$WTnt%in%c("g")&OverviewDFBach$makesCpG==0]
SelCoeffZanZanSynAG<-OverviewDFZanini$EstSelCoeffZan[
    OverviewDFZanini$TypeOfSite=="syn"&OverviewDFZanini$WTnt%in%c("a")&OverviewDFBach$makesCpG==0]

wilcox.test(
    SelCoeffZanZanSynGA, 
    SelCoeffZanZanSynAG, 
    alternative = "greater", paired = FALSE)

###ADD ZANINI MUT RATES AND DO TEST AGAIN> 

#####Same for Lehman
q=10

png("../Output/EstSelCoeffLehmanDataNov27GtoA.png",width=12,height=7.5,units="in",res=100)
par(mfrow=c(1,1))
maxnuc=984
par(mar = c(3,5,1,2))
selcoeffcolumn = which(names(OverviewDFLehman)=="EstSelCoeff")
plot(OverviewDFLehman$num[40:maxnuc],OverviewDFLehman[40:maxnuc,selcoeffcolumn],
     log="y", ylab="Estimated Selection Coefficient (cost)",cex.lab=1.3,
     xaxt="n",yaxt="n", xlab="",
     col="darkgrey",t="n",pch=".", ylim=c((1/q)*3.5*10^-4,1),xlim=c(40,979))
axis(1,at=c(3*seq(15,95,by=20)-1,296+30),labels=c(seq(15,95,by=20),""))
axis(1,at=3*seq(109,349,by=20)-1,labels=seq(109-99,349-99,by=20))
axis(2,at=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),labels=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),las=1,line=0,tick=FALSE)
eaxis(side = 2, at = 10^((-0):(-(5))),label=rep("",6))

#color Protease region grey
rect(0, 0.00001, 297.5, 2, density = NULL, angle = 45,col="grey70",border = NA)
for(i in 1:5){abline(h = 1:10 * 10^(-i), col = "gray41")}

for (i in 40:maxnuc){
    c=0; co = 1
    #    if (OverviewDFLehman$TypeOfSite[i]=="stop"&OverviewDFLehman$WTnt[i]%in%c("g","c")) {c=1;p=21}
    #    if (OverviewDFLehman$TypeOfSite[i]=="syn"&OverviewDFLehman$WTnt[i]%in%c("g","c")) {c=cols[3];p=21; c=3}
    if (OverviewDFLehman$TypeOfSite[i]=="syn"&OverviewDFLehman$WTnt[i]%in%c("g")&OverviewDFBach$makesCpG[i]==0) {c=cols[5];p=21}
    if (OverviewDFLehman$TypeOfSite[i]=="syn"&OverviewDFLehman$WTnt[i]%in%c("a")&OverviewDFBach$makesCpG[i]==0) {c=3;p=21; print(i)}
    #    if (OverviewDFLehman$TypeOfSite[i]=="nonsyn"&OverviewDFLehman$WTnt[i]%in%c("c","g")) {c=cols[4];p=21; c=3}
    #    if (OverviewDFLehman$TypeOfSite[i]=="nonsyn"&OverviewDFLehman$WTnt[i]%in%c("a","t")) {c=cols[5];p=21; c=3}
    if (c!=0) points(OverviewDFLehman$num[i],OverviewDFLehman[i,selcoeffcolumn],pch=p,col=co,
                     bg=rgb(red=col2rgb(c)[1]/255,
                            green=col2rgb(c)[2]/255,
                            blue=col2rgb(c)[3]/255,
                            maxColorValue = 1,alpha=0.8),cex=2)
}
abline (h = median(OverviewDFLehman$EstSelCoeff[
    OverviewDFLehman$TypeOfSite=="syn"&OverviewDFLehman$WTnt%in%c("g")&OverviewDFBach$makesCpG==0]), lty=1, col=2, lwd=3)
abline (h = median(OverviewDFLehman$EstSelCoeff[
    OverviewDFLehman$TypeOfSite=="syn"&OverviewDFLehman$WTnt%in%c("a")&OverviewDFBach$makesCpG==0]), lty=2, col=cols[6],lwd=3)

#Add "Protease" and "RT" words
rect(0, (1/q)*0.000001, 1200, (1/q)*3.5*10^-4, density = NULL, angle = 45,col=1,border = NA)
text(55*3,(1/q)*2.9*10^-4,"PROTEASE",col="white")
rect(297.5, (1/q)*0.000001, 1200, (1/q)*3.5*10^-4, density = NULL, angle = 45,col="grey40",border = NA)
text(220*3,(1/q)*2.9*10^-4,"REVERSE TRANSCRIPTASE",col="white")

#Add legend
legpos=296;legposV=0.4
rect(legpos*3, 0.4*legposV, (legpos+42.5)*3, 1.7*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+5)*3,legposV/0.7,pch=21,bg=1,col=1,cex=2)
text((legpos+9)*3,legposV/0.7,"Nonsense",adj=0)
points((legpos+5)*3,legposV,pch=21,bg=cols[4],col=1,cex=2)
text((legpos+9)*3,legposV,"Non-syn, C/G",adj=0)
points((legpos+5)*3,legposV*0.7,pch=21,bg=cols[5],col=1,cex=2)
text((legpos+9)*3,legposV*0.7,"Non-syn, A/T",adj=0)
points((legpos+5)*3,legposV*0.49,pch=21,bg=cols[3],col=1,cex=2)
text((legpos+9)*3,legposV*0.49,"Synonymous",adj=0)

dev.off()

wilcox.test(
    OverviewDFLehman$EstSelCoeff[
        OverviewDFLehman$TypeOfSite=="syn"&OverviewDFLehman$WTnt%in%c("g")&OverviewDFBach$makesCpG==0], 
    OverviewDFLehman$EstSelCoeff[
        OverviewDFLehman$TypeOfSite=="syn"&OverviewDFLehman$WTnt%in%c("a")&OverviewDFBach$makesCpG==0], 
    alternative = "greater", paired = FALSE)


#Look at median for Bacheler data
median(OverviewDFBach$EstSelCoeff[
    OverviewDFBach$TypeOfSite=="syn"&OverviewDFBach$WTnt%in%c("a","t")&OverviewDFBach$makesCpG==1])

median(OverviewDFBach$EstSelCoeff[
    OverviewDFBach$TypeOfSite=="syn"&OverviewDFBach$WTnt%in%c("a","t")&OverviewDFBach$makesCpG==0])
