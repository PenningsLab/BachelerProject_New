source("GLM_Analysis/code/prepareDataForGLM.R")
source("GLM_Analysis/code/helperFunctionsForF2.R")

#Substantially more complicated model with structural elements
#fullmodel.int <- glm(minor ~ t + c + g + bigAAChange + inRT + t*nonsyn + c*nonsyn + g*nonsyn + shape + CpG + CpG*t  + CpG*nonsyn + CpG*nonsyn*t + helix*nonsyn + beta*nonsyn + coil*nonsyn,  family = "binomial", data = datFitModel[datFitModel$res == 0 & datFitModel$stop == 0,])

#Run model 
fullmodel.int <- glm(minor ~ inRT + shape + t + c + g + CpG + CpG*t  + t*nonsyn + c*nonsyn + g*nonsyn + nonsyn*CpG + t:nonsyn:CpG + bigAAChange,  family = "binomial", data = datFitModel[datFitModel$res == 0 & datFitModel$stop == 0,])
sumOfModel <- summary(fullmodel.int)

#Create a table for latex 
require(xtable)
#This doesn't work yet: #write(xtable(sumOfModel, digits = 3), "ModelTable.txt")
xtable(sumOfModel, digits = 3)
print(xtable(sumOfModel, digits = 3),type="html",file="GLM_Analysis/out/SumOfModel1.html")

#fullmodel.small <- glm(minor ~ t + c + g + bigAAChange + t*nonsyn + c*nonsyn + g*nonsyn + CpG + CpG*t  + CpG*nonsyn + CpG*nonsyn*t,  family = "binomial", data = datFitModel[datFitModel$res == 0 & datFitModel$stop == 0,])
#sumOfModel <- summary(fullmodel.small)

modcoef <- sumOfModel$coef

coef.vals <- modcoef[,1]
coef.pSE.vals <- coef.vals + modcoef[,2]
coef.mSE.vals <- coef.vals - modcoef[,2]

inPRrows <- intersect(which(suppDat$TypeOfSite != "overlap"), 1:297)
inRTrows <- intersect(which(suppDat$TypeOfSite != "overlap"), 298:984)

#Pleuni: not sure what this does
tmpDat.0 <- DataFrameOfData()
tmpDat.0[,8] <- 0 #shape
tmpDat.1 <- DataFrameOfData()
tmpDat.1[,8] <- 1 #shape

mean(exp(tmpDat.0 %*% coef.vals))
mean(exp(tmpDat.1 %*% coef.vals))

#syn vs nonsyn
nonsynrows <- which(suppDat$TypeOfSite == "nonsyn")
synrows <- which(suppDat$TypeOfSite == "syn")
mean(exp(DataFrameOfData()[nonsynrows,] %*% coef.vals))
mean(exp(DataFrameOfData()[synrows,] %*% coef.vals))
mean(suppDat$colMeansTs0[nonsynrows])
mean(suppDat$colMeansTs0[synrows])


#Let's do some fresh analysis on this for Marion.
#What's the average shape coefficient?
avShape <- mean(suppDatShape$SHAPE)
avShape
inRTval <- 1
mus <- c(1.11e-05, 1.11e-05, 2.41e-05, 5.48e-05)

# Point 1:
#synonymous CpG forming mutations
CpGSyn<-mus/exp(makeDataFrameToModify.withSHAPEandinRT(0,1,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

#synonymous non-CpG forming mutations
NonCpGSyn<-mus/exp(makeDataFrameToModify.withSHAPEandinRT(0,0,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

#Magnitude changes
magchanges<-(mus/exp(makeDataFrameToModify.withSHAPEandinRT(0,1,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals))/(mus/exp(makeDataFrameToModify.withSHAPEandinRT(0,0,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals))

cat(        "\n\nPOINT 1\n" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=FALSE,sep="\n")


cat("Using model-predicted frequencies and known mutation rates, we find that CpG-creating synonymous mutations are ",file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE)
cat(round(mean(magchanges[1:2])),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE)
cat(" times more costly (selection coefficient appr. ",file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE)
cat(round(mean(CpGSyn[1:2]),3),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE)
cat( ") than non-CpG-creating synonymous mutations (selection coefficient ~",file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE)
cat(round(mean(NonCpGSyn[1:2]),5),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE)
cat(" )).\n",file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE)


cat("\nMore detail \nA-G mutations\n",file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE)
cat(round(mean(magchanges[1]),2),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(mean(CpGSyn[1]),4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(mean(NonCpGSyn[1]),6),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")

cat("\nC-T mutations\n",file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE)
cat(round(mean(magchanges[2]),2),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(mean(CpGSyn[2]),4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(mean(NonCpGSyn[2]),6),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")

# Point 2: Pleuni: not sure which point this is. 
AGGA<-mus/exp(makeDataFrameToModify.withSHAPEandinRT(0,0,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

cat(        "\n\nPOINT 2\n" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")

cat("\n\nIndeed, the estimated selection coefficients based on model predictions suggested that synonymous G to A mutations are" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(AGGA[4]/AGGA[1],2),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat("times as costly as non-CpG-forming A to G mutations (",file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(AGGA[4],4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat("vs",file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(AGGA[1],4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(").",file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")

cat("\n Note that p-values come from analyseAndFigures-Bacheler.Rmd",file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")

#  Point 3:
#Non-synonymous (changes AA group versus doesn’t change AA group) - among non-CpG forming mutations. 

#Does not change AA group:
NotDrastic<-mus/exp(makeDataFrameToModify.withSHAPEandinRT(1,0,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

#Does change AA group:
Drastic<-mus/exp(makeDataFrameToModify.withSHAPEandinRT(1,0,1, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

#Magnitude change:
MagChange<-(mus/exp(makeDataFrameToModify.withSHAPEandinRT(1,0,1, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals))/(mus/exp(makeDataFrameToModify.withSHAPEandinRT(1,0,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals))
cat(        "\n\nPOINT 3\n" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")

cat("\n\nIn general, mutations that led to a drastic amino acid change were found at lower frequency than mutations that did not ($p < 0.001$).",
    file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")

cat("\nFor example, A to G mutations that result in a drastic amino acid change are roughly" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(MagChange[1],2),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat("times more costly than A to G mutations that do not (",file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(Drastic[1],4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
    cat(        "vs" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
    cat(round(NotDrastic[1],4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
    cat(       ").",file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n") 
cat("We observed similar fold changes for the other possible transitions.",file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")

#Pleuni: continue here
#  Point 4:
#Nonsynonymous, does not change AA group, does not create new CpG:
noNewCpG <- mus/exp(makeDataFrameToModify.withSHAPEandinRT(1,0,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

#Nonsynonymous, does not change AA group, does create new CpG:
NewCpG <- mus/exp(makeDataFrameToModify.withSHAPEandinRT(1,1,0, avShape, inRTval)[,rownames(modcoef)] %*% coef.vals)

#Magnitude change
NewCpG/noNewCpG

cat(        "\n\nPOINT 4\n" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")

cat("\nThere was also an effect of whether or not a non-synonymous mutation created a 
CpG site ($p < 0.001$ for both A-G and T-C mutations). 
The difference in frequencies suggests that, 
among mutations that do not lead to a drastic amino acid change, 
    A-G mutations that create a CpG site are approximately" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(NewCpG[1]/noNewCpG[1],2),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(        "times more costly than those that do not " ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(NewCpG[1],4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(        "vs" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(noNewCpG[1],4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")

cat(        "\n\nSimilarly, 
among mutations that do not lead to a drastic amino acid change, 
            T-C mutations that create a CpG site are approximately " ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(NewCpG[2]/noNewCpG[2],2),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat( "times more costly than those that do not" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(NewCpG[2],4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
    cat(        "vs" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")

    cat(round(noNewCpG[2],4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")


#  Point 5
#Non synonymous, non-CpG forming, does not change AA group:
noNewCpG

cat(        "\n\nPOINT 5\n" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")

cat(        "\nWe estimated that, among non-synonymous mutations that do not involve a drastic amino acid change 
or create a CpG site, 
            C-T mutations are \n" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")

cat(round(noNewCpG[3]/noNewCpG[1],2),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(        "times more costly than A-G mutations " ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(noNewCpG[3],4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
    cat(        "vs" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
    cat(round(noNewCpG[1],4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
        cat(        ", and G-A mutations are " ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(noNewCpG[4]/noNewCpG[2],2),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(        "times more costly than A-G mutations" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
cat(round(noNewCpG[4],4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
    cat(        "vs" ,file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")
    cat(round(noNewCpG[2],4),file = "GLM_Analysis/out/GLMResultsText.txt", append=TRUE,sep="\n")


#Pleuni continue here. 
    
    #geeeeee
library(plotrix)
require(RColorBrewer)
pdf("GLM_Analysis/graphs/modeled_freqs_May2017_2.pdf", width = 12, height = 7)
cols <- brewer.pal(4, "Set2")
layout(matrix(1:2, nrow = 1))
par(mar = c(4, 4.5, 1.5, 1))
makePlot.axisbreak(main = "Synonymous Sites")
plotVals(0, 1, 0, cols[1], .1)
plotVals(0, 0, 0, cols[2], -.1)
plotDat(0, 1, 0, cols[1], .1)
plotDat(0, 0, 0, cols[2], -.1)
abline(v = 1:3 + .5, col = "black")
#legend("topleft", c("CpG-forming", "non-CpG-forming"), col = cols[1:2], pch = 16, bg = "white")
legend("bottomright", c("No drastic AA change (non-CpG-forming)", "No drastic AA change (CpG-forming)", "Drastic AA change (non-CpG-forming)",  "Drastic AA change (CpG-forming)"), col = cols[c(2,1,3,4)], pch = 16, bg = "white" )
makePlot.axisbreak(main = "Non-synonymous Sites")
plotDat(1, 0, 1, cols[3], .1)
plotDat(1, 0, 0, cols[2], -.3)
plotDat(1, 1, 0, cols[1], -.1)
plotDat(1, 1, 1, cols[4], .3)
plotVals(1, 0, 1, cols[3], .1)
plotVals(1, 0, 0, cols[2], -.3)
plotVals(1, 1, 0, cols[1], -.1 )
plotVals(1, 1, 1, cols[4], .3 )
abline(v = 1:3 + .5, col = "black")
dev.off()

pdf("GLM_Analysis/graphs/modeled_sels_May2017.pdf", width = 12, height = 7)
cols <- brewer.pal(4, "Set2")
layout(matrix(1:2, nrow = 1))
par(mar = c(4, 4.5, 1.5, 1))
makePlot.svals(main = "Synonymous Sites")
plotVals.svals(0, 1, 0, cols[1], .1)
plotVals.svals(0, 0, 0, cols[2], -.1)
plotDat.svals(0, 1, 0, cols[1], .1)
plotDat.svals(0, 0, 0, cols[2], -.1)
abline(v = 1:3 + .5, col = "black")
                                        #legend("topleft", c("CpG-forming", "non-CpG-forming"), col = cols[1:2], pch = 16, bg = "white")
legend("topleft", c("No drastic AA change (non-CpG-forming)", "No drastic AA change (CpG-forming)", "Drastic AA change (non-CpG-forming)",  "Drastic AA change (CpG-forming)"), col = cols[c(2,1,3,4)], pch = 16, bg = "white" )
#legend("topleft", c("Same AA group (non-CpG-forming)", "Same AA group (CpG-forming)", "Changes AA group (non-CpG-forming)",  "Changes AA group (CpG-forming)"), col = cols[c(2,1,3,4)], pch = 16, bg = "white" )
makePlot.svals(main = "Non-synonymous Sites")
plotDat.svals(1, 0, 1, cols[3], .1)
plotDat.svals(1, 0, 0, cols[2], -.3)
plotDat.svals(1, 1, 0, cols[1], -.1)
plotDat.svals(1, 1, 1, cols[4], .3)
plotVals.svals(1, 0, 1, cols[3], .1)
plotVals.svals(1, 0, 0, cols[2], -.3)
plotVals.svals(1, 1, 0, cols[1], -.1 )
plotVals.svals(1, 1, 1, cols[4], .3 )
abline(v = 1:3 + .5, col = "black")
dev.off()



#Print xtable for the model
require(xtable)
xtable(sumOfModel, digits = 3)


#length(intersect(intersect(which(makesCpG == 1), which(suppDat$TypeOfSite == "nonsyn")), which(suppDat$WTnt == "a")))


#Ok, let's the our real data values for Pleuni
#Pleuni continue here. 

prepMat <- matrix(data = NA, nrow = length(1:nrow(suppDat)), ncol = length(coef.vals))
for(i in 1:nrow(suppDat)){
    intercept <- 1
    TCG <- c(0,0,0)
    toReplace <- which(suppDat[i,]$WTnt == c("t", "c", "g"))
    if(length(toReplace) > 0){
        TCG[toReplace] <- 1
    }
    changeval <- bigChange[i]
    inRTval <- as.numeric(suppDat[i,]$unit == "RT")
    nonsynVal <- as.numeric(suppDat[i,]$TypeOfSite == "nonsyn")
    shapeVal <- suppDatShape[i,]$SHAPE    
    cpgval <- makesCpG[i]
    TCGns <- TCG*nonsynVal
    tCpG <- TCG[1] * cpgval
    nonsynCpG <- cpgval*nonsynVal
    nonsynCpGandT <- cpgval*nonsynVal*TCG[1]
    prepMat[i,] <-  c(intercept, TCG, changeval, inRTval, nonsynVal, shapeVal, cpgval, TCGns, tCpG, nonsynCpG, nonsynCpGandT)
}

colnames(prepMat) <- names(coef.vals)
t(t(prepMat[1:10,]))
names(coef.vals)

obsFreqs <- suppDat[-(1:40),]$colMeansTs0
modFreqs <- exp(prepMat %*% coef.vals)[-(1:40)]
cor.test(obsFreqs,modFreqs)

#Let's try looking at synonymous mutations at ancestrally As

layout(matrix(1:2, nrow = 1))
relInds <- intersect(which(suppDat$WTnt == "a"), which(suppDat$TypeOfSite == "syn"))
cpginds <- which(makesCpG[relInds] == 1)
colsp <- rep("black", length(relInds))
colsp[cpginds] <- "red"
plot(suppDat[relInds,]$colMeansTs0, exp(prepMat[relInds,] %*% coef.vals), col = colsp)
relInds <- intersect(which(suppDat$WTnt == "t"), which(suppDat$TypeOfSite == "syn"))
cpginds <- which(makesCpG[relInds] == 1)
colsp <- rep("black", length(relInds))
colsp[cpginds] <- "red"
plot(suppDat[relInds,]$colMeansTs0, exp(prepMat[relInds,] %*% coef.vals), col = colsp)

#Maybe nonsyn mutations (big change v no?)
require(RColorBrewer)
par(mar = c(4, 4, 1, 8))
plot(0, type = "n", xlim = c(0, .5), ylim = c(0, .041), xlab = "Observed frequency", ylab = "Modeled frequency")
abline(0, 1)
nucord <- c("a", "t", "c", "g")
pcols <- brewer.pal(8, "Paired")
for(i in  1:length(nucord)){
    relInds <- intersect(which(suppDat$WTnt == nucord[i]), which(suppDat$TypeOfSite == "nonsyn"))
    cpginds <- which(makesCpG[relInds] == 1)
    bigaachanges <- which(bigChange[relInds] == 1)
    pchp <- rep(1, length(relInds))
    pchp[cpginds] <- 2
    colsp <- rep(pcols[2*i - 1], length(relInds))
    colsp[bigaachanges] <- pcols[2*i]
    points(suppDat[relInds,]$colMeansTs0, exp(prepMat[relInds,] %*% coef.vals), col = colsp, pch = pchp, xlim = c(0, .025))
}
legend("bottomright", c("A", "A (big AA)", "T", "T (big AA)", "C", "C (big AA)", "G", "G (big AA)", "non-CpG-form", "CpG-form"), col = c(pcols, "black", "black"), pch = c(rep(16, 8), 1, 2), xpd = NA)
#legend(.043, .02, c("A", "A (big AA)", "T", "T (big AA)", "C", "C (big AA)", "G", "G (big AA)", "non-CpG-form", "CpG-form"), col = c(pcols, "black", "black"), pch = c(rep(16, 8), 1, 2), xpd = NA)

modelfits <- prepMat %*% coef.vals
outmat <- cbind(suppDat, exp(modelfits))
write.table(outmat, "../dat/dat.plus.modelfits.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)


