#setwd("/Users/pleuni/Documents/Git/bachelerProject")
source("Rscripts/baseRscript.R")

OverviewDFZanini <- read.table("Output/OverviewSelCoeffZanini.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

head(OverviewDFZanini)

vec<-c()
vec<-c(vec, mean(OverviewDFZanini$colMeansTsZanini[seq(from=1,to= nrow(OverviewDFZanini), by=3)]))
vec<-c(vec, mean(OverviewDFZanini$colMeansTsZanini[seq(from=2,to= nrow(OverviewDFZanini), by=3)]))
vec<-c(vec, mean(OverviewDFZanini$colMeansTsZanini[seq(from=3,to= nrow(OverviewDFZanini), by=3)]))

vec<-c(vec, mean(OverviewDFZanini$colMeansTsZanini[OverviewDFZanini$TypeOfSite=="nonsyn"]))
vec<-c(vec, mean(OverviewDFZanini$colMeansTsZanini[OverviewDFZanini$TypeOfSite=="syn"]))

barplot(vec, names.arg = c("first", "second", "third", "nonsyn", "syn"), main="mean freq of transition mutations in Pol in Zanini data")
