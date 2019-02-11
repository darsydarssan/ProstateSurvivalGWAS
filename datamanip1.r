###############################################################
#Purpose: Data manipulation
###############################################################

#Clear work space
rm(list=ls())

getwd()
setwd()

#library
library(dataPreparation)
library(genetics)

#Read prognestic data
clin <- read.csv(file="QLDpatientswithOncoID.csv", header=TRUE, sep=",")
#only keep the vaiables we need
clin1 <- clin[,c("Onc_ID", "CaCo", "DateDiag", "DateLastF", "AgeDiag", "VitalStatus")]

#Read the PC data
pcdat <- read.csv(file="QLD_release_Nov_2015.csv", header=TRUE, sep=",")
#only keep the variables we need
pcdat1 <- pcdat[,c("Onc_ID", "Geno_Ancestry", "pc1_euro",	"pc2_euro",	"pc3_euro",	"pc4_euro",	"pc5_euro",	"pc6_euro",	"pc7_euro")]

##Read SNP data - remember we have 54 geno data
for (j in 1:54){
lhs <- paste("geno",j,sep="")
rhs <- paste("read.csv('","geno",j,".csv')",sep="")
al <- paste(paste(lhs,rhs,sep="<-"),collapse=";")
eval(parse(text=al))

#Merge the data geno and clin
lhs1 <- paste("geno","clin1",sep="")
rhs1 <- paste("merge(geno",j,", clin1, ", "by='Onc_ID', ", "all=T)", sep="")
al1 <- paste(paste(lhs1,rhs1,sep="<-"),collapse=";")
eval(parse(text=al1))

#Merge all three data
dat <- merge(genoclin1, pcdat1, by="Onc_ID", all=T)

##Data manipulation
#Remove non-europian ancestry
dateuro <- dat[which(dat$Geno_Ancestry=='European'), ]
#Remove the patients who's vital status is unknown
dateuronovit <- dateuro[which(dateuro$VitalStatus!='U'), ]
#Remove the control patients
dateuronovitcases <- dateuronovit[which(dateuronovit$CaCo==1),]

#Get rid of the SNPs which has more than 90% missing data
#Names of the SNPs we got rid of
dropsnps <- c(names(which(colMeans(dateuronovitcases=="--") > 0.9)==T), 
names(which(colMeans(dateuronovitcases=="") > 0.9)==T),
names(which(colMeans(is.na(dateuronovitcases)) > 0.9)==T))
dropsnploc <- which(colnames(dateuronovitcases) %in% c(dropsnps))
d <- dateuronovitcases[,-dropsnploc]

#Replace all the "--" in SNPs by NA
d[d=="--"] <- NA
#make the SNP as factor
factsnp <- lapply(d[,-c(1,(ncol(d)-12):ncol(d))], unlist, use.names=FALSE) 
#apply the genotype function
genosnp <- lapply(factsnp,genotype,sep="")  
#summary of genotype function to get allelle frequency
sumgeno <- lapply(genosnp, summary) 
#name of allelle without NA
namallele <- lapply(sumgeno, function(x){row.names(x$"allele.freq")[which(row.names(x$"allele.freq") != "NA")] })

#encode as per allele frequency
for(i in 1:dim(d[,-c(1,(ncol(d)-12):ncol(d))])[2]){
  if (length(namallele[[i]])==1){
    d[,i+1] <- 0
  } else if (length(namallele[[i]])==2){
    d[,i+1] <- ifelse(unlist(d[,i+1])==paste0(namallele[[i]][1],namallele[[i]][1]),0,
                         ifelse(unlist(d[,i+1])==paste0(namallele[[i]][1],namallele[[i]][2]) |
                                  unlist(d[,i+1])==paste0(namallele[[i]][2],namallele[[i]][1]),1,
                                ifelse(unlist(d[,i+1])==paste0(namallele[[i]][2],namallele[[i]][2]),2,
                                       d[,i+1])))
  }
}

#Getrid of the SNP they have same phenotype among patients
conssnpnam <- whichAreConstant(dataSet=d[,-c(1,(ncol(d)-12):ncol(d))], verbose=F)
conssnpnames <- names(d)[conssnpnam+1]
d <- d[,-c(conssnpnam+1)]

wscvd <- paste("write.csv(d, file='d",j,".csv', row.names=F)", sep="")
eval(parse(text=wscvd))
wscvdropsnp90 <- paste("write.csv(dropsnps, file='dropsnps",j,".csv', row.names=F)", sep="")
eval(parse(text=wscvdropsnp90))
wscvconssnp <- paste("write.csv(conssnpnames, file='conssnpnames",j,".csv', row.names=F)", sep="")
eval(parse(text=wscvconssnp))
rmgeno <- paste("rm(geno",j,")",sep="")
eval(parse(text=rmgeno))
rm(d)
}
