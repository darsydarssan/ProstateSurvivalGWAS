#####################################################################
#Purpose: Adjusted cox model p-values
####################################################################

#Clear the work space
rm(list=ls())

#directory
getwd()
setwd("")

#library
library(rms)

#Read the data
for (j in 1:54){
lhs <- c("dat")
rhs <- paste("read.csv(file='d",j,".csv', header=TRUE, sep=',')", sep="")
al <- paste(paste(lhs,rhs,sep="<-"),collapse=";")
eval(parse(text=al))

#Survival time calculation in months
dat$timeend <- (as.Date(dat$DateLastF) - as.Date(dat$DateDiag))*0.032854

#Make the SNPs as factor
dat[2:(which(names(dat)=="CaCo")-1)] <- lapply(dat[2:(which(names(dat)=="CaCo")-1)], factor)
dat$AgeDiag <- suppressWarnings(as.numeric(as.character(dat$AgeDiag)))

#Data distribution for rms package
d <- suppressWarnings(datadist(dat)); options(datadist='d')

#survival function
survival <- Surv(time=dat$timeend, event=dat$VitalStatus, type='right')
#adjusted model formula
cphform <- lapply(colnames(dat[2:(which(names(dat)=="CaCo")-1)]), 
function(x){as.formula(paste("survival", paste(c(c(x), c("AgeDiag", "pc1_euro", "pc2_euro",	
"pc3_euro",	"pc4_euro",	"pc5_euro",	"pc6_euro",	"pc7_euro")), sep="", collapse= " + "), sep = " ~ "))})
#fitted model
model <- lapply(cphform, cph, data=dat, method="breslow", surv=TRUE, x=TRUE, y=TRUE)

#p-value of the SNP the variables
anpval <- lapply(model,function(x) if(x$'fail'==F){anova(x)[1,3]} else {NA})

#Pvalue data
snppavaldat <- cbind(names(dat[2:(which(names(dat)=="CaCo")-1)]), unlist(anpval))
colnames(snppavaldat) <- c("snp", "pval")

#list of snps with p-values
wscvd <- paste("write.csv(snppavaldat, file='adjustedpvalo",j,".csv', row.names=F)", sep="")
eval(parse(text=wscvd))
rm(dat, d, cphform, anpval, model, snppavaldat)
}
