########################################
#Purpose: SNP unadjusted model fit
########################################

#Clear the work space
rm(list=ls())

#directory
getwd()
setwd()

#library
library(rms)

#Read the data
for (j in 2:54){
lhs <- paste("d","at",sep="")
rhs <- paste("read.csv(file='d",j,".csv', header=TRUE, sep=',')", sep="")
al <- paste(paste(lhs,rhs,sep="<-"),collapse=";")
eval(parse(text=al))

#Survival time calculation in months
dat$timeend <- (as.Date(dat$DateLastF) - as.Date(dat$DateDiag))*0.032854

#Make the SNPs as factor
dat[2:(which(names(dat)=="CaCo")-1)] <- lapply(dat[2:(which(names(dat)=="CaCo")-1)], factor)

#Data distribution for rms package
d <- datadist(dat); options(datadist='d')

#survival function
survival <- Surv(time=dat$timeend, event=dat$VitalStatus, type='right')
#unadjusted model formula
cphform <- lapply(colnames(dat[2:(which(names(dat)=="CaCo")-1)]), function(x){as.formula(paste("survival", c(x), sep = " ~ "))})
#fitted model
model <- lapply(cphform, cph, data=dat, method="breslow", surv=TRUE, x=TRUE, y=TRUE)

#p-value of the SNP the variables
anpval <- lapply(model,function(x) if(x$'fail'==F){anova(x)[1,3]} else {NA})

#Pvalue data
snppavaldat <- cbind(names(dat[2:(which(names(dat)=="CaCo")-1)]), unlist(anpval))
colnames(snppavaldat) <- c("snp", "pval")

#list of snps with p-values
wscvd <- paste("write.csv(snppavaldat, file='unadjustedpvalo",j,".csv', row.names=F)", sep="")
eval(parse(text=wscvd))
}
