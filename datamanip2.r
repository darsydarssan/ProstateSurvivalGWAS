###############################################################
#Purpose: Prostate Cancer SNP Analysis - Data manipulation
###############################################################

#Clear work space
rm(list=ls())

getwd()
setwd()
##Read constant snp data - remember we have 54 files
samephenormall <- snpmis90all <- c()
for (j in 1:54){
lhs <- paste("conssnpnames",j,sep="")
rhs <- paste("read.csv('","conssnpnames",j,".csv')",sep="")
al <- paste(paste(lhs,rhs,sep="<-"),collapse=";")
eval(parse(text=al))
samephenormall <- rbind(samephenormall, eval(parse(text=lhs)))

lhs1 <- paste("dropsnps",j,sep="")
rhs1 <- paste("read.csv('","dropsnps",j,".csv')",sep="")
al1 <- paste(paste(lhs1,rhs1,sep="<-"),collapse=";")
eval(parse(text=al1))
snpmis90all <- rbind(snpmis90all, eval(parse(text=lhs1)))
}

#list of snps with constant phenotypes
write.csv(samephenormall, file="samephenormall.csv", row.names=FALSE)
#list of snps wih more than 90% missing
write.csv(snpmis90all, file="snpmis90all.csv", row.names=FALSE)

