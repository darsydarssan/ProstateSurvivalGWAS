###################################################
#Purpose:unadjusted significant SNPs p-values
###################################################

#Identify the SNP factor which has p-value less than the threshold. Find HRs CIs and p-values of those SNPs

#Clear the work space
rm(list=ls())

#directory
getwd()
setwd()

#library
library(rms)
library(qqman)
library(Jmisc)
library(genetics)
library(plyr)
library(magrittr)
library(officer)
library(flextable)

#Read the data
sigsnpall <- NULL
for (j in 1:54){
lhs <- paste("d","at",sep="")
rhs <- paste("read.csv(file='unadjustedpvalo",j,".csv', header=TRUE, sep=',')", sep="")
al <- paste(paste(lhs,rhs,sep="<-"),collapse=";")
eval(parse(text=al))

#p-value threshold
#significant snps #491354 SNPs-tests, 0.05 is the type I error.
sigsnp <- dat[which(dat$pval < (0.05/491354)),] 
if(dim(sigsnp)[1] != 0){sigsnp <- addCol(sigsnp,data=j)}
if(dim(sigsnp)[1] != 0){sigsnpall <- rbind(sigsnpall, sigsnp)}
}

#Getrid of list for the data set each snp coming from
datneed <- unname(unlist(sigsnpall$data))

coefall <- NULL
for (i in unique(datneed)){
lhs <- paste("d","at",sep="")
rhs <- paste("read.csv(file='d",i,".csv', header=TRUE, sep=',')", sep="")
al <- paste(paste(lhs,rhs,sep="<-"),collapse=";")
eval(parse(text=al))

#subset the orginal data where the snps are available
snpnam <- subset(sigsnpall, data==i)
datsnp <- dat[c("Onc_ID", as.character(snpnam$snp), "DateDiag", "DateLastF", "VitalStatus")]
#Survival time calculation in months
datsnp$timeend <- (as.Date(datsnp$DateLastF) - as.Date(datsnp$DateDiag))*0.032854

#run the model on those snps
#Make the SNPs as factor
datsnp[2:(which(names(datsnp)=="DateDiag")-1)] <- lapply(datsnp[2:(which(names(datsnp)=="DateDiag")-1)], factor)

##Data distribution for rms package
d <- suppressWarnings(datadist(datsnp)); options(datadist='d')

#survival function
survival <- Surv(time=datsnp$timeend, event=datsnp$VitalStatus, type='right')
#unadjusted model formula
cphform <- lapply(colnames(datsnp[2:(which(names(datsnp)=="DateDiag")-1)]), function(x){as.formula(paste("survival", c(x), sep = " ~ "))})
#fitted model
model <- lapply(cphform, cph, data=datsnp, method="breslow", surv=TRUE, x=TRUE, y=TRUE)

coef1 <- t(t(unlist(lapply(model, function(x){x$coef}))))
se1 <- t(t(unlist(lapply(model, function(x){sqrt(diag(x$var))}))))
pval1 <- t(t(unlist(lapply(model, function(x){pnorm(abs(x$coef/sqrt(diag(x$var))),lower.tail=F)*2})))) 
coef1 <- cbind(coef1, se1, pval1)
coef1 <- addCol(coef1, data=i)
coefall <- rbind(coefall, coef1)
}
colnames(coefall) <- c("Coef", "stderr", "Pval", "data")

#Find the significant allelle name
#This could have been done in the previous loop but we do it here first and put it in the above loop later if necessary
snpnamall <- levdatall <- genotypeall <- freqall <- NULL
for (i in unique(datneed)){
lhs <- paste("d","at",sep="")
rhs <- paste("read.csv(file='d",i,".csv', header=TRUE, sep=',')", sep="")
al <- paste(paste(lhs,rhs,sep="<-"),collapse=";")
eval(parse(text=al))

lhs1 <- paste("ge","no",sep="")
rhs1 <- paste("read.csv(file='geno",i,".csv', header=TRUE, sep=',')", sep="")
al1 <- paste(paste(lhs1,rhs1,sep="<-"),collapse=";")
eval(parse(text=al1))

snpdat <- subset(sigsnpall, data==i)
datsnp <- dat[c("Onc_ID", as.character(snpdat$snp))]
genosnp <- geno[c("Onc_ID", as.character(snpdat$snp))]

if(dim(datsnp)[2]>2){
levdat <- lapply(datsnp[,-1], function(x){levels(factor(x))})
snpnam <- rep(names(datsnp)[-1], lengths(levdat))
} else {
levdat <- levels(factor(datsnp[,-1]))
snpnam <- rep(names(datsnp)[-1], length(levdat))
}

datsnp1 <- datsnp[,1:2]
colnames(datsnp1) <- c("Onc_ID", "noneed")
datgeno <- merge(genosnp, datsnp1, by="Onc_ID", all.y=TRUE)
datgeno <- datgeno[,-c(dim(datgeno)[2])]

if (dim(datgeno)[2] > 2){
typing <- lapply(datgeno[,-1],genotype,sep="")
sumgeno <- lapply(typing, summary) 
genotype <- unlist(lapply(sumgeno, function(x){row.names(x$genotype.freq)[row.names(x$genotype.freq)!="-/-"]}), use.names=F)
freq <- unlist(lapply(sumgeno, function(x){x$genotype.freq[which(row.names(x$genotype.freq)!="-/-"),2]}), use.names=F)
} else {
typing <- genotype(datgeno[,-1],sep="")
sumgeno <- summary(typing)
genotype <- row.names(sumgeno$genotype.freq)[row.names(sumgeno$genotype.freq)!="-/-"]
freq <- sumgeno$genotype.freq[which(row.names(sumgeno$genotype.freq)!="-/-"),2]
}

snpnamall <- c(snpnamall, snpnam)
levdatall <- c(levdatall, unlist(levdat, use.names=F))
#levgenoall <- c(levgenoall, unlist(levgeno, use.names=F)) - cannot use this
genotypeall <- c(genotypeall, genotype)
freqall <- c(freqall, freq)
}

finaldat <- as.data.frame(cbind(snpnamall, levdatall, genotypeall, freqall))

#Merge finaldat and coefall
coefall1 <- cbind(row.names(coefall), coefall)
colnames(coefall1)[1] <- "snpnamall1"
rownames(coefall1) <- NULL
coefall1 <- as.data.frame(coefall1)
coefall1$snpnamall <- substr(coefall1$snpnamall1,1,nchar(coefall1$snpnamall1)-2)
coefall1$coefnum <- substr(coefall1$snpnamall1,nchar(coefall1$snpnamall1), nchar(coefall1$snpnamall1))
coefall1$mnam <- paste(coefall1$snpnamall, coefall1$coefnum, sep="")
finaldat$mnam <- paste(finaldat$snpnamall, finaldat$levdatall, sep="")
fincoef <- merge(finaldat, coefall1, by="mnam", all.x=T)
fincoef <- fincoef[c("snpnamall.x", "genotypeall", "freqall", "Coef", "stderr", "Pval")]
names(fincoef)[names(fincoef)=="snpnamall.x"] <- "snpnam"

#Add chromosome details
chrdat <- read.csv(file='oncoarray_with_consortium_annotation_snpmatched.csv', header=T, sep=',')
names(chrdat)[names(chrdat)=="Onc_SNPName"] <- "snpnam"
chrdat <- chrdat[c("snpnam", "Chr")]
chrdat1 <- merge(fincoef, chrdat, by="snpnam", all.x=T)

#Merge with the overall pvalue
names(sigsnpall)[names(sigsnpall)=="snp"] <- "snpnam"
names(sigsnpall)[names(sigsnpall)=="pval"] <- "PvalOverall"
finchrdat <- merge(chrdat1, sigsnpall, by="snpnam", all.x=T)

#Find HR and CI
finchrdat$hr <- round(exp(as.numeric(finchrdat$Coef)),2)
finchrdat$lci <- round(exp(as.numeric(finchrdat$Coef) - 1.96*as.numeric(finchrdat$stderr)),2)
finchrdat$uci <- round(exp(as.numeric(finchrdat$Coef) + 1.96*as.numeric(finchrdat$stderr)),2)
finchrdat$hrci <- paste("[",finchrdat$lci," ",finchrdat$uci,"]",sep="")

#Reorganize the table for the output
unadjsnppval <- data.frame(ChrNo=finchrdat$Chr, SNPName=finchrdat$snpnam, Genotype=finchrdat$genotypeall, Frequency=round(as.numeric(as.character(finchrdat$freqall)),4), 
                              HR=finchrdat$hr, CI=finchrdat$hrci, pval=finchrdat$PvalOverall)

#Reaed the report we are already writing
report <- read_docx(path=".docx")

#Export to word
double_format <- function(x){sprintf("%.2f", x)}
pval_format <- function(x){sprintf("%.2e", x)}
unadjsnppvalout <- regulartable(data=unadjsnppval)
unadjsnppvalout <- align(unadjsnppvalout, align="left", part="all") #this is not working but leave for a while
unadjsnppvalout <- set_formatter(unadjsnppvalout, "HR"=double_format, "pval"=pval_format)
unadjsnppvalout <- theme_box(unadjsnppvalout)

#Export to my word report
report <- report %>% 
  cursor_reach(keyword="Unadjustedpval") %>%
  body_add_flextable(value=unadjsnppvalout, align="left", split=T)

#Printout the word output  
print(report,
      target = ".docx")










