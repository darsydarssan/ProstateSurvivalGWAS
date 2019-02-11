###################################################
#Purpose: Manhatton plots for adjusted p-values
###################################################

#Clear the work space
rm(list=ls())

#directory
getwd()
setwd()

#library
library(officer)
library(magrittr)
library(rms)
library(qqman)
library(Jmisc)
library(flextable)

#Read the data
adjpall <- NULL
for (j in 1:54){
  lhs <- paste("d","at",sep="")
  rhs <- paste("read.csv(file='adjustedpvalo",j,".csv', header=TRUE, sep=',')", sep="")
  al <- paste(paste(lhs,rhs,sep="<-"),collapse=";")
  eval(parse(text=al))
  
  adjpall <- rbind(adjpall, dat)
}

#Link the chromosome data
#Read chr data
chrdat <- read.csv(file='oncoarray_with_consortium_annotation_snpmatched.csv', header=T, sep=',')
names(chrdat)[names(chrdat)=="Onc_SNPName"] <- "snp"
chrdat <- chrdat[c("snp", "Chr")]

#Merge the chr dat with pval dat
pvalchrdat <- merge(adjpall, chrdat, by="snp", all.x=T)
#Frequency of SNPs in each chromosome
freqchrunadjp <- as.data.frame(table(pvalchrdat$Chr))
#Get rid of the chromosomes X, Y, XY, MT.
#pvalchrdat <- subset(pvalchrdat, Chr %in% c(1:23))
#Get rid of NA p-values
pvalchrdat <- subset(pvalchrdat, !is.na(pval))
pvalchrdat$pval <- ifelse(pvalchrdat$pval==0,1e-16,pvalchrdat$pval) # R put any value above 1e-320 become zero, however to have a resonable plot we use the value of 1e-16,
#-thesea are some p-values show-up in the top of the plot now. Not that important on the purpose of plot.
#Chromosome variable factor to numeric
pvalchrdat$CHR <- as.numeric(ifelse(as.character(pvalchrdat$Chr)=="MT",23,ifelse(as.character(pvalchrdat$Chr)=="X",24,
ifelse(as.character(pvalchrdat$Chr)=="XY",25,ifelse(as.character(pvalchrdat$Chr)=="Y",26,as.character(pvalchrdat$Chr))))))
pvalchrdat <- pvalchrdat[order(pvalchrdat$CHR),]
pvalchrdat$BP <- sequence(rle(pvalchrdat$CHR)$lengths)

#The manhattan plot
names(pvalchrdat)[names(pvalchrdat)=="pval"] <- "P"
names(pvalchrdat)[names(pvalchrdat)=="snp"] <- "SNP"
#Getrid of outliers
#pvalchrdat1 <- subset(pvalchrdat, P>1e-320)
manhatadjust <- tempfile(fileext = "manhatadjust.png")
png(filename=manhatadjust, width=5, height = 6, units= 'in', res=300)
manhattan(pvalchrdat, cex=0.5, cex.axis = 0.7, yaxt="n", suggestiveline=F, genomewideline = -log10(0.05/491354))
axis(2, seq(0,16,4), cex.axis=0.7)
dev.off()


#Reaed the report we are already writing
report <- read_docx(path=".docx")

#Export to my word report
report <- report %>% 
  cursor_reach(keyword="manhatadjust") %>%
  body_add_img(src = manhatadjust, width = 5, height = 6, style = "Normal")

#Printout the word output  
print(report,
      target = ".docx")


