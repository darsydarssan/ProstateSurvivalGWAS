#####################################################################
#Purpose: SNP data split
#####################################################################

##Read SNP data
geno <- read.delim(file ="QLD_onco_geno.txt",header=TRUE, sep = "\t", dec = ".")
names(geno)[names(geno) == 'SNP'] <- 'Onc_ID'
##geno data description
dimgeno <- dim(geno) #dimension of the data - answer 4618, 533633
patidgeno <- geno[,1] #patient id
namesgeno <- names(geno) #SNP names

#Split the data into 10000 variables (SNPs) per data set
geno1 <- geno[,c(1:10000)]
lhs1 <- paste("geno",2:53,sep="")
rhs1 <- paste("geno[,c(1,",seq(10001,,by=10000,length(lhs1)),":",seq(20000,,by=10000,length(lhs1)),")]",sep="")
al1 <- paste(paste(lhs1,rhs1,sep="<-"),collapse=";")
eval(parse(text=al1))
geno54 <- geno[,c(1,530001:533633)]

#Export the splitted data sets into csv
for (i in 1:54){
write.csv(eval(as.symbol(paste0("geno",i,sep=""))), file=paste0("geno", i, ".csv", sep=""), row.names=FALSE)
}