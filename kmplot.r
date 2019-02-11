############################################################################
#Purpose:Kaplan-Mier Plot and Univariate analysis
############################################################################

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

#K-M pLot
#Read in one of those data
dat <- read.csv(file='d1.csv', header=T, sep=',')

#Calculate the time variable
dat$timeend <- (as.Date(dat$DateLastF) - as.Date(dat$DateDiag))*0.032854

#Data distribution for rms package
d <- suppressWarnings(datadist(dat)); options(datadist='d')

#survival function
survival <- Surv(time=dat$timeend, event=dat$VitalStatus, type='right')
fit <- npsurv(survival~1)
kmplot <- tempfile(fileext = "kmplot.png")
png(filename=kmplot, width=5, height = 6, units= 'in', res=300)
plot(fit, xlim=c(0,220), ylim=c(0,1.0), mark.time=F, col=c(2) ,conf.int=F, xlab='Time (days since diagnosed)',
     ylab = 'Probability of survival', lty=1, lwd=3,
     cex.axis=0.8, cex.lab=0.8, mgp=c(1.2,0.5,0))
dev.off()

#Reaed the report we are already writing
report <- 
read_docx(".docx")

#Export to my word report
report <- report %>% 
  cursor_reach(keyword="kmplot") %>%
  body_add_img(src = kmplot, width = 5, height = 6, style = "Normal")

#Find the median survival time
medsuv <- print(fit)

#Make age numeric
dat$AgeDiag <- suppressWarnings(as.numeric(as.character(dat$AgeDiag)))

#Univariate cox model for patient characteristics
modelage <- cph(survival~AgeDiag, data=dat, method="breslow",
              surv=TRUE, x=TRUE, y=TRUE)
modelpc1 <- cph(survival~pc1_euro, data=dat, method="breslow",
                surv=TRUE, x=TRUE, y=TRUE)
modelpc2 <- cph(survival~pc2_euro, data=dat, method="breslow",
                surv=TRUE, x=TRUE, y=TRUE)
modelpc3 <- cph(survival~pc3_euro, data=dat, method="breslow",
                surv=TRUE, x=TRUE, y=TRUE)
modelpc4 <- cph(survival~pc4_euro, data=dat, method="breslow",
                surv=TRUE, x=TRUE, y=TRUE)
modelpc5 <- cph(survival~pc5_euro, data=dat, method="breslow",
                surv=TRUE, x=TRUE, y=TRUE)
modelpc6 <- cph(survival~pc6_euro, data=dat, method="breslow",
                surv=TRUE, x=TRUE, y=TRUE)
modelpc7 <- cph(survival~pc7_euro, data=dat, method="breslow",
                surv=TRUE, x=TRUE, y=TRUE)

#Creating the table for univariate cox model output
covariate <- c("Age at diagnosis", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7")
hr <- round(exp(c(modelage$coef, modelpc1$coef, modelpc2$coef, modelpc3$coef,
                  modelpc4$coef, modelpc5$coef, modelpc6$coef, modelpc7$coef)),2)
lcl <- round(exp(c(modelage$coef-1.96*sqrt(modelage$var), modelpc1$coef-1.96*sqrt(modelpc1$var),
                   modelpc2$coef-1.96*sqrt(modelpc2$var), modelpc3$coef-1.96*sqrt(modelpc3$var),
                   modelpc4$coef-1.96*sqrt(modelpc4$var), modelpc5$coef-1.96*sqrt(modelpc5$var),
                   modelpc6$coef-1.96*sqrt(modelpc6$var), modelpc7$coef-1.96*sqrt(modelpc7$var))),2)
ucl <- round(exp(c(modelage$coef+1.96*sqrt(modelage$var), modelpc1$coef+1.96*sqrt(modelpc1$var),
         modelpc2$coef+1.96*sqrt(modelpc2$var), modelpc3$coef+1.96*sqrt(modelpc3$var),
         modelpc4$coef+1.96*sqrt(modelpc4$var), modelpc5$coef+1.96*sqrt(modelpc5$var),
         modelpc6$coef+1.96*sqrt(modelpc6$var), modelpc7$coef+1.96*sqrt(modelpc7$var))),2)
pval <- round(c(anova(modelage)[1,3], anova(modelpc1)[1,3], anova(modelpc2)[1,3], anova(modelpc3)[1,3], 
                anova(modelpc4)[1,3], anova(modelpc5)[1,3], anova(modelpc6)[1,3], anova(modelpc7)[1,3]),2)
          
univartab <- cbind.data.frame(covariate, hr, lcl, ucl, pval)
colnames(univartab) <- c("Covariate", "HR", "95% LCL", "95% UCL", "P-Value")

#Export to word
double_format <- function(x){sprintf("%.2f", x)}
univartabout <- regulartable(data=univartab)
univartabout <- align(univartabout, align="left", part="all") #this is not working but leave for a while
univartabout <- set_formatter(univartabout, "HR"=double_format, "95% LCL"=double_format,
                              "95% UCL"=double_format, "P-Value"=double_format)
univartabout <- theme_box(univartabout)

#Export to my word report
report <- report %>% 
  cursor_reach(keyword="TableUniVariate") %>%
  body_add_flextable(value=univartabout, align="left", split=T)

#Printout the word output  
print(report,
      target = ".docx")


