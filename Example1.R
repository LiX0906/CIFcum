#########################################################
## R code for the data application analysis in 5.1 ######
#########################################################

rm(list=ls())
library(ggplot2)
library(cmprsk)
library(plyr)
source("Function_analysis.R")
OPdata=read.table("example1.csv",sep=",",header=TRUE)
OPdata_long=read.table("example1_long.csv",sep=",",header=TRUE)


###########################################################################################################  
## SELECT CANDIDATE CONFOUNDING FIXED COVARIATES USING UNIVARIATE FINE-GRAY MODEL

OPdata_last=ddply(OPdata,.(id), tail,1) 
sexfit=crr(OPdata_last$fup,OPdata_last$event,OPdata_last$sex,failcode=1,cencode=0)
racefit=crr(OPdata_last$fup,OPdata_last$event,OPdata_last$race,failcode=1,cencode=0)
agefit=crr(OPdata_last$fup,OPdata_last$event,OPdata_last$age,failcode=1,cencode=0)
disabilityfit=crr(OPdata_last$fup,OPdata_last$event,OPdata_last$disability,failcode=1,cencode=0)
elixhauserfit=crr(OPdata_last$fup,OPdata_last$event,OPdata_last$elixhauser,failcode=1,cencode=0)
LISfit=crr(OPdata_last$fup,OPdata_last$event,OPdata_last$LIS,failcode=1,cencode=0)


###########################################################################################################  
## FIT THE SUBDISTRIBUTION HAZARD REGRESSION MODEL WITH DIFFERENT MAXIMUM LATENCY TIME WINDOW Q

alpha=se_alpha=rep(NA,5)
est.weight=matrix(NA,5,210)
surv.base=timehaz=list()
window=c(90,120,150,180,210)

for (q in length(window)){

fit=crr_WCE(OPdata_long, lag=window[q], constrained = T,  id = "id", event = "event", start = "start", stop = "stop", expos = "dose", covariates =c("age","elixhauser","disability","LIS"), censor="IPCW")  
alpha[q]=fit$est.alpha
se_alpha[q]=sqrt(fit$vcovmat.alpha)
est.weight[q,1:window[q]]=fit$weight
surv.base[[q]]=fit$surv.base
timehaz[[q]]=fit$timehaz

}


###########################################################################################################  
##   PLOT ESTIMATED WEIGHT FUNCTION

par(mfrow=c(1,1),mar=c(5,5,3,3))
plot(est.weight[1,],type="l",xaxt='n',yaxt='n',ylim=c(-0.008,0.055),ylab="Estimated weight function",xlab="Days before last follow-up",cex=1.7,cex.lab=1.7)
for (q in 1:5){
  lines(est.weight[q,1:window[q]],col=q,lwd=2,lty=q)
}
axis(1, at=seq(0,210,30),cex.axis=1.5) 
axis(2, at=seq(-0.01,0.05,0.01),cex.axis=1.5) 
legend(170,0.053, legend=c("90 Days", "120 Days","150 Days","180 Days","210 Days"),col=1:5, lwd=2,cex=1.5,lty=1:5)
abline(h=0)



###########################################################################################################
## SIMULATE OPIOID USE PATTERNS 

maxQ=210
scenario0 <- rep(0, maxQ)  
scenario1 <- rep(c(rep(15,60),rep(0,5),rep(10,10),rep(0,30)),2)                    # exposure pattern as a function of days before last follow up
scenario2 <- rep(c(rep(20,50),rep(0,55)),2)    
scenario3 <- rep(c(rep(50,20),rep(0,85)),2)      
scenario4 <- rep(c(rep(0,60),rep(85,10),rep(0,20),rep(10,15)),2)


###########################################################################################################
# PLOT ESTIMATED CUMULATIVE INCIDENCE FUNCTION

cif1=cif2=cif3=cif4=rep(NA,5)

for (q in 1:5){
  cif1[q]=1-(1-surv.base[[q]][window[q]])^exp(alpha[q]*sum(est.weight[q,1:window[q]] * rev(scenario1[1:window[q]])))
  cif2[q]=1-(1-surv.base[[q]][window[q]])^exp(alpha[q]*sum(est.weight[q,1:window[q]] * rev(scenario2[1:window[q]])))
  cif3[q]=1-(1-surv.base[[q]][window[q]])^exp(alpha[q]*sum(est.weight[q,1:window[q]] * rev(scenario3[1:window[q]])))
  cif4[q]=1-(1-surv.base[[q]][window[q]])^exp(alpha[q]*sum(est.weight[q,1:window[q]] * rev(scenario4[1:window[q]])))
}



par(mfrow=c(2,1),bty="l",mar=c(5,5,3,3))

plot(seq(1:maxQ),scenario1,type="l",xaxt='n',yaxt='n',lwd=2,lty=1,ylab="Dose",xlab="Days from first opioid use",cex.axis=1.1, cex=1.5,cex.lab=1.3,ylim=c(0,110))
lines(seq(1:maxQ),scenario2,type="l",lwd=2,lty=2,xlab="",ylab="Dose",cex.lab=1.7,ylim=c(0,40),col="red")
lines(seq(1:maxQ),scenario3,type="l",lwd=2,lty=3,xlab="",ylab="Dose",cex.lab=1.7,ylim=c(0,40),col="green")
lines(seq(1:maxQ),scenario4,type="l", lwd=2,lty=4,xlab="Days from first opioid use",ylab="Dose",col="blue",cex.lab=1.3,ylim=c(0,40))
xlabel <- seq(0,210,by=30)
axis(1,at=xlabel,cex.axis=1.1)
axis(2,las=2)
legend(0,110,legend =c("Pattern 1","Pattern 2","Pattern 3","Pattern 4"),lty=c(1,2,3,4),lwd=c(2,2,2,2),col=c("black","red","green","blue"))



plot(cif1,lwd=1.5,pch=19, xaxt='n',yaxt='n',xlab="Days from first opioid use",ylab="Cumulative incidence",cex.axis=1.1,cex=1.7,cex.lab=1.3,ylim=c(0.02,0.08))
points(cif2,lwd=1.5,col="red", pch=7)
points(cif3,lwd=1.5,col="green", pch=24)
points(cif4,lwd=1.5,col="blue", pch=8)
xlabel2=seq(90,210,by=30)
axis(1,at=1:5,label=xlabel2,cex.axis=1.1)
axis(2,las=2,at=seq(0.02,0.08,0.01))
legend(1,0.08, legend=c("Pattern 1", "Pattern 2","Pattern 3","Pattern 4"),pch=c(19,7,24,8),
       col=1:4,lwd=1.5)


