rm(list=ls())
library(ggplot2)
######################################################################################
############### SPECIFY PARAMETERS  ##################################################
######################################################################################


m=100                                      # maximum follow up days                                        # number of subjects
lag=100                                    # Maximum latency time window

source("/Users/lixingyuan/Dropbox/Simulation/Paper 2/P2_function.R")

var3d <- c("6","6","6")
var3d2 <- c(6,6,6)
lag3d <- paste("lag",c(40,40,40),sep="")
lag3d2 <- c(30,30,30)


var3d <- c("8","7","8")
var3d2 <- c(8,7,8)
lag3d <- paste("lag",c(20,20,15),sep="")
lag3d2 <- c(20,20,15)

ylimhigh <- vector("list",3)
ylimhigh[[1]] <- rep(1.3,4)
ylimhigh[[2]] <- rep(1.6,4)
ylimhigh[[3]] <- rep(c(1.5,1.6),c(2,2))



######################################################################################
# Figure 1. Plot true 3D cumulative effect of exposure #################################
######################################################################################


# ARGUMENTS FOR 3D PLOTS
arg3D <- list(x=seq(0,10,0.25),y=0:100,ticktype="detailed",theta=230,
              ltheta=200,phi=30,lphi=30,xlab="Exposure", ylab="Lag",zlab="RR",
              shade = 0.75,r=sqrt(3),d=5,cex.axis=1.5,cex.lab=1.5,border=grey(0.3),
                          col=grey(0.99))


png(filename="/Users/lixingyuan/Dropbox/meeting/P2_simtrue2.png",width=1000, height=450)
par(mfrow=c(1,3),mar=c(2,2,2,2))
for (j in 1:3){ 
  d3 <- do.call(persp,modifyList(arg3D,list(z=exp(effsimlist[[j]]))))
  title(rownames(combsim)[j],cex.main=2)
  lines (trans3d(x=var3d2[j],y=0:lag,z=exp(effsimlist[[j]])[var3d[j],],pmat=d3),
         col=1,lwd=2)
  lines (trans3d(x=0:40/4,y=lag3d2[j],z=exp(effsimlist[[j]])[,lag3d[j]],pmat=d3),
         col=1,lwd=2)
}
dev.off()




#############################################################################################
# Figure 2. PLOT ESTIMATED EXPOSURE EFFECT CURVE AND LAG EFFECT CURVE(30 SAMPLES)   ########
#############################################################################################

# Linear-Constant

j=1
source("/Users/xingyuanli/Dropbox/Simulation/Paper 2/P2_function.R")
load("/Users/xingyuanli/Dropbox/Simulation/Output/Paper 2/2018.12.11 REML no df/B=250/aictotsample_logit_j1_n600_REML.Rdata")

png(filename="/Users/xingyuanli/Dropbox/meeting/P2_simest_logit_scenario1.png",width=900, height=600)

layout(matrix(1:15,5),heights=c(1,rep(0.8,4)))
par(mfrow=c(2,3),mar=c(5,5,5,5))
yname=c("RR","OR","SHR")

for (k in 1:3){ 
  
  # lag effect 
  
  plot(0:lag,exp(effsimlist[[j]][var3d[j],]),type="n",xlab="Lag",ylab=paste(yname[k]),
       frame.plot=FALSE,cex=1.5,cex.lab=1.5,ylim=c(0.7,1.25))
  for(m in aictotsample[[k]][51:80]) lines(0:lag,exp(m[var3d[j],]),
                                    lwd=1.5,col=grey(0.8))
  #abline(h=1)
  lines(0:lag,exp(effsimlist[[j]][var3d[j],]),lwd=1.5)
  legend("top",c("True","Estimated"),lty=c(1,1),cex=1.5,
         col=c(1,grey(0.7)),horiz=T,bty="n",x.intersp=0.1,inset=0.03)
  title(names(aictotsample)[[k]],cex.main=2)
}  


for (k in 1:3){   
  # dose effect
  
  plot(seq(0,10,0.25),exp(effsimlist[[j]][,lag3d[j]]),type="n",xlab="Exposure", ylab=paste(yname[k]),
       frame.plot=FALSE,cex=1.5,cex.lab=1.5,ylim=c(1,1.12))
  for( m in aictotsample[[j]][1:30]) lines(seq(0,10,0.25),exp(m[,lag3d[j]]),
                                     lwd=1.5,col=grey(0.8))  
  #abline(h=1)
  lines(seq(0,10,0.25),exp(effsimlist[[j]][,lag3d[j]]),lwd=1.5)
  legend("top",c("True","Estimated"),lty=c(1,1),cex=1.5,
         col=c(1,grey(0.7)),horiz=T,bty="n",x.intersp=0.2,inset=0.03)
  
}

dev.off()

# Exponential-Decay

j=2
source("/Users/xingyuanli/Dropbox/Simulation/Paper 2/P2_function.R")
load("/Users/xingyuanli/Dropbox/Simulation/Output/Paper 2/2018.12.11 REML no df/B=250/aictotsample_logit_j2_n600_REML.Rdata")

png(filename="/Users/xingyuanli/Dropbox/meeting/P2_simest_logit_scenario2.png",width=900, height=600)

layout(matrix(1:15,5),heights=c(1,rep(0.8,4)))
par(mfrow=c(2,3),mar=c(5,5,5,5))
yname=c("RR","OR","SHR")

for (k in 1:3){ 
  
  # lag effect 
  
  plot(0:lag,exp(effsimlist[[j]][var3d[j],]),type="n",xlab="Lag",ylab=paste(yname[k]),
       frame.plot=FALSE,cex=1.5,cex.lab=1.5,ylim=c(1,1.4))
  for(m in aictotsample[[k]][51:80]) lines(0:lag,exp(m[var3d[j],]),
                                           lwd=1.5,col=grey(0.8))
  #abline(h=1)
  lines(0:lag,exp(effsimlist[[j]][var3d[j],]),lwd=1.5)
  legend("top",c("True","Estimated"),lty=c(1,1),cex=1.5,
         col=c(1,grey(0.7)),horiz=T,bty="n",x.intersp=0.1,inset=0.03)
}  


for (k in 1:3){   
  # dose effect
  
  plot(seq(0,10,0.25),exp(effsimlist[[j]][,lag3d[j]]),type="n",xlab="Exposure", ylab=paste(yname[k]),
       frame.plot=FALSE,cex=1.5,cex.lab=1.5,ylim=c(1,1.4))
  for( m in aictotsample[[j]][1:30]) lines(seq(0,10,0.25),exp(m[,lag3d[j]]),
                                           lwd=1.5,col=grey(0.8))  
  #abline(h=1)
  lines(seq(0,10,0.25),exp(effsimlist[[j]][,lag3d[j]]),lwd=1.5)
  legend("top",c("True","Estimated"),lty=c(1,1),cex=1.5,
         col=c(1,grey(0.7)),horiz=T,bty="n",x.intersp=0.2,inset=0.03)
  
}

dev.off()

## Plateau-Peak

j=3
source("/Users/xingyuanli/Dropbox/Simulation/Paper 2/P2_function.R")
load("/Users/xingyuanli/Dropbox/Simulation/Output/Paper 2/2018.12.11 REML no df/B=250/aictotsample_logit_j3_n600_REML.Rdata")

png(filename="/Users/xingyuanli/Dropbox/meeting/P2_simest_logit_scenario3.png",width=900, height=600)

layout(matrix(1:15,5),heights=c(1,rep(0.8,4)))
par(mfrow=c(2,3),mar=c(5,5,5,5))
yname=c("RR","OR","SHR")

for (k in 1:3){ 
  
  # lag effect 
  
  plot(0:lag,exp(effsimlist[[j]][var3d[j],]),type="n",xlab="Lag",ylab=paste(yname[k]),
       frame.plot=FALSE,cex=1.5,cex.lab=1.5,ylim=c(1,1.25))
  for(m in aictotsample[[k]][51:80]) lines(0:lag,exp(m[var3d[j],]),
                                           lwd=1.5,col=grey(0.8))
  #abline(h=1)
  lines(0:lag,exp(effsimlist[[j]][var3d[j],]),lwd=1.5)
  legend("top",c("True","Estimated"),lty=c(1,1),cex=1.5,
         col=c(1,grey(0.7)),horiz=T,bty="n",x.intersp=0.1,inset=0.03)
}  


for (k in 1:3){   
  # dose effect
  
  plot(seq(0,10,0.25),exp(effsimlist[[j]][,lag3d[j]]),type="n",xlab="Exposure", ylab=paste(yname[k]),
       frame.plot=FALSE,cex=1.5,cex.lab=1.5,ylim=c(1,1.2))
  for( m in aictotsample[[j]][1:30]) lines(seq(0,10,0.25),exp(m[,lag3d[j]]),
                                           lwd=1.5,col=grey(0.8))  
  #abline(h=1)
  lines(seq(0,10,0.25),exp(effsimlist[[j]][,lag3d[j]]),lwd=1.5)
  legend("top",c("True","Estimated"),lty=c(1,1),cex=1.5,
         col=c(1,grey(0.7)),horiz=T,bty="n",x.intersp=0.2,inset=0.03)
  
}

dev.off()



########################################################
# Figure 3. PLOT AMSE AND COVERAGE   ###################
########################################################

link=c("cloglog","logit","log")
scenario=c("j1","j2","j3")
samplesize=c("150","200","250","400","600")
AMSE=coverage=list(matrix(NA, 9, 5), matrix(NA, 9, 5),matrix(NA, 9, 5))
names(AMSE)=names(coverage)=link=c("cloglog","logit","log")

for (k in 1:3){
  for (j in 1:3){
    for (n in 1:5){
result=read.table(paste("/Users/lixingyuan/Dropbox/Simulation/Output/Paper 2/2018.12.11 REML no df/B=250/result_",link[k],"_",scenario[j],"_n",samplesize[n],"_REML.txt",sep=""))
AMSE[[k]][((3*j-2):(3*j)),n]=result[,1]
coverage[[k]][((3*j-2):(3*j)),n]=result[,2]
    }
  }
}

# AMSE
png(filename="/Users/lixingyuan/Dropbox/meeting/P2_AMSE_2.png",width=900, height=700)

par(mfrow=c(3,3),mar=c(4,4,4,4))

for (k in 1:3){
plot(AMSE[[k]][1,],type="l",col=1,xlab="Sample size",xaxt="n",yaxt="n",ylab="AMSE",cex=1.5,cex.lab=1.5,lwd=2,ylim=c(0,1),main="Constant-Linear")
axis(1, at=c(1:5),labels=samplesize,cex.lab=1.3)
axis(2, las=1,cex=1.3)
lines(AMSE[[k]][2,],col=2,lwd=2)
lines(AMSE[[k]][3,],col=3,lwd=2)
legend("topright",c("Log link","Logit link","C-loglog link"),lty=c(1,1,1),
       col=c(1,2,3),bty="n",x.intersp=0.2,inset=0.03)

plot(AMSE[[k]][4,],type="l",col=1,xlab="Sample size",xaxt="n",yaxt="n",ylab="AMSE",cex=1.5,cex.lab=1.5,lwd=2,ylim=c(0,0.7),main="Decay-Exponential")
axis(1, at=c(1:5),labels=samplesize)
axis(2, las=1)
lines(AMSE[[k]][5,],col=2,lwd=2)
lines(AMSE[[k]][6,],col=3,lwd=2)
legend("topright",c("Log link","Logit link","C-loglog link"),lty=c(1,1,1),
       col=c(1,2,3),bty="n",x.intersp=0.2,inset=0.03)

plot(AMSE[[k]][7,],type="l",col=1,xlab="Sample size",xaxt="n",yaxt="n",ylab="AMSE",cex=1.5,cex.lab=1.5,lwd=2,ylim=c(0,3),main="Peak-Plateau")
axis(1, at=c(1:5),labels=samplesize) 
axis(2, las=1)
lines(AMSE[[k]][8,],col=2,lwd=2)
lines(AMSE[[k]][9,],col=3,lwd=2)
legend("topright",c("Log link","Logit link","C-loglog link"),lty=c(1,1,1),
       col=c(1,2,3),bty="n",x.intersp=0.2,inset=0.03)
}
dev.off()

# coverage
png(filename="/Users/lixingyuan/Dropbox/meeting/P2_coverage_2.png",width=900, height=700)

par(mfrow=c(3,3),mar=c(4,4,4,4))

for (k in 1:3){
  plot(coverage[[k]][1,],type="l",col=1,xlab="Sample size",xaxt="n",yaxt="n",ylab="Coverage probability",cex=1.5,cex.lab=1.5,lwd=2,ylim=c(0.88,1),main="Constant-Linear")
  axis(1, at=c(1:5),labels=samplesize) 
  axis(2, las=1)
  lines(coverage[[k]][2,],col=2,lwd=2)
  lines(coverage[[k]][3,],col=3,lwd=2)
  abline(h=0.95,col="grey",lty=2)
  legend("topright",c("Log link","Logit link","C-loglog link"),lty=c(1,1,1),
         col=c(1,2,3),bty="n",x.intersp=0.2,inset=0.03)
 
  
  plot(coverage[[k]][4,],type="l",col=1,xlab="Sample size",ylab="Coverage probability",xaxt="n",yaxt="n",cex=1.5,cex.lab=1.5,lwd=2,ylim=c(0.88,1),main="Decay-Exponential")
  axis(1, at=c(1:5),labels=samplesize) 
  axis(2, las=1)
  lines(coverage[[k]][5,],col=2,lwd=2)
  lines(coverage[[k]][6,],col=3,lwd=2)
  abline(h=0.95,col="grey",lty=2)
  legend("topright",c("Log link","Logit link","C-loglog link"),lty=c(1,1,1),
         col=c(1,2,3),bty="n",x.intersp=0.2,inset=0.03)
  
  plot(coverage[[k]][7,],type="l",col=1,xlab="Sample size",ylab="Coverage probability",xaxt="n",yaxt="n",cex=1.5,cex.lab=1.5,lwd=2,ylim=c(0.85,1),main="Peak-Plateau")
  axis(1, at=c(1:5),labels=samplesize)
  axis(2, las=1)
  lines(coverage[[k]][8,],col=2,lwd=2)
  lines(coverage[[k]][9,],col=3,lwd=2)
  abline(h=0.95,col="grey",lty=2)
  legend("topright",c("Log link","Logit link","C-loglog link"),lty=c(1,1,1),
         col=c(1,2,3),bty="n",x.intersp=0.2,inset=0.03)
}
dev.off()


############################################################
# Figure 4. PLOT ESTIMATE COVERAGE SURFACE  (coverage.txt) ####
############################################################

