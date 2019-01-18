rm(list=ls())

# LOAD PACKAGES
library(PermAlgo)
library(survival)
library(dlnm)
library(plyr)
library(WCE)
library(kmi)
library(mgcv)

################################################################################
# SPECIFY PARAMETERS
B=3                                      # Monte carlo replicates
m=100                                         # maximum follow up days
n=400                                         # number of subjects
lag=100                                    # Maximum effective time window
qn <- qnorm(0.975)                            # nominal value
nsample=50
prob1=0.8
xrange=seq(0,10,by=0.25)                   # prediction range for x

edf=aic=matrix(NA,B,3)
colnames(edf)=colnames(aic)=c("Log link","Logit link","C-loglog link")
covtemp=rep(list(matrix(0,length(xrange),(lag+1))),3)  # Temporarily store the coverage
names(covtemp) <- c("Log link","Logit link","C-loglog link")
aictotsample <- vector("list",3)
names(aictotsample) <- c("Log link","Logit link","C-loglog link")


#source("/Users/lixingyuan/Dropbox/Simulation/Paper 2/P2_function.R")
#source("C:/Users/xil143/Dropbox/Simulation/Paper 2/P2_function.R")
source("/Users/xingyuanli/Dropbox/Simulation/Paper 2/P2_function.R")
#source("/Users/XingyuanLi/Dropbox/Simulation/Paper 2/P2_function.R")



time <- proc.time()

############################################################################
##################### DATA GENERATION ######################################
############################################################################

j=1
# GENERATE EXPOSURE TIME-DEPENDENT EXPOSURE MATRIX
expeffmat <- TDexp(n,m)



# COMPUTE THE CUMULATIVE EFFECT AT EACH TIME FOR EACH SUBJECT
cumeffmat <- t(apply(expeffmat,1,fcumeffexp,lag=c(0,lag),
                     f1=combsim[j,1],f2=combsim[j,2]))

# EXPOSURE MATRICES
Xmat=matrix(ncol=2, nrow=n*m)
Xmat <- cbind(Cov1=rep(runif(n,0,10), each=m), WCE=as.vector(t(cumeffmat)))   


    
for (i in 1:B){


# SET THE SEED
seed <- 12345 + i
set.seed(seed)




###############################################################################  
# GENERATE EVENT TIMES FROM PERMUTATION ALGORITHM: C-LOG-LOG LINK
###############################################################################

#if(FALSE){
#a=matrix(NA,30,8)
#for (j in 1:30){
#### GENERATE COMPETING RISKS DATA
eventRandom<-rep(NA,n)
v<- sapply(1:2, function(i){runif(n)})
failure <-2-rbinom(n,1,prob1)                                                            # Probability of failure on cause 1 
eventRandom[failure==1]=round((-log(1-(1/prob1)*v[failure==1,1]*prob1))*200)+1            # Latent main event time
failure[is.na(eventRandom)] <- 2
eventRandom[failure==2] <- round(-log(1-v[failure==2,2])*200)+1                           # Latent competing event time
censorRandom <- round(runif(n,80,m*2),0)                                                  # censoring time                                      


# PERMUTATION ALGORITHM TO MATCH EVENT TIMES AND COVARIATE MATRICES
data <- PAlgo_SH(n, m, Xmat, XmatNames=c( "Cov1","WCE"),eventRandom = eventRandom, censorRandom=censorRandom, 
                 failureCause =failure, betas1=c(1,1), betas2=c(log(1.1),log(1.1)))


#data$Event2=data$Event==1
#test=gam(Event2~adm.cens.exit+WCE+Cov1,data,family=binomial(link=logit),method="REML")
#test2=coxph(Surv(Start,Stop,Event==1)~WCE+Cov1,data=data)
#test3=gam(Event2~adm.cens.exit+WCE+Cov1,data,family=poisson,method="REML")
#test4=gam(Event2~s(adm.cens.exit,bs="cr",k=5)+WCE+Cov1,data,family=binomial(link=cloglog),method="REML")
#a[j,]=c(test$coefficients[3:4],test2$coefficients[1:2],test3$coefficients[3:4],test4$coefficients[2:3])}
#table(data$Event)
#apply(a,2,mean)



#if(FALSE){

###################################################################
# GENERATE EVENT TIMES FROM BINOMIAL MODEL: LOGIT LINK AND LOG LINK
###################################################################

#a=matrix(NA,100,2)
#for (j in 1:100){
t=rep(1:m,n)
baserisk=-13.3+0.01*t   
failure <-2-rbinom(n,1,prob1)  # Probability of failure on cause 1 

# Generate main event time

# log link
p=exp(baserisk[rep(failure==1,each=m)]+Xmat[rep(failure==1,each=m),2]+Xmat[rep(failure==1,each=m),1])

# logit link
#p=exp(baserisk[rep(failure==1,each=m)]+Xmat[rep(failure==1,each=m),2]+Xmat[rep(failure==1,each=m),1])/(1+exp(baserisk[rep(failure==1,each=m)]+Xmat[rep(failure==1,each=m),2]+Xmat[rep(failure==1,each=m),1]))   # logit link
u=runif(sum(failure==1)*m)
data=data.frame(cbind(Id=rep(seq(1:n)[failure==1],each=m),Stop=rep(seq(1:m),sum(failure==1)),Xmat[rep(failure==1,each=m),],p,u,Event=u<=p))
drops <- c("p","u")
data=data[ , !(names(data) %in% drops)]
mainRandom=ddply(data[data$Event==1,],.(Id),function(x) head(x,1))[,c("Id","Stop")]
colnames(mainRandom)=c("Id","Fup")
data=merge(data,mainRandom,by="Id",all=T)
data=data[order(data$Id,data$Stop),]
data$Fup[is.na(data$Fup)]=m

# Generate competing event time
data_cr=data.frame(cbind(Id=rep(seq(1:n)[failure==2],each=m),Stop=rep(seq(1:m),sum(failure==2)),Xmat[rep(failure==2,each=m),],Fup=rep(round(rexp(sum(failure==2),0.1),0)+1,each=m)))
data_cr$Event[data_cr$Stop<data_cr$Fup]=0
data_cr$Event[data_cr$Stop==data_cr$Fup]=2
data=rbind(data,data_cr)
data$Start=data$Stop-1 

# Generate non-informative censoring time
data$adm.cens.exit=rep(round(runif(n,80,2*m),0)+1,each=m)
data$Fup[data$adm.cens.exit<data$Fup]=data$adm.cens.exit[data$adm.cens.exit<data$Fup]
data <- subset(data, data$Stop <= data$Fup)
data$adm.cens.exit[! data$Event==2] = data$Stop[! data$Event==2]
data=data[order(data$Id,data$Stop),]

#test
#data$Event2=data$Event==1
#test=gam(Event2~s(adm.cens.exit,bs="cr",k=5)+WCE+Cov1,data, family=binomial(link=logit),method="REML")
#summary(test)

#table(data$Event)
#a[j,]=test$coefficients[2:3]}
#apply(a,2,mean)
#}

############################################################################
##################### DATA ANALYSIS  #######################################
############################################################################


# RC DATA: IMPUTE CENSORING TIMES FOR SUBJECTS HAVING A COMPETING EVENT
#imp<- kmi(Surv(data$Start, data$Stop, ifelse(data$Event !=0,1,0)) ~ 1, data = data,        # multiple imputation to recover potential censoring time for competing event
#               etype = data$Event, id = data$Id, failcode = 1, nimp = 5)

#data_imp <- lapply(seq_along(imp$imputed.data), function(i) {
#  datan <- imp$original.data
#  datan[, "Stop_impute"] <- imp$imputed.data[[i]][, 1]
#  datan[, "Event_impute"] <- imp$imputed.data[[i]][, 2]
#  datan })
  


##############################################################################################################
# ADMINISTRATIVE CENSORING
##############################################################################################################


# DATA RE-FORMATTING
#data_short <- data[cumsum(tapply(data$Id,data$Id,length)),c("Id","Event","Fup","Cov1","adm.cens.exit")]  # last observation of each ID
#data_short <- data[data$Stop==data$Fup,c("Id","Event","Fup","Cov1","adm.cens.exit")]

# SPLIT THE DATA
#ftime <- sort(unique(data_short$adm.cens.exit[data_short$Event==1]))
#dataspl <- survSplit(Surv(adm.cens.exit,Event==1)~., data_short, cut=ftime, start="Start")
#dataspl <- dataspl[order(dataspl$Id),]


data$event=data$Event==1

# CREATE THE MATRIX Q OF LAGGED EXPOSURES
Qlag <- do.call(rbind, lapply(seq(nrow(data)),
                              function(i) exphist(expeffmat[data$Id[i],],data$adm.cens.exit[i],c(0,lag))))


# CROSS BASIS MATRIX
Dvar <- crossbasis(Qlag,lag=lag,argvar=list(fun="ps"),arglag=list(fun="ps"))
#Dvar <- crossbasis(Qlag,lag=lag,argvar=list(fun="cr",knots=c(0,2.5,5,7.5,10)),arglag=list(fun="cr",knots=c(0,20,40,70,100,120)))

# PENALTY MATRIX
DvarPen <- cbPen(Dvar)
  
  
# FIT THE MODEL

modellist[[1]]  <- gam(event~s(adm.cens.exit,bs="cr",k=5)+Dvar+Cov1,data,family=poisson,paraPen=list(Dvar=DvarPen),method="REML")
modellist[[2]] <- gam(event~s(adm.cens.exit,bs="cr",k=5)+Dvar+Cov1,data,family=binomial(link=logit),paraPen=list(Dvar=DvarPen),method="REML")
modellist[[3]]  <- gam(event~s(adm.cens.exit,bs="cr",k=5)+Dvar+Cov1,data,family=binomial(link=cloglog),paraPen=list(Dvar=DvarPen),method="REML")  


# STORE ESTIMATION AND PREDICTION RESULTS  

    for (k in 1:3){
cpaic1 <- crosspred(Dvar,modellist[[k]],from=0,to=10,by=0.25,cen=0)

# STORE EDF AND AIC OF THE FITTED SPLINE
edf[i,k]=sum(modellist[[k]]$edf[grep("Dvar",names(modellist[[k]]$coefficients))])
aic[i,k]=modellist[[k]]$aic

#RC
#aic= mean(unlist(lapply(seq_along(imp.data$imputed.data), function(i) {result[[i]]$aic})))
#edf=mean(unlist(lapply(seq_along(imp.data$imputed.data), function(i) {sum(result[[i]]$edf[grep("Dvar",names(result[[i]]$coefficients))])})))

# STORE COVERAGE FOR THE WHOLE SURFACE
covtemp[[k]]<- covtemp[[k]]+ (effsimlist[[j]] >= cpaic1$matfit-qn*
                                    cpaic1$matse & effsimlist[[j]] <= cpaic1$matfit+qn*cpaic1$matse)


# STORE SAMPLES
  aictotsample[[k]] <- c(aictotsample[[k]],list(cpaic1$matfit))
  }
 }


proc.time()-time

#write.table(aictotsample, file="/Users/xil143/Dropbox/Simulation/aictotsample.txt")
#write.table(effsimlist, file="/Users/xil143/Dropbox/Simulation/effsimlist.txt")

# CALCULATE AMSE
AMSE=rep(0,3)
for (k in 1:3){
AMSE[k]=mean(unlist(lapply(seq(1:B),function(i) mean((aictotsample[[k]][[i]]-effsimlist[[j]])^2)/mean(effsimlist[[j]]))))
}

AMSE
unlist(lapply(covtemp,mean))/B
apply(edf,2,mean)
apply(aic,2,mean)


