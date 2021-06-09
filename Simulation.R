rm(list=ls())
library(survival)
library(parallel)

source("Function_data generation.R")
source("Function_analysis.R")


################################################################################
# SPECIFY PARAMETERS

B=3                                        # Monte carlo replicates
lag=120                                       # Maximum effective time window
m=365                                         # maximum follow up days
n=500                                         # number of subjects
beta1=c(1.5,4)                                # alpha
prob1=0.8                                    # probability of type 1 failure
elapsed=seq(0,m)
censor_list=c("AC","IPCW","MI")

# CREATE LISTS AND VECTORS TO SAVE RESULTS
edf=CP=lapply(1:3,matrix,data=NA,nrow=4,ncol=B)
est.weight_ac=est.weight_ipcw=est.weight_mi=lapply(1:4, matrix, data= NA, nrow=B, ncol=lag)
names(edf)=names(CP)=censor_list


# SIMULATE WEIGHT FUNCTION (STANDARDIZED)
weight_mat=rbind(
  weight1=ifelse(elapsed<lag,1/lag,0),                                            # Constant 
  weight2=(1-1/(1+exp(-(elapsed-70)/12)))/70,                                     # S-shape decay
  weight3=ifelse(elapsed<lag,dexp(0:lag,0.035),0),                                # Exponential decay
  weight4=ifelse(elapsed<lag,dnorm(0:lag,30,10)/2+ dnorm(0:lag,85,25)/1.82,0))    # Peak

# GENERATE EXPOSURE TIME-DEPENDENT EXPOSURE MATRIX
dose <- TDexp(n,m)

time0=proc.time()

for (q in 3){ 
  ######################################################################################################################### 
  
  # CALCULATE CUMULATIVE EFFECT OF TIME-DEPENDENT EXPOSURE AT EACH TIME POINTS
  ww=weight_mat[q,]
  
  WCE=matrix(NA,n,m)
  WCE[,1]=rep(0,n)
  WCE[,2]=dose[,1]*rev(ww[1])
  for (u in 3:m){
    WCE[,u]=dose[,1:u-1] %*% rev(ww[1:u-1])
  }
  
  # EXPOSURE MATRICES
  Xmat=matrix(ncol=2, nrow=n*m)
  Xmat <- cbind(rep(rnorm(n,5,1), each=m), as.vector(t(WCE)))   # binary fixed baseline covariates
  
  
  
  #########################################################################################################################  
  
  for (i in 1:B){
    
    # SET THE SEED
    seed <- 1000 + i
    set.seed(seed)
    
    #### GENERATE COMPETING RISKS DATA
    
    
    v<- sapply(1:2, function(i){runif(n)})
    failure <-2-rbinom(n,1,prob1)                                                            # Probability of failure on cause 1 
    eventRandom<-rep(NA,n)
    eventRandom[failure==1]=round((-log(1-(1/prob1)*v[failure==1,1]*prob1))*200)+1            # Latent main event time
    failure[is.na(eventRandom)] <- 2
    eventRandom[failure==2] <- round(-log(1-v[failure==2,2])*200)+1                           # Latent competing event time
    censorRandom <- round(runif(n,1,m*2),0)                                                  # censoring time                                      
    
    
    # SAMPLE THE SURVIVAL TIMES AND EVENTS: PERMUTATION ALGORITHM 
    simdata <- PAlgo_SH(n, m, Xmat, XmatNames=c( "cov1","WCE"),
                        eventRandom = eventRandom, censorRandom=censorRandom, failureCause , betas1=c(log(beta1[1]),log(beta1[2])), betas2=c(log(1.1),log(1.1)))     
    
    
    # MERGE WITH ORIGINAL DATA
    simdata2=data.frame(cbind(id=as.vector(rep(seq(1:n),each=m)),stop=as.vector(rep(seq(1:m),n)),dose=as.vector(t(dose))))
    new=merge(simdata,simdata2)
    simdata=new[order(new$id,new$stop),]
    
    
    #########################################################################################################################  
    # FIT THE MODEL 
    
    fit=mclapply(1:3,function(k) crr_WCE(simdata, lag=lag, constrained = F,  id = "id", event = "event", start = "start", stop = "stop", fup="fup",adm.cens="adm.cens.exit",expos = "dose", covariates ="cov1",censor=censor_list[k]))
    
    
    # SAVE RESULTS
    for (k in 1:3){
      if (is.null(fit[[k]])==TRUE){est.weight_ac[[k]][i,]=CP[[k]][q,i]=edf[[k]][q,i]=NA} else{
        
        # COVERAGE PROBABILITY
        CP[[k]][q,i]=sum(fit[[k]]$weight-weight_mat[q,][1:lag]<=1.96*fit[[k]]$se.weight & fit[[k]]$weight-weight_mat[q,][1:lag]>=-1.96*fit[[k]]$se.weight)/(lag+1)   
        
        # EDF
        edf[[k]][q,i]=fit[[k]]$edf
        
        # ESTIMATED WEIGHT
        est.weight_ac[[q]][i,]=fit[[1]]$weight
        est.weight_ipcw[[q]][i,]=fit[[2]]$weight
        est.weight_mi[[q]][i,]=fit[[3]]$weight
        
        # STANDARD ERROR OF WEIGHT
        #se.weight[[k]][i,]=fit[[k]]$se.weight
        
      }
    }
  }
}


proc.time()-time0



## OUTPUT RESULTS AND PLOT WEIGHT FUNCTION

lapply(edf,function(x) apply(x,1,mean))
lapply(CP,function(x) apply(x,1,mean))
lapply(CP,function(x) apply(x,1,median))


par(mfrow=c(2,2),mar=c(5,5,3,2))
for (q in 1:4){
  plot(weight_mat[q,1:lag],type="l",lwd=3,ylim=c(-0.02,0.04),xlab="Days before last follow-up",ylab="Estimated weight function")
  for (b in 1:B){lines(est.weight_ac[[q]][b,],col="grey")}
  for (b in 1:B){lines(est.weight_ipcw[[q]][b,],col="grey")}
  for (b in 1:B){lines(est.weight_mi[[q]][b,],col="grey")}
}
