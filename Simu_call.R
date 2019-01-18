rm(list=ls())
library(survival)

time0=proc.time()

source("Function_data generation.R")
source("Function_analysis.R")


################################################################################
# SPECIFY PARAMETERS

B=100                                          # Monte carlo replicates
lag=120                                       # Maximum effective time window
m=365                                         # maximum follow up days
n=500                                         # number of subjects
beta1=c(1.5,4)
prob1=0.8                                    # probability of type 1 failure

################################################################################
# CREATE LISTS AND VECTORS TO SAVE RESULTS
MSE_ac=MSE_ipcw=MSE_mi=
  edf_ac=edf_ipcw=edf_mi=
    aic_ac=aic_ipcw=aic_mi=rep(NA,B)
result_ac=result_ipcw=result_mi=matrix(NA,4,3)
est.weight_ac=est.weight_ipcw=est.weight_mi=matrix(NA,B,lag)
output_ac=output_ipcw=output_mi=list() 

################################################################################
# SIMULATE WEIGHT FUNCTION (STANDARDIZED)
elapsed=seq(0,m)
weight1=ifelse(elapsed<lag,1/lag,0)                                            # Constant: Constant weight 
weight2=(1-1/(1+exp(-(elapsed-70)/12)))/70                                     # S-shape decay
weight3<- ifelse(elapsed<lag,dexp(0:lag,0.035),0)                          # Exponential decay
weight4 <-ifelse(elapsed<lag,dnorm(0:lag,30,10)/2+ dnorm(0:lag,85,25)/1.82,0)  # shape 4, e.g.
weight_mat=rbind(weight1,weight2,weight3,weight4)


for (q in 1:4){
  
######################################################################################################################### 
  # GENERATE EXPOSURE TIME-DEPENDENT EXPOSURE MATRIX
  dose <- TDexp(n,m)
  
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
  seed <- 10000 + i
  set.seed(seed)
    
  #### GENERATE COMPETING RISKS DATA
  
  eventRandom<-rep(NA,n)
  v<- sapply(1:2, function(i){runif(n)})
  failure <-2-rbinom(n,1,prob1)                                                            # Probability of failure on cause 1 
  eventRandom[failure==1]=round((-log(1-(1/prob1)*v[failure==1,1]*prob1))*200)+1            # Latent main event time
  failure[is.na(eventRandom)] <- 2
  eventRandom[failure==2] <- round(-log(1-v[failure==2,2])*200)+1                           # Latent competing event time
  censorRandom <- round(runif(n,1,m*2),0)                                                  # censoring time                                      
  
  
  # SAMPLE THE SURVIVAL TIMES AND EVENTS: PERMUTATION ALGORITHM 
  simdata <- PAlgo_SH(n, m, Xmat, XmatNames=c( "cov1","WCE"),
                   eventRandom = eventRandom, censorRandom=censorRandom, failureCause =failure, betas1=c(log(beta1[1]),log(beta1[2])), betas2=c(log(1.1),log(1.1)))     


  # MERGE WITH ORIGINAL DATA
  simdata2=data.frame(cbind(id=as.vector(rep(seq(1:n),each=m)),stop=as.vector(rep(seq(1:m),n)),dose=as.vector(t(dose))))
  new=merge(simdata,simdata2)
  simdata=new[order(new$id,new$stop),]
  

#########################################################################################################################  
  # FIT THE MODEL 
 
  ## AC
  fit_ac=crr_WCE(simdata, lag=lag, constrained = T,  id = "id", event = "event", start = "start", stop = "stop", adm.cens="adm.cens.exit",expos = "dose", covariates ="cov1",censor="AC")  # ADMINISTRATIVE CENSORING

  ## IPCW
  fit_ipcw=crr_WCE(simdata,lag=lag, constrained = T,  id = "id", event = "event", start = "start", stop = "stop", expos = "dose", covariates ="cov1",censor="IPCW") # RIGHT CENSORING: Inverse probability censoring weight

  ## MI
  fit_mi=crr_WCE(simdata,lag=lag, constrained = T,  id = "id", event = "event", start = "start", stop = "stop", expos = "dose", covariates ="cov1",censor="MI") # RIGHT CENSORING: Multiple imputation
  
#########################################################################################################################  
  # SAVE RESULTS
  
  # AC
  # ESTIMATED WEIGHT
  est.weight_ac[i,]=fit_ac$weight
  
  # MSE
  MSE_ac[i]=sqrt(sum((fit_ac$WCEmat-log(beta1[2])*weight_mat[q,][1:lag])^2)/(lag+1))    # relative MSE: compare to a reference
  
  # EDF
  edf_ac[i]=fit_ac$edf
  
  # AIC
  aic_ac[i]=fit_ac$aic
 
  # IPCW 
  if (is.null(fit_ipcw)==TRUE){est.weight_ipcw[i,]=MSE_ipcw[i]=edf_ipcw[i]=aic_ipcw[i]=NA} else{
    #ESTIMATED WEIGHT
    est.weight_ipcw[i,]=fit_ipcw$weight

    # MSE
    MSE_ipcw[i]=sqrt(sum((fit_ipcw$WCEmat-log(beta1[2])*weight_mat[q,][1:lag])^2)/(lag+1))    # relative MSE: compare to a reference
    
    # EDF
    edf_ipcw[i]=fit_ipcw$edf
    
    # AIC
    aic_ipcw[i]=fit_ipcw$aic
    }

  # MI
  if (is.null(fit_mi)==TRUE){est.weight_mi[i,]=MSE_mi[i]=edf_mi[i]=aic_mi[i]=NA} else{
    #ESTIMATED WEIGHT
    est.weight_mi[i,]=fit_mi$weight
    
    # MSE
    MSE_mi[i]=sqrt(sum((fit_mi$WCEmat-log(beta1[2])*weight_mat[q,][1:lag])^2)/(lag+1))    # relative MSE: compare to a reference
    
    # EDF
    edf_mi[i]=fit_mi$edf
    
    # AIC
    aic_mi[i]=fit_mi$aic
  }

  
}

result_ac[q,]=c(mean(na.omit(aic_ac)),mean(na.omit(edf_ac)),mean(na.omit(MSE_ac)))   
result_ipcw[q,]=c(mean(na.omit(aic_ipcw)),mean(na.omit(edf_ipcw)),mean(na.omit(MSE_ipcw)))   
result_mi[q,]=c(mean(na.omit(aic_mi)),mean(na.omit(edf_mi)),mean(na.omit(MSE_mi)))   

output_ac[[q]]=est.weight_ac
output_ipcw[[q]]=est.weight_ipcw
output_mi[[q]]=est.weight_mi
}



proc.time()-time0


result_ac
result_ipcw
result_mi
colnames(result_ac)=colnames(result_ipcw)=colnames(result_mi)=c("AIC","EDF","AMSE")
rownames(result_ac)=rownames(result_ipcw)=rownames(result_mi)=c("weight1","weight2","weight3","weight4")

par(mfrow=c(2,2),mar=c(5,5,3,2))
for (q in 1:4){
plot(weight_mat[q,1:lag],type="l",lwd=3)
for (b in 1:B){lines(output_ac[[q]][b,],col="green")}
for (b in 1:B){lines(output_ipcw[[q]][b,],col="red")}
for (b in 1:B){lines(output_mi[[q]][b,],col="blue")}
}

saveRDS(result_ac,file="/Users/xingyuanli/Dropbox/simulation/Output/20190110/result_ac_constrained.Rdata")
saveRDS(result_ipcw,file="/Users/xingyuanli/Dropbox/simulation/Output/20190110/result_ipcw_constrained.Rdata")
saveRDS(result_mi,file="/Users/xingyuanli/Dropbox/simulation/Output/20190110/result_mi_constrained.Rdata")

saveRDS(output_ac,file="/Users/xingyuanli/Dropbox/simulation/Output/20190110/output_ac_constrained.Rdata")
saveRDS(output_ipcw,file="/Users/xingyuanli/Dropbox/simulation/Output/20190110/output_ipcw_constrained.Rdata")
saveRDS(output_mi,file="/Users/xingyuanli/Dropbox/simulation/Output/20190110/output_mi_constrained.Rdata")


test=readRDS("/Users/xingyuanli/Dropbox/simulation/Output/20190110/result_ac.Rdata")


