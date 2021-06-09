##############################################
### Internal function for data analysis   ####
##############################################
library(survival)
library(kmi)
library(mitools)
library(plyr)
library(splines)
library(mgcv)
library(zoo)
library(dplyr)         


########################################################################################################################
# FUNCTION TO CALCULATE WCE BASIS (D) FOR A GIVEN INDIVIDUAL FOR ALL RELEVANT RISK SETS 

Wcecalc <- function(ev, dose, stop, pbasis, lag){
  fup <- length(dose)
  myev <- ev[ev<=stop[fup]] 
  if (length(myev)>0)    {
    linesfor1 <- matrix(NA, ncol=dim(pbasis)[2], nrow=length(myev))                   
    for (i in 1:length(myev)){
      vec <- dose[stop <= myev[i]] 
      pos <- length(vec)
      if (pos<lag) {
        vec <- c(rep(0, (lag-length(vec))), vec)} else {
          pos <- length(vec)
          vec <- vec[(pos-lag+1):pos]} 
      linesfor1[i,] <- rev(vec)%*%pbasis}
  }   else {linesfor1 <- rep(NA, dim(pbasis)[2])}
  linesfor1}



########################################################################################################################
# FUNCTION FOR SUBDISTRIBUTION HAZARDS REGRESSION MODEL WITH CUMULATIVE WEIGHTED EXPOSURE 


crr_WCE <- function(data, lag, constrained,  id, event , start, stop, fup, adm.cens, expos, covariates = NULL, censor) {
  
  maxTime <- max(data[,stop])+1
  if (lag > maxTime)  stop("ERROR: lag must be smaller than the longest follow-up time")  
  coef.alpha = var.alpha = aic = penalty = cumhaz.base = timehaz = surv.base = covbeta = covse =list()
  WCEmat=var.WCEmat=coef.weight=var.weight=rep(NA,lag)
  if (is.null(covariates) == F){covbeta = covrse <- rep(0, length(covariates))}         
  
  n.events <-   sum(data[,event]==1)
  unique_fup <- sort(unique(data[,fup]))
  
  
  ## Using smoothCon function in mgcv package to construct p-spline basis matrix and penalty matrix
  
  temp=seq(1:lag)
  pb=smoothCon(s(temp,bs="ps",k=12),data=data.frame(temp))
  pbasis=pb[[1]]$X      # model design matrix
  S=pb[[1]]$S[[1]]      # penalty matrix
  
  ##################### constrained model ###############
  
  
  if (constrained==T){
    S[nrow(S),ncol(S)]=0.1          # additional penalty on the last coefficient
    
    #pbasis=pbasis[,1:(ncol(pbasis)-2)]    # set last two coefficient=0
    #S=S[1:(ncol(pbasis)),1:(ncol(pbasis))] 
    
  }   
  
  penalty[[1]]=S
  uid <- unique(data[,id])
  
  
  # Calculate new time-dependent covariates Dvar (B %*% X)
  kal <- data.frame(do.call("rbind", lapply(1:length(uid), function(i) Wcecalc(data[data[,id]==uid[i],stop], data[data[,id]==uid[i],expos],data[data[,id]==uid[i],stop],pbasis, lag))))
  kal <- kal[is.na(kal[,1])==FALSE,]
  names(kal) <- paste("Dvar", 1:dim(kal)[2], sep="")   
  kal=as.matrix(kal)
  data=cbind(data,kal)
  
  # Fit the model
  
  ## administrative censored data
  
  if (censor=="AC"){
    data$event2=as.numeric(data[,event]==1)
    formula=as.formula(paste("event2~", paste(c("s(log(adm.cens.exit),k=5)", "kal" ,covariates ),collapse="+")))
    co=gam(formula,data,family=binomial(link=cloglog),paraPen=list(kal=penalty),method="REML",random= ~1|id)
    coef = co$coefficients
    vcov = co$Vp
    
    # Output and save estimated parameters Phi, gamma, aic, effective degree of freedom
    
    coef.Phi = coef[grep("Dvar",names(co$coefficients))]
    vcov.Phi = vcov[grep("Dvar",names(co$coefficients)),grep("Dvar",names(co$coefficients))]
    
    if (is.null(covariates) == F){
      covbeta = coef[c(covariates)]
      covse = vcov[names(co$coefficients) %in% covariates,names(co$coefficients) %in% covariates] }
    
    aic=co$aic
    edf=sum(co$edf[grep("Dvar",names(co$coefficients))])
  }
  
  ## right censored data: IPCW
  
  if (censor=="IPCW")  {
    
    
    # expand dataset and calculate weights after competing event
    max=max(data[,fup])
    data_cr=data[data[,event]==2,]
    data_cr$Fup_new[data_cr[,event]==2]=max
    data_cr_impute <- survSplit(Surv(Fup_new, event) ~., data_cr,
                                cut=seq(1:max), start="start",end="stop")
    data_cr_impute[,event]=0
    #data_cr_impute <- subset(data_cr_impute, data_cr_impute[,stop] > data_cr_impute[,fup])
    data_cr_impute <- subset(data_cr_impute, data_cr_impute$stop > data_cr_impute$fup)
    
    
    
    # IPCW weight
    data_last=ddply(data,.(id), tail,1) 
    km <- survfit(Surv(fup, event == 0)~1, data=data_last)          ## Kaplan-Meier estimate of G
    #cox <- coxph(Surv(fup, event == 0)~Cov1, data=data_last)       ## Cox proportional hazard model estimate of G
    #survfit(cox)$surv
    
    G=data.frame(time=km$time,w=km$surv)
    G=merge(G,data.frame(time=seq(1:max)),all=T)
    
    w=na.locf(G$w)
    if (length(w)<length(G$w)) {w=c(rep(1,(length(G$w)-length(w))),w) } 
    data_cr_impute$weight=rep(NA,nrow(data_cr_impute))
    
    wmat=w %*% t(1/w)
    for (j in 1:length(unique(data_cr[,id]))){
      if (data_cr[j,fup]< max){
        data_cr_impute$weight[data_cr_impute[,id]==data_cr[j,id]]=1/wmat[data_cr[j,fup],((data_cr[j,fup]+1):max)]
      } }
    data_cr_impute=data_cr_impute[!(data_cr_impute$weight==0),]
    
    data_ipcw=rbind(cbind(data,weight=rep(1,nrow(data))),data_cr_impute)
    data_ipcw=data_ipcw[order(data_ipcw[,id],data_ipcw[,start]),]
    data_ipcw[data_ipcw[,event]==2,event]=0
    
    kal_ipcw=as.matrix(data_ipcw[grep("Dvar",names(data_ipcw))])
    
    formula=as.formula(paste("event~", paste(c("s(log(stop),k=5)", "kal_ipcw" ,covariates ),collapse="+")))
    co=gam(formula,data_ipcw,weight=data$weight,family=binomial(link=cloglog),paraPen=list(kal_ipcw=penalty),method="REML",random = ~1|id)
    
    coef = co$coefficients
    vcov = co$Vp   # Bayesian posteriorvcovariance matrix
    
    # Output and save estimated parameters Phi, gamma, aic, effective degree of freedom
    
    coef.Phi = coef[grep("Dvar",names(co$coefficients))]
    vcov.Phi = vcov[grep("Dvar",names(co$coefficients)),grep("Dvar",names(co$coefficients))]
    
    if (is.null(covariates) == F){
      covbeta = coef[c(covariates)]
      covse = vcov[names(co$coefficients) %in% covariates,names(co$coefficients) %in% covariates] }
    
    aic=co$aic
    edf=sum(co$edf[grep("Dvar",names(co$coefficients))])
  }
  
  ## right censored data: Multiple imputation
  
  if (censor=="MI"){
    imp.data <-  tryCatch(kmi(Surv(start, stop, event != 0) ~ 1, data = data,                                           # KM estimate for G
                              etype = event, id = id, failcode = 1, nimp = 5), error=function(e) NULL) 
    
    #imp.data <-  tryCatch(kmi(Surv(start, stop, event != 0) ~ age+male+race1+race2+race3+charlson, data = data,        # Cox model estimate for G
    #                          etype = event, id = id, failcode = 1, nimp = 5), error=function(e) NULL) 
    
    if (is.null(imp.data)==FALSE){
      
      info <- imp.data$info
      formula=as.formula(paste("event2~", paste(c(info[1], "kal",covariates),collapse="+")))
      result <- lapply(seq_along(imp.data$imputed.data), function(i) {
        datan <- imp.data$original.data
        datan[, info[1]] <- imp.data$imputed.data[[i]][, 1]
        datan$event2=as.numeric(datan$event==1)
        datan=datan[order(datan$id,datan$stop),]
        formula=as.formula(paste("event2~", paste(c("s(log(stop),k=5)", "kal" ,covariates ),collapse="+")))
        co=gam(formula,datan,family=binomial(link=cloglog),paraPen=list(kal=penalty),method="REML",random=~1|id)
        co})
      
      # Output and save estimated parameters Phi, gamma, aic, effective degree of freedom  
      
      coef = MIcombine(result)$coefficients
      vcov = MIcombine(result)$variance
      
      coef.Phi = coef[grep("Dvar",names(MIcombine(result)$coefficients))]     ## average coefficient estimate from multiple imputed datasets
      vcov.Phi = vcov[grep("Dvar",names(MIcombine(result)$coefficients)),grep("Dvar",names(MIcombine(result)$coefficients))]  ## combine from Rubin's rule
      
      if (is.null(covariates) == F){
        covbeta = coef[c(covariates)]
        covse =  vcov[names(coef) %in% covariates,names(coef) %in% covariates] }
      
      aic= mean(unlist(lapply(seq_along(imp.data$imputed.data), function(i) {result[[i]]$aic})))
      edf=mean(unlist(lapply(seq_along(imp.data$imputed.data), function(i) {sum(result[[i]]$edf[grep("Dvar",names(result[[i]]$coefficients))])})))
    }else {coef=NULL}}
  
  
  
  if (is.null(coef)==F){ 
    
    coef.alpha <- as.numeric(t(coef.Phi) %*% apply(pbasis[,1:length(coef.Phi)],2,sum))        # estimated alpha 
    var.alpha <-t(apply(pbasis[,1:length(coef.Phi)],2,sum)) %*% vcov.Phi %*% apply(pbasis[,1:length(coef.Phi)],2,sum)   ## estimated variance of alpha
    
    
    weight <- pbasis[,1:length(coef.Phi)] %*% coef.Phi/coef.alpha         # estimated weight 
    # estimated variance-covariance matrix of weight and point-wise standard error
    h=(diag(length(coef.Phi)) *coef.alpha-as.matrix(coef.Phi) %*% apply(pbasis[,1:length(coef.Phi)],2,sum))/coef.alpha^2
    vcov.weight <- pbasis[,1:length(coef.Phi)] %*% t(h) %*% vcov.Phi %*% h %*% t(pbasis[,1:length(coef.Phi)]) 
    pointwise.se.weight <- sqrt(diag(vcov.weight))
    
    
    WCEmat=pbasis[,1:length(coef.Phi)] %*% coef.Phi    ## alpha*weight
    var.WCEmat=pbasis[,1:length(coef.Phi)] %*% vcov.Phi %*% t(pbasis[,1:length(coef.Phi)])
    
    # cumulative baseline subdistribution hazard and baseline CIF
    m=maxTime
    data_new=cbind(rep(1,m),seq(1:m), matrix(0,m,(length(coef)-2)))
    cumhaz.base=cumsum(1-exp(-exp(data_new %*% coef)))
    #cumhaz.base=cumsum(predict(co, data.frame(data_new),type="response"))
    surv.base=1-exp(-cumhaz.base)
    timehaz=seq(1:length(surv.base))
    
    
    if (is.null(covariates)){
      est <- list( WCEmat = WCEmat, var.WCEmat=var.WCEmat, est.Phi = coef.Phi, vcovmat.Phi = vcov.Phi, weight=weight, se.weight=pointwise.se.weight, timehaz=timehaz, surv.base=surv.base, est.alpha=coef.alpha, vcovmat.alpha=var.alpha, constrained = constrained, covariates = NULL, ne = n.events, aic = aic, edf=edf)} else {
        est <- list( WCEmat = WCEmat, var.WCEmat=var.WCEmat, est.Phi = coef.Phi, vcovmat.Phi = vcov.Phi, weight=weight, se.weight=pointwise.se.weight, timehaz=timehaz, surv.base=surv.base, est.alpha=coef.alpha, vcovmat.alpha=var.alpha, beta.covariates=covbeta, se.covariates= covse, covariates = covariates, constrained = constrained, ne = n.events, aic = aic, edf=edf)}
    return(est)}else return (est=NULL)
}
