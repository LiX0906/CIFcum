#################################################################################
####### Internal function for data generation  ##################################
##################################################################################


########################################################################################################################
# FUNCTION TO CALCULATE PARTIAL LIKELIHOOD

partialHazards <-
  function (t, v, covArray, betas) 
  {
    if (length(betas) > 1) {return(exp(covArray[v, t, ] %*% betas))} else {return(exp(covArray[v, t, ] * betas))}
  }


########################################################################################################################
# FUNCTION FOR PERMUTATION SAMPLING BASED ON PARTIAL LIKELIHOOD

PermuteCovariateVectors <-
  function (t, d, count, I, covArray, betas1, betas2) 
  {
    n = sum(count[I])
    p <- integer(n)
    v <- seq(n)
    ip = 1
    for (k in I) {
      if (d[k]==1) {
        J = sample(length(v), size = count[k], replace = FALSE, prob = partialHazards(t[k], v, covArray, betas1))
      }
      else if (d[k]==2){
        J = sample(length(v), size = count[k], replace = FALSE, prob = partialHazards(t[k], v, covArray, betas2))
      }
      else {
        J = sample(length(v), size = count[k], replace = FALSE, 
                   prob = NULL)
      }
      p[seq(from = ip, along.with = J)] <- v[J]
      ip <- ip + count[k]
      v = v[-J]
    }
    return(p)
  }


########################################################################################################################
# FORMAT GENERATED DATA

form.data <-
  function (ordered.info, m, n, Xmat, XmatNames = NULL, censorTime) 
  {
    data <- data.frame(cbind(sort(rep(seq(1, n), m)), 
                             rep(0, m * n), rep(ordered.info$obs.t, each = m), rep(seq(0, (m-1)), n), 
                             rep(seq(1, m), n), Xmat,rep(censorTime[ordered.info$id.tuples],each=m)))
    if (is.null(XmatNames) == TRUE) {
      XmatNames <- paste("X", seq(1, dim(Xmat)[2]), sep = "")
    }
    colnames(data) <- c("id", "event", "fup", "start", "stop", 
                        XmatNames,"adm.cens.exit")
    data$event[data$stop == data$fup & rep(ordered.info$d, each = m) == 
                 1] <- 1
    data$event[data$stop == data$fup & rep(ordered.info$d, each = m) == 
                 2] <- 2
    data <- subset(data, data$stop <= data$fup)
    data$adm.cens.exit[data$event !=2]=data$stop[data$event !=2]
    return(data)
  }



################################################################################
# FUNCTIONS TO GENERATE TIME-DEPENDENT EXPOSURE PROFILE (#sub=n,follow-up=m)

TDexp <- function(n,m){
  Dose=matrix(0,n,m)
  for(i in seq(n)){
    start <- round(runif(1,1,m),0) # individual start date
    duration <- 7 + 7*rpois(1,3) # in days
    dose <- round(runif(1,0,10),1)      # Uniform distributed dosage
    vec <- c(rep(0, start-1), rep(dose, duration))
    while (length(vec)<=m){
      intermission <- 14 + 7*rpois(1,3) # in days
      duration <- 7 + 7*rpois(1,3) # in days
      dose <- round(runif(1,0,10),1)      
      vec <- append(vec, c(rep(0, intermission), rep(dose, duration)))}
    Dose[i,]=vec[1:m]}
  rownames(Dose) <- seq(n)
  return(Dose)
}

#########################################################################################################################  
# FUNCTION TO SAMPLE THE SURVIVAL TIMES AND EVENTS: PERMUTATION ALGORITHM FOR COMPETING RISKS BASED ON SUBDISTRIBUTION HAZARDS MODEL
  
  PAlgo_SH <-
    function (numSubjects, maxTime, Xmat, XmatNames = NULL, eventRandom, 
              censorRandom, failureCause, betas1, betas2) 
    {
      if (length(betas1) ==1 ) {Xmat = matrix(Xmat, ncol=1)}
      nc <- dim(Xmat)[2]   #The number of covariates
      I <- rep(seq(numSubjects), each = maxTime, times = nc)
      J <- rep(seq(maxTime), times = nc * numSubjects)
      K <- rep(seq(nc), each = maxTime * numSubjects)
      covArray <- array(0, dim = c(numSubjects, maxTime, nc))
      covArray[cbind(I, J, K)] <- Xmat         ## re-organize covariates in array form
      
      survivalTime <- as.integer(eventRandom)
      failureCause <- as.integer(failure)      ## event type: 1=main event, 2=competing event
      censorTime <- as.integer(censorRandom)
      
      survivalTime_max <-apply(cbind(survivalTime, maxTime),     
                               MARGIN = 1, FUN = min)
      notMax <- ifelse(survivalTime < maxTime, 1, 0)
      failureCause_max <-ifelse(notMax,failureCause,0)
        
      observedTime <- apply(cbind(survivalTime, censorTime, maxTime),     
                            MARGIN = 1, FUN = min)
      notCensored <- ifelse(survivalTime <= apply(cbind(censorTime,       
                                                        maxTime), MARGIN = 1, FUN = min), 1, 0)  # censor at censor time or maxtime
      failureCause <- ifelse(notCensored,failureCause,0)
      
      I <- count <- integer(0)
      
      
      # Match by permutation algorithm
      
      I <-c(order(survivalTime_max)[order(survivalTime_max) %in% which(failureCause_max==1)],
            order(survivalTime_max)[order(survivalTime_max) %in% which(failureCause_max==2)],
            order(survivalTime_max)[order(survivalTime_max) %in% which(failureCause_max==0)])  # Rank all event times with failure=1 in ascending order, followed by all event time with failure=2 in ascending order, followed by all event time at maxtime (censored)
      count = rep(1, times = length(I))
      p = PermuteCovariateVectors(survivalTime_max, failureCause_max, count, I, covArray, betas1, betas2) 
      
      
      tuples = cbind(obs.t = observedTime[I], d = failureCause[I], id.tuples = I)
      
      info = data.frame(cbind(tuples, cov.id = p))
      ordered.info = info[order(info$cov.id), ]
      return(form.data(ordered.info, maxTime, numSubjects, Xmat, 
                       XmatNames,censorTime))
    }
  
  
