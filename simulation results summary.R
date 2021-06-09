setwd("/Users/lixingyuan/Dropbox/Paper 1 submission/Supplemental materials")


# unconstrained

#q1
q1=readRDS("edf_unconstrained_q1.RData")
lapply(q1, function(x) apply(x,1,na.omit(mean)))
test=q1$MI[1,]
mean(na.omit(test))

q1=readRDS("CP_unconstrained_q1.RData")
lapply(q1, function(x) apply(x,1,na.omit(mean)))
test=q1$MI[1,]
mean(na.omit(test))

q1=readRDS("CP_unconstrained_q1.RData")
lapply(q1, function(x) apply(x,1,na.omit(median)))
test=q1$MI[1,]
median(na.omit(test))

# q2
q1=readRDS("edf_unconstrained_q2.RData")
lapply(q1, function(x) apply(x,1,na.omit(mean)))
test=q1$MI[2,]
mean(na.omit(test))

q1=readRDS("CP_unconstrained_q2.RData")
lapply(q1, function(x) apply(x,1,na.omit(mean)))
test=q1$MI[2,]
mean(na.omit(test))

q1=readRDS("CP_unconstrained_q2.RData")
lapply(q1, function(x) apply(x,1,na.omit(median)))
test=q1$MI[2,]
median(na.omit(test))


# q3
q1=readRDS("edf_unconstrained_q3.RData")
lapply(q1, function(x) apply(x,1,na.omit(mean)))
test=q1$MI[3,]
mean(na.omit(test))

q1=readRDS("CP_unconstrained_q3.RData")
lapply(q1, function(x) apply(x,1,na.omit(mean)))
test=q1$MI[3,]
mean(na.omit(test))

q1=readRDS("CP_unconstrained_q3.RData")
lapply(q1, function(x) apply(x,1,na.omit(median)))
test=q1$MI[3,]
median(na.omit(test))



# q4
q1=readRDS("edf_unconstrained_q4.RData")
lapply(q1, function(x) apply(x,1,na.omit(mean)))
test=q1$MI[4,]
mean(na.omit(test))

q1=readRDS("CP_unconstrained_q4.RData")
lapply(q1, function(x) apply(x,1,na.omit(mean)))
test=q1$MI[4,]
mean(na.omit(test))

q1=readRDS("CP_unconstrained_q4.RData")
lapply(q1, function(x) apply(x,1,na.omit(median)))
test=q1$MI[4,]
median(na.omit(test))



## constrained


edf=readRDS("edf_constrained.RData")
lapply(edf, function(x) apply(x,1,na.omit(mean)))
test=edf$MI
apply(test,1,function(x){mean(na.omit(x))})

CP=readRDS("CP_constrained.RData")
lapply(CP, function(x) apply(x,1,na.omit(mean)))
test=CP$MI
apply(test,1,function(x){mean(na.omit(x))})


lapply(CP, function(x) apply(x,1,na.omit(median)))
apply(test,1,function(x){median(na.omit(x))})

