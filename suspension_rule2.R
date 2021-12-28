#####################################
## Subject Accrual Suspension Rules
## April 03, 2019
## Aim: Simulate suspension rule and 
##   sample size
#####################################


#### Part1: Sample Size Simulation ####
nSims=100000
out<-NULL
bounds<-c(rep(100,5), rep(5,3), rep(6,2), rep(7,2), rep(8,3), rep(9, 3), rep(10,2), rep(11,3), rep(12,2), rep(13,3),rep(14, 3),
          rep(15,3), rep(16,2), rep(17,3),rep(18,3),rep(19,3),rep(20,3),rep(21,2),rep(22,3),rep(23,3),rep(24,3),rep(25,3),
          rep(26,3), rep(27,3),rep(28,3), rep(29,2),rep(30,3),rep(31,3),31)
## large bound value 100 is used in the beginning so that no cross bound can happen in the beginning

for (t in 1:6){
  p<-t/10
  n.t<-NULL
  sus.t<-NULL
  for (i in 1:nSims){
    tox<-rbinom(80, 1, p)  ## 1: toxicity observed
    tox.cum<-cumsum(tox)  ## count of toxicity observations till index of a vector
    n.i<- ifelse(length(which(tox.cum>=bounds))>0, min(which(tox.cum>=bounds)), 80)  
    ## sample size at which more toxcities observed than the bounds to suspend study
    sus.i<- ifelse(n.i<80, 1, 0)  ## suspend or not
    
    n.t<-c(n.t, n.i)
    sus.t<-c(sus.t, sus.i)
  }
  
  pr<-mean(sus.t)
  n_25=quantile(n.t, 0.25)
  n_50=quantile(n.t, 0.50)
  n_75=quantile(n.t, 0.75)
  out.i<-c(p, pr, n_25, n_50, n_75)
  out.i<-rbind(NULL, out.i)
  out<-rbind(out, out.i)
}

out
colnames(out)<-c("True Toxicity Rate","Suspension Prob.","N, 25% Percentile","N, 50% Percentile", "N, 75% Percentile")
row.names(out)<-NULL
write.csv(out,'C:/Users/qzhao/OneDrive - Nektar Therapeutics/documents/Suspension Rules/Outputs/Table8.csv')



######## Part2: Replicate Suspension Rule #######
## replicate bounds list in Table 7
bds=NULL

n=1
bd.prev=1
while (n<=80){
  if (n<=5){
    bds<-c(bds, NA)
  }else{
    post.p<-0  ## need to reset post.p, otherwise while loop cannot continue
    bd=bd.prev-1
    while (post.p<0.95 & bd<n){
      bd=bd+1   ## put bd = bd+1 before post.p to make sure that the final bd will not violate post.p rule
      post.p<-1-pbeta(0.3, 0.6+bd, n-bd+1.4)
    }
    bd.prev<-bd  ## record bd corresponding to n=n-1
    bds<-c(bds, bd)
  }

  n<-n+1
}
bds
sus.bds<-data.frame(cbind(c(1:80), bds))
colnames(sus.bds)<-c("N Treated","Suspension Boundaries")
row.names(sus.bds)<-NULL
write.csv(sus.bds,'C:/Users/qzhao/OneDrive - Nektar Therapeutics/documents/Suspension Rules/Outputs/Table7.csv')



## Wrong Approach:
bds=NULL

n=1
bd = bd.prev=1
while (n<=80){
  if (n<=5){
    bds<-c(bds, NA)
  }else{
    post.p<-0
    #bd=bd.prev-1
    while (post.p<0.95 & bd<n){
      
      post.p<-1-pbeta(0.3, 0.6+bd, n-bd+1.4)
      bd=bd+1
    }
    #bd.prev<-bd  ## record bd corresponding to n=n-1
    bds<-c(bds, bd)
  }
  
  n<-n+1
}


############################################
## Evaluate Suspension Rules via Simulations
############################################

## For one cohort: p is a single value
out.fn <- function(tox, p, bounds){
  
  n <- length(tox)
  tox.cum<-cumsum(tox)  ## count of toxicity observations till index of a vector
  n.i<- ifelse(length(which(tox.cum>=bounds[1:n]))>0, min(which(tox.cum>=bounds[1:n])), n)  
  ## sample size at which more toxcities observed than the bounds to suspend study
  sus.i<- ifelse(n.i<n, 1, 0)  ## suspend or not
  
  fls.neg <- !sus.i & p>0.3
  fls.pos <- sus.i & p<=0.3
  n.exp <- ifelse(p>0.3, n.i, 0)
  
  out <- c(fls.neg, fls.pos, n.exp)
  return(out)
}

## For mulitple cohorts: p.v is a single vector
out.fn2 <- function(tox, p.v, bounds){
  
  n <- length(tox)
  n.coh <- length(p.v) 
  tox.cum<-cumsum(tox)  ## count of toxicity observations till index of a vector
  n.i<- ifelse(length(which(tox.cum>=bounds[1:n]))>0, min(which(tox.cum>=bounds[1:n])), n)  
  ## sample size at which more toxcities observed than the bounds to suspend study
  sus.i<- ifelse(n.i< n, 1, 0)  ## suspend or not
  
  fls.neg <- sum(!sus.i & p.v>0.3)
  fls.pos <- sum(sus.i & p.v<=0.3)
  
  rem <- (n.i %% n.coh)
  n.exp <- sum((n.i %/% n.coh)*(p.v>0.3)) + ifelse(rem==0, 0, sum((n.i %% n.coh)*(p.v[1:rem]>0.3)))
  
  out <- c(fls.neg, fls.pos, n.exp)
  return(out)
  
}

nSims <- 10000
p.v <- c(rep(0.2, 2), rep(0.5,3))
res1 <- NULL
res.pooled <- NULL
for (i in c(1:nSims)){
  set.seed(i)
  mtx <- NULL
  out1 <- rep(0, 3)
  for (p in p.v) {
    x<-rbinom(20, 1, p)
    out <- out.fn(x, p, bounds)
    out1 <- out1 + out
    mtx <- rbind(mtx, x)
  }
  res1 <- rbind(res1, out1)
  out.pooled <- out.fn2(as.vector(mtx), p.v, bounds)
  res.pooled <- rbind(res.pooled, out.pooled)
}


colMeans(res1)
colMeans(res.pooled)
as.vector(mtx)






