#####################################
## Subject Accrual Suspension Rules
## Aim: Simulate suspension rule and 
##   sample size
#####################################

## Reference:

######## Part1: Derive Suspension Boundary #######
## Suspension Rule: P[p_tox|Data]>0.95

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
write.csv(sus.bds,'./Outputs/Table_bounds.csv')


#### Part2: Evaluate Sample Size via Simulation ####

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
write.csv(out,'./Outputs/Table_sim.csv')





