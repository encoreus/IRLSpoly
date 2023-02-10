# Code of 'IRLS for Polyserial and Polychoric Correlation Coefficients'

# library(devtools)
# install_github('encoreus/IRLSpoly')
library(IRLSpoly)
library(polycor)
library(psych)

################################################################################
# Code for Figure 1
# The track of the estimating process of the tetrachoric correlation coefficient
set.seed(1997)
rho<-0.5
N <- 100
a = 0
b = 0
# Generate simulated tetrachoric pairs
X = gen_rho(N,rho)
x1 = X[1,]
x2 = X[2,]
X = rbind(ifelse(x1>=0,2,1),ifelse(x2>=0,2,1))
r0<-cor(X[1,],X[2,])
res = esti_polychoric(X)

EX <-res$Ex
EY = res$Ey
iter = res$iter

# The regression line of the first iteration
lm0 <- lm(c(EY[1,1],EY[2,1])~0 + c(EX[1,1],EX[2,1]))
# The regression line of the last iteration
lm1 <- lm(c(EY[1,iter],EY[2,iter])~0 + c(EX[1,iter],EX[2,iter]))
# Draw the simulated points
plot(x1,x2,pch=19,xlim=c(-1,1),ylim=c(-1,1),xlab=expression(Z[1]),ylab=expression(Z[2]),col='gray')
# Points the track.
points(cbind(EX[1,1:iter],EY[1,1:iter]),pch='*',col='red',type='b',lty=3)
points(cbind(EX[2,1:iter],EY[2,1:iter]),pch=3,col='blue',type='b',lty=3)
abline(a=0,b=r0,col='red')
abline(a=0,b=rho,col='green')
abline(lm1)
legend("topleft", inset=0.03, text.width=0.5, x.intersp=0.3, c("pearson correlation","true correlation", "tetrachoric correlation"),
       lty=c(1, 1, 1), cex=0.8,col=c("red", "green", 'black'))
grid()

################################################################################
# Code for Table 1
# Write a function to repeatedly sample and estimate N times for specific threshold a and rho.
# Use both IRLS and MLE(package polycor) methods
sim_polyserial= function(N,n,rho,a){
  rho0 = rep(0,N)
  std0 = rep(0,N)
  rho1 = rep(0,N)
  std1 = rep(0,N)
  for (i in 1:N) {
    X = gen_polyseries(n,rho,a)
    out0 = esti_polyserial(X)
    out1 = polyserial(X[1,],X[2,],ML=T,std.err=T)
    rho0[i] = out0$rho
    std0[i] = out0$std
    rho1[i] = out1$rho
    std1[i] = sqrt(out1$var[1,1])
    if(i%%100==0){print(i)}
  }
  return(list('est_rho'=rho0,
              'ML_rho'=rho1,
              'est_std'=std0,
              'ML_std'=std1))
}
# In this paper, setting:
# N = 30, 50, 100, 500, 1000
# rho = 0, 0.2, 0.4, 0.6, 0.8
# s1 = 2, 3, 5, 7
# Each category has an equal probability to fall in
# Change the sample sizes N, true coefficient rho and the number of categories to get different parts in Tabel 1
# With N and s1 increasing, the time cost of MLE will be huge.
set.seed(1998)
n = 1000
N = 50
rho = 0.2
s1 = 2
a1 = qnorm(1:s1/s1)[1:(s1-1)]
out1 = sim_polyserial(n,N,rho,a1)

# Based on the 1000 times MC simulation with IRLS and MLE, we calculate:
# mean of the 1000 rhohats
# standard deviation of the 1000 rhohats
# mean bias of the 1000 rhohats
# rmse bias of the 1000 rhohats
# mean of the 1000 estimated standard deviation

re = c(mean(out1$est_rho),
       sqrt(var(out1$est_rho)),
       mb(out1$est_rho,rho),
       rmse(out1$est_rho,rho),
       mean(out1$est_std),
       mean(out1$ML_rho),
       sqrt(var(out1$ML_rho)),
       mb(out1$ML_rho,rho),
       rmse(out1$ML_rho,rho))

# To get the whole Table1, you should traverse all N,rho,s1. The code followed will cost a long time so I annotate it.

# out_polyserial = list()
# for (N in c(30, 50, 100, 500, 1000)) {
#   for (rho in c(0, 0.2, 0.4, 0.6, 0.8)) {
#     for (s1 in c(2, 3, 5, 7)) {
#       out1 = sim_polyserial(n,N,rho,a1)
#       re = c(mean(out1$est_rho),
#              sqrt(var(out1$est_rho)),
#              mb(out1$est_rho,rho),
#              rmse(out1$est_rho,rho),
#              mean(out1$est_std),
#              mean(out1$ML_rho),
#              sqrt(var(out1$ML_rho)),
#              mb(out1$ML_rho,rho),
#              rmse(out1$ML_rho,rho))
#       tag1 = paste(as.character(N),as.character(rho),as.character(s1))
#       out_polyserial[[tag1]] = re
#     }
#   }
# }

################################################################################
# Code for Table 2
# Write a function to repeatedly sample and estimate N times for specific threshold a and rho.
# Use both IRLS and MLE(package polycor) methods
sim_polychoric= function(N,n,rho,a,b){
  rho0 = rep(0,N)
  std0 = rep(0,N)
  rho1 = rep(0,N)
  std1 = rep(0,N)
  for (i in 1:N) {
    X = gen_polychoric(n,rho,a,b)
    out0 = try(esti_polychoric(X),T)
    if(typeof(out0)=='character'){
      rho0[i] = NA
      std0[i] = NA
    }else{
      rho0[i] = out0$rho
      std0[i] = out0$std
    }
    out1 = try(polychor(X[1,],X[2,],ML=T,std.err=T),T)
    if(typeof(out1)=='character'){
      rho1[i] = NA
      std1[i] = NA
    }else{
      rho1[i] = out1$rho
      std1[i] = sqrt(out1$var[1,1])
    }
    if(i%%100==0){print(i)}
  }
  return(list('est_rho'=rho0,
              'ML_rho'=rho1,
              'est_std'=std0,
              'ML_std'=std1))
}
# In this paper, setting:
# N = 30, 50, 100, 500, 1000
# rho = 0, 0.2, 0.4, 0.6, 0.8
# s1 = 2, 3, 5, 7
# Each category has an equal probability to fall in.
# Change the sample sizes N, true coefficient rho and the number of categories to get different parts in Tabel 2.
# With N and s1 increasing, the time cost of MLE will be huge.
# You can change the parameter of function 'polychor' in line 141  to 'F' to accelerate the estimation.


set.seed(1998)
n = 1000
N = 50
rho = 0.2
s1 = 2
a1 = qnorm(1:s1/s1)[1:(s1-1)]
out1 = sim_polychoric(n,N,rho,a1,a1)

re = c(mean(out1$est_rho),
       sqrt(var(out1$est_rho)),
       mb(out1$est_rho,rho),
       rmse(out1$est_rho,rho),
       mean(out1$est_std),
       mean(out1$ML_rho),
       sqrt(var(out1$ML_rho)),
       mb(out1$ML_rho,rho),
       rmse(out1$ML_rho,rho))

# To get the whole Table 2, you should traverse all N,rho,s1. The code followed will cost a long time so I annotate it.

# out_polychoric = list()
# for (N in c(30, 50, 100, 500, 1000)) {
#   for (rho in c(0, 0.2, 0.4, 0.6, 0.8)) {
#     for (s1 in c(2, 3, 5, 7)) {
#       out1 = sim_polychoric(n,N,rho,a1,a1)
#       re = c(mean(out1$est_rho),
#              sqrt(var(out1$est_rho)),
#              mb(out1$est_rho,rho),
#              rmse(out1$est_rho,rho),
#              mean(out1$est_std),
#              mean(out1$ML_rho),
#              sqrt(var(out1$ML_rho)),
#              mb(out1$ML_rho,rho),
#              rmse(out1$ML_rho,rho))
#       tag1 = paste(as.character(N),as.character(rho),as.character(s1))
#       out_polychoric[[tag1]] = re
#     }
#   }
# }

################################################################################
# Code for Table 3
# Comparison of running time of three methods
# IRLS, MLE, 2-step MLE
# Both polyserial and polychoric
# True coefficient rho=0.4
# The number of categories s1=2,3,5,7
# Sample size N=500
# Replicated times n=1000

set.seed(2001)
N = 500
n = 1000
rho = 0.4

## Running time of polyserial
s1 = 2
a1 = qnorm(1:s1/s1)[1:(s1-1)]
tt1 = 0
tt2 = 0
tt3 = 0
for (i in 1:n) {
  X = gen_polyseries(N,rho,a1)

  ptm <- proc.time()
  out1 = try(esti_polyserial(X))
  ptm = proc.time()-ptm
  tt1 = tt1 + ptm

  ptm <- proc.time()
  out2 = try(polyserial(X[1,],X[2,],ML=T,std.err=T))
  ptm = proc.time()-ptm
  tt2 = tt2 + ptm

  ptm <- proc.time()
  out3 = try(polyserial(X[1,],X[2,],ML=F,std.err=T))
  ptm = proc.time()-ptm
  tt3 = tt3 + ptm
  if(i%%100==0){print(i)}
}
tt1
tt2
tt3

## Running time of polychoric
s1 = 2
a1 = qnorm(1:s1/s1)[1:(s1-1)]
tt1 = 0
tt2 = 0
tt3 = 0
for (i in 1:n) {
  X = gen_polychoric(N,rho,a1,a1)

  ptm <- proc.time()
  out1 = try(esti_polychoric(X))
  ptm = proc.time()-ptm
  tt1 = tt1 + ptm

  ptm <- proc.time()
  out2 = try(polychor(X[1,],X[2,],ML=T,std.err=F))
  ptm = proc.time()-ptm
  tt2 = tt2 + ptm

  ptm <- proc.time()
  out3 = try(polychor(X[1,],X[2,],ML=F,std.err=F))
  ptm = proc.time()-ptm
  tt3 = tt3 + ptm
  if(i%%100==0){print(i)}
}
tt1
tt2
tt3

################################################################################
# Code for Figure 3
# SPLOM for part of bif dataset in package psych
pairs_panels1(bfi[,1:5])
pairs_panels1(bfi[,1:5], MLE = TRUE)

################################################################################
# Code for Table 4
# The running time of calculate full bfi dataset in 3 methods.
# This will take a long time.
l1 = 27
t1 = 0
t2 = 0
t3 = 0
for (i in 2:l1) {
  for (j in 1:(i-1)) {
    df = na.omit(bfi[c(i,j)])
    df = t(as.matrix(df))

    ptm1 <- proc.time()
    test3 = try(esti_polychoric(df))
    ptm1 = proc.time()-ptm1
    t1 = t1 + ptm1

    ptm2 <- proc.time()
    test3 = try(polychor(df[1,],df[2,],ML=T,std.err = T))
    ptm2 = proc.time()-ptm2
    t2 = t2 + ptm2

    ptm3 <- proc.time()
    test3 = try(polychor(df[1,],df[2,],ML=F,std.err = T))
    ptm3 = proc.time()-ptm3
    t3 = t3 + ptm3
  }
  print(i)
}
t1
t2
t3

################################################################################
# Code for Figure 4
# SPLOM for part of built-in dataset Parenteral Nutrition
names1 = c('clinical_stages','dietary_status','NRS-2002','nausea','weakness')
pairs_panels1(Parenteral_nutrition[names1])
pairs_panels1(Parenteral_nutrition[names1],MLE = TRUE)

################################################################################
# Code for Table 5
# The running time of calculate full Parenteral Nutrition dataset in 3 methods.
# This will take a long time.
l0 = c(-1,-4,-9,-13)
PN = Parenteral_nutrition[l0]
l1 = ncol(PN)
t1 = 0
t2 = 0
t3 = 0
for (i in 2:l1) {
  for (j in 1:(i-1)) {
    df = na.omit(PN[c(i,j)])
    df = t(as.matrix(df))

    ptm1 <- proc.time()
    test3 = try(esti_polychoric(df))
    ptm1 = proc.time()-ptm1
    t1 = t1 + ptm1

    ptm2 <- proc.time()
    test3 = try(polychor(df[1,],df[2,],ML=T,std.err = T))
    ptm2 = proc.time()-ptm2
    t2 = t2 + ptm2

    ptm3 <- proc.time()
    test3 = try(polychor(df[1,],df[2,],ML=F,std.err = T))
    ptm3 = proc.time()-ptm3
    t3 = t3 + ptm3
  }
  print(i)
}
t1
t2
t3
