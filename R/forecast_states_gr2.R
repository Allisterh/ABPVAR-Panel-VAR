#Forecasting performance using training data from May 2007 to September 2019
#ABPVAR model with grouping 2
setwd("C:/Users/nchak/OneDrive/Desktop/R_codes")
source("function_states.R")
load("data_US_states_gr2.dat")

library(matrixcalc)
library(dplyr)
library(MASS)
library(Matrix)
library(glmnet)
library(clusterGeneration)
library(fBasics)
library(base)
library(SuppDists)
library(LaplacesDemon)
library(mvnfast)
library(snow)
library(scoringRules)
library("CVglasso")
index <- 3
horizon <- 6
nobs <- 150 - index
k1 <- 0
k2 <- 15
num <- 5
npop <- 5
horizon <-10
N <- 1000
p <- (3*k1)+k2
p1 <- k1+k2
ngr <- p1 + (p1*(p1-num))
vec <- permutation(k1,k2)$vec
vec1 <- permutation(k1,k2)$vec1
vec2<- permutation(k1,k2)$vec2
M <- matrix(c(1:npop))
alpha <- 1
beta <- 1
lambda <- rep(1,ngr)
Q <- 0.3*diag(npop)
Q1 <- 0.3*diag(p)
Q2 <- 0.3*diag(num)
df <- 1
theta <- rep(0.5,npop)
iteration <- 2000
burn <- 1000
#=======================================Initialize====================================
sigma <- genPositiveDefMat(p, covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),  alphad=1, eta=1, rangeVar=c(0.1,1), lambdaLow=1, ratioLambda=2)$Sigma
sim_mat <- genPositiveDefMat(npop, covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),  alphad=1, eta=1, rangeVar=c(0.1,1), lambdaLow=1, ratioLambda=2)$Sigma
U_mat <- genPositiveDefMat(num, covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),  alphad=1, eta=1, rangeVar=c(0.1,1), lambdaLow=1, ratioLambda=2)$Sigma

tau <- rgamma(ngr,1.5,1.5)


for(i in 1:npop){
  for(j in 1: p)
  {
    lf[[i]][j,] <- lf[[i]][j,] - mean(lf[[i]][j,])
   
  }
}
gamma4_list <- list()
s1_list <- list();s2_list <- list();s3_list <- list()
data <- list();forecast_true <- list();d_A22 <- list()

for(j in 1:npop){
  data[[j]] <- lf[[j]][,1:(nobs + 1)]
  forecast_true[[j]] <- lf[[j]][,(nobs + 1 + 1):(nobs + 1 + horizon)]
  b4 <- numeric((k2*k2))
  s_func <- S_matrix(data[[j]],nobs)
  s1_list[[j]] <- s_func$s1
  s2_list[[j]]<- s_func$s2
  s3_list[[j]] <- s_func$s3
  A22_func <- para_A22(s1_list[[j]],s2_list[[j]],sigma)
  gamma4_list[[j]] <- A22_func$gamma_A22
  d_A22[[j]] <- A22_func$d_A22
}#check


#===================================Gibbs sampling===================================

A_sim <- list()
A_final <- list()

for(j in 1:npop){
  A_sim[[j]] <- array(0, c(p1, p1, (iteration-burn)))
  A_final[[j]] <- matrix(0,p1,p1)
  
}
sigma_sim <- array(0, c(p, p, (iteration-burn)))
sigma_final <- matrix(0,p,p)
sim_mat_sim <- array(0, c(npop, npop, (iteration-burn)))
sim_mat_final <- matrix(0,npop,npop)

mat <- list()
for(i in 1:p1){
  mat[[i]] <- tau[i]*U_mat
}
tau2 <- tau[-c(1:p1)]
d_tau <- bdiag(bdiag(mat),diag(tau2))    

count <- 0

for(ind in 1:iteration){
  
  #========covariance matrix for prior of A for all pop==========
  cov_A_temp <- kronecker(sim_mat,d_tau)
  #===============permutation for cov_A================================
  perm <- NULL
  for(m in 1:npop)
  {
    t1 <- NULL
    for(i in 1:num){
      t2 <- NULL
      for(j in 1:p1){
        t2 <- c(t2,((m-1)*(p1^2))+((j-1)*num)+i)
      }
      t1 <- c(t1,t2)
    }
    t3 <- NULL
    for(i in 1: (p1*(p1-num))){
      t3 <- c(t3,((m-1)*(p1^2))+(p1*num)+i)
    }
    perm1 <- c(t1,t3)
    perm <- c(perm,perm1)
  }
  #=========================================================================
  cov_A <- cov_A_temp[perm,]
  cov_A <- cov_A[,perm]
vec_A22 <- list()
for(m in 1:npop)
{
  vec_A22[[m]] <- c((((m-1)*(k2^2))+1):(((m-1)*(k2^2))+(k2^2)))
}

sig_A22 <- list()
for(j in 1:npop)
{
  sig_A22[[j]] <- as.matrix(cov_A[vec_A22[[j]],])
  sig_A22[[j]] <- as.matrix(sig_A22[[j]][,vec_A22[[j]]])
}

b4_all <- NULL

for(j in 1:npop)
{
  cov_A22 <- gamma4_list[[j]] + solve((sig_A22[[j]]))
  b4 <- rmvn(1,solve((cov_A22))%*%d_A22[[j]],solve((cov_A22)))
  b4_all <- c(b4_all,b4)
}
b4 <- matrix(b4_all,nrow=k2*k2,ncol=npop)

#=================drawing tau==================
A_list <- list()
A <- matrix(0,nrow=k2^2,ncol=npop)
#=======================A matrix for all population=======================
for(j in 1:npop){
  A_list[[j]] <- matrix(b4[,j],nrow=k2,ncol=k2)
  A[,j] <- vec(A_list[[j]]) #dim of A: (k1+k2)^2 by npop
}
A_tau <- matrix(0,nrow=p1^2,ncol=npop)
for(j in 1:npop){
  sub1 <- A_list[[j]][,1:num]
  sub2 <- A_list[[j]][,-c(1:num)]
  mm <- NULL
  for(i in 1:p1){
    mm <- c(mm,vec(sub1[i,]))
  }
  A_tau[,j] <- c(mm,vec(sub2))
}

#=====================drawing tau,lambda=================================
mat <- A_tau %*% solve(sim_mat) %*% t(A_tau)
bl <- NULL
for ( i in 1: p1)
{
  bl[[i]] <- list()
  
  for (j in 1: p1)
  {
    bl[[i]][[j]] <- mat[(((i-1)*num)+1):(i*num),(((j-1)*num)+1):(j*num)]
  }
}

for(j in 1:p1){
  lambda[j] <- sqrt(rgamma(1,(((num*npop)+1)/2+alpha),(tau[j]/2 + beta)))
  nu <- lambda[j]/sqrt(tr(solve(U_mat)%*%bl[[j]][[j]]))
  tau[j] <- 1/rinvGauss(1, nu,lambda[j]^2 )
 }

for(j in 1:(p1*(p1-num))){
  lambda <- sqrt(rgamma(1,((npop+1)/2+alpha),(tau[p1+j]/2 + beta)))
   nu <- lambda/sqrt(mat[(num*p1)+j,(num*p1)+j])
  tau[p1+j] <- 1/rinvGauss(1, nu,lambda^2 )
}

#=====================updating U_mat===============================
sum <- matrix(0,nrow=num,ncol=num)
for(j in 1: p1){
  sum <- sum + (1/tau[j])*bl[[j]][[j]]
}
U_nu <- df + (p1*npop)
U_S <- sum + Q2
U_mat <- rinvwishart(U_nu, U_S)
#=====================updating d_tau===============================
mat <- list()
for(i in 1:p1){
  mat[[i]] <- tau[i]*U_mat
}
tau2 <- tau[-c(1:p1)]
d_tau <- bdiag(bdiag(mat),diag(tau2))
#================drawing similarity matrix===============================
sim_nu <- df + (p1^2)
sim_S <- (t(A_tau) %*% solve(d_tau) %*% A_tau) + Q
sim_mat <- rinvwishart(sim_nu, sim_S)
#================updating w,gamma_list,d=======================

cl <- makeCluster(npop, type = "SOCK")
j <- 1:npop
clusterExport(cl,"j")

upd_gamma <- parApply(cl, M,1, rep_gamma,s1_list,s2_list,sigma)
upd_gamma[[1]]

stopCluster(cl)

################################################

gamma4_list <- list()
d_A22 <- list()

for(j in 1:npop){
  gamma4_list[[j]] <- upd_gamma[[j]]$gamma_A22
  d_A22[[j]] <- upd_gamma[[j]]$d_A22
  
}


#================drawing error covariance:sigma==========================
para <- matrix(0,nrow=p,ncol=p)  
for(j in 1: npop){
  para <- para + (s3_list[[j]] - s2_list[[j]]%*% t(A_list[[j]]) - t(s2_list[[j]]%*% t(A_list[[j]])) + A_list[[j]]%*%(s1_list[[j]]%*%t(A_list[[j]]) ))  
  
}

sigma <- rinvwishart(nu=((nobs*npop)+df), S=(para + Q1))

if (ind > burn) 
{
  count <- count + 1
  sigma_sim [ , , count] <- as.matrix(sigma)
  sigma_final <- sigma_final + sigma
  sim_mat_sim [ , , count] <- sim_mat
  sim_mat_final <-sim_mat_final + sim_mat

  for(j in 1: npop){
  A_sim[[j]][, , count] <- A_list[[j]]
  A_final[[j]] <- A_final[[j]] + A_list[[j]]
  }
 
 
}
}

for(j in 1:npop){
  A_final[[j]] = A_final[[j]]/count
  
}
sigma_final = sigma_final/count
sim_mat_final = sim_mat_final/count

#============================== forecast ===================================
crps_mat <- list() ; logs_mat <- list()
ferr_med_low <- matrix(0,nrow=horizon,ncol=npop)
for(j in 1:npop)
{
pred <- matrix(0,nrow=((3*k1)+k2),ncol=nobs+1)
err <- matrix(0,nrow=((3*k1)+k2),ncol=nobs+1)
for ( t in 2:(nobs+1))
{
  pred[,t] <- as.vector((A_final[[j]] %*% data[[j]][,t-1]))
  err[,t] <-data[[j]][,t] - pred[,t] 
  
}
lasso <- CVglasso(X = t(err), S = NULL, nlam = 10, lam.min.ratio = 0.01,
                  lam = NULL, diagonal = FALSE, path = FALSE, tol = 1e-04,
                  maxit = 10000, adjmaxit = NULL, K = 5, crit.cv = c("loglik", "AIC",
                                                                     "BIC"), start = c("warm", "cold"), cores = 1)

cov_mat <- lasso$Sigma
pred_sim <- array(0, c(p, horizon, (iteration-burn)))
for(i in 1:(iteration-burn)){
  pred_y <- list()
  pred_y[[1]] <- mvrnorm(1,(A_sim[[j]][,,i] %*% data[[j]][,nobs+1]),cov_mat)
  for ( t in 2:horizon)
  {
    pred_y[[t]] <- mvrnorm(1,(A_sim[[j]][,,i] %*% pred_y[[t-1]]),cov_mat)
    
  }
  pred_y <- as.matrix(do.call(cbind,pred_y))
  
  pred_sim [ , , i] <- pred_y
}
mean <- apply(pred_sim,c(1,2),mean)
sd <- apply(pred_sim,c(1,2),sd)
crps_mat[[j]] <- matrix(0,nrow = ((3*k1)+k2),ncol=horizon)
logs_mat[[j]] <- matrix(0,nrow = ((3*k1)+k2),ncol=horizon)
for ( i in 1:horizon)
{
  crps_mat[[j]][,i] <- crps(y = forecast_true[[j]][,i], family = "normal", mean = mean[,i], sd = sd[,i])
  logs_mat[[j]][,i] <- logs(y = forecast_true[[j]][,i], family = "normal", mean = mean[,i], sd = sd[,i])
}#CRPS and LPS values provided by ABPVAR model

forecast <- forecast_cred(pred_sim,forecast_true[[j]],horizon)
ferr_med_low[,j] <- forecast$ferr_med_low

}
#Final CRPS and LPS values are obtained after averaging over variables and replicates.

#PMSE values from the ABPVAR model are computed averaging the forecasting error over estimation samples
#and then taking square-root


# Forecasting performance by the benchmark random walk model

library(StMoMo)
err_rw <- matrix(0,nrow=horizon,ncol=npop)

for(j in 1:npop){
  rw <- mrwd(data[[j]])
  drift <- rw$drift
  
  forecast <- list()
  forecast[[1]] <- as.vector(drift) + data[[j]][,nobs+1] 
  for ( t in 2:horizon)
  {
    forecast[[t]] <- as.vector(drift) + forecast[[t-1]]
    
  }
  forecast <- do.call(cbind,forecast)
  forecast <- as.matrix(forecast) 
  for ( h in 1:horizon)
  {
    err_rw[h,j] <- sum((forecast_true[[j]][,h][((3*k1)+1):((3*k1)+k2)]-forecast[,h][((3*k1)+1):((3*k1)+k2)])^2)/k2
  }#forecast error from the random walk model
  
}
#PMSE values from the random walk model are computed averaging the forecasting error 
#over estimation samples and then taking square-root

#The relative PMSE values are computed normalizing the PMSE values from the ABPVAR model 
#by that of the random walk model
