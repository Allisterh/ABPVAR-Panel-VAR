source("functions.R")
load("data_state_georgia.dat")

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
index <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(index+1)
horizon <- 10
nobs <- 150 - index
npop <- 7
N <- 1000
k1 <- 0
k2 <- 15
p <-  k2
p1 <- k2
alpha <- 1
beta <- 1
lambda <- rep(1,(k1+k2)^2)
Q <- 0.3*diag(npop)
Q1 <- 0.3*diag(p)
df <- 1
theta <- rep(0.5,npop)
iteration <- 2000
burn <- 1000

#=======================================Input=============================================
vec <- permutation(k1,k2)$vec
vec1 <- permutation(k1,k2)$vec1
vec2<- permutation(k1,k2)$vec2
M <- matrix(c(1:npop))
#=======================================Initialize====================================
sigma <- genPositiveDefMat(p, covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),  alphad=1, eta=1, rangeVar=c(0.1,1), lambdaLow=1, ratioLambda=2)$Sigma
sim_mat <- genPositiveDefMat(npop, covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),  alphad=1, eta=1, rangeVar=c(0.1,1), lambdaLow=1, ratioLambda=2)$Sigma
tau <- rgamma((p1^2),1.5,1.5)


for(i in 1:npop){
  for(j in 1: p)
  {
    lf[[i]][j,] <- lf[[i]][j,] - mean(lf[[i]][j,])
   
  }
}
#======================Initialize...Gamma matrices for_all_pop==================#change to for lopp
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


#===================================Gibbs=====================================

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

count <- 0

for(ind in 1:iteration){
 
#========covariance matrix for prior of A for all pop==========
cov_A <- kronecker(sim_mat,as.matrix(diag(tau)))

#=================Extracting the part for A12 and A22 for all pop==========
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

#===========================for drawing A11 and A12 for all pop=======================

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
#=====================drawing tau,lambda=================================
mat <- A %*% solve(sim_mat) %*% t(A)

for(j in 1:(k1+k2)^2){
  lambda[j] <- sqrt(rgamma(1,((npop+1)/2+alpha),(tau[j]/2 + beta)))
  nu <- lambda[j]/sqrt(mat[j,j])
  tau[j] <- 1/rinvGauss(1, nu,lambda[j]^2 )
  #lambda[j] <- 1
  #ch <- c(ch,nu)
}
#sqrt(rgamma(1,((npop+1)/2+alpha),(tau[j]/2 + beta)))
#================drawing similarity matrix===============================
sim_nu <- df + (k1+k2)^2
sim_S <- t(A) %*% (as.matrix(diag(1/tau))) %*% A + Q
sim_mat <- rinvwishart(sim_nu, sim_S)
#================drawing theta and updating w,gamma_list,d=======================

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
print(ind)
print(proc.time()[3])
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
}

forecast <- forecast_cred(pred_sim,forecast_true[[j]],horizon)
ferr_med_low[,j] <- forecast$ferr_med_low

}

ferr_med_low

