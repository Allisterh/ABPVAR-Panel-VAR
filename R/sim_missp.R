setwd("C:/Users/nchak/OneDrive/Desktop/R_codes")
library(matrixcalc)
library(clusterGeneration)
library(dplyr)
library(MASS)
library(Matrix)
library(glmnet)
library(fBasics)
library(base)
library(SuppDists)
library(LaplacesDemon)
library(mvnfast)
library(snow)
library(scoringRules)
library("CVglasso")
source("function_states.R")
################  Generating true VAR transition matrices (with mild heterogeneity)
k1 <- 0
k2 <- 30
p1 <- k1+k2
p <- (3*k1)+k2
target <- 0.7
hetr <- 1.5
nobs <- 300
npop <- 6
d_tau_1 <- truth_common(k1,k2,target)*10
tt <- truth_gen_1(k1,k2,d_tau_1)
Am_true <- list()
for(i in 1:5){
  Am_true[[i]] <- tt 
}
Am_true[[6]] <- tt*0.1

############## Generating true VAR transition matrices (with strong heterogeneity)

d_tau_1 <- truth_common(k1,k2,target)*10
block_structure_matrix <- function(K, blocks = 2, noise = 0.05) {
  A <- matrix(0, K, K)
  block_size <- floor(K / blocks)
  
  for (b in 0:(blocks - 1)) {
    i_start <- b * block_size + 1
    i_end <- min((b + 1) * block_size, K)
    block <- matrix(runif((i_end - i_start + 1)^2, 0.5, 1.0),
                    nrow = i_end - i_start + 1)
    A[i_start:i_end, i_start:i_end] <- block
  }
  
  A <- A + matrix(rnorm(K * K, mean = 0, sd = noise), K, K)
  return(A)
}
A1 <- block_structure_matrix(30)

tt1 <- truth_gen_1(k1,k2,A1)
tt2 <- truth_gen_1(k1,k2,d_tau_1)

A_true <- list()
for(i in 1:5){
  A_true[[i]] <- tt1
}
A_true[[6]] <- tt2*0.5

##### Estimation using data generated from transition matrices with strong heterogeneity
#For the misspecified case that includes the 6th entity
k1 <- 0
k2 <- 30
p1 <- k1+k2
p <- (3*k1)+k2
npop <- 6
nobs <- 300
theta_true <- rep(0.5,npop)
horizon <-10
N <- 1000
alpha <- 1
beta <- 1
lambda <- rep(1,(k1+k2)^2)
Q <- 0.3*diag(npop)
Q1 <- 0.3*diag(p)
df <- 1
theta <- rep(0.5,npop)
iteration <- 2000
burn <- 1000
initial_Theta <- genPositiveDefMat(npop)$Sigma

#=======================================Input=============================================
vec <- permutation(k1,k2)$vec
vec1 <- permutation(k1,k2)$vec1
vec2<- permutation(k1,k2)$vec2
M <- matrix(c(1:npop))
#======================================================
sigma_true <- Gen_Sigma(k2)$sigma_true
desired_SNR <- 2
norms <- sapply(A_true, function(mat) norm(mat, type = "F"))
which_min <- which.min(norms)
k <- norm(A_true[[which_min]], type="F")/norm(desired_SNR*sigma_true, type="F")
sigma_true <- k * sigma_true
#=======================================Initialize====================================
sigma <- genPositiveDefMat(p, covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),  alphad=1, eta=1, rangeVar=c(0.1,1), lambdaLow=1, ratioLambda=2)$Sigma
sim_mat <- initial_Theta
tau <- rgamma((p1^2),1.5,1.5)

gamma4_list <- list()
s1_list <- list();s2_list <- list();s3_list <- list()
data <- list();forecast_true <- list();d_A22 <- list()

for(j in 1:npop){
  tt <- data_gen(A_true[[j]],k1,k2,horizon,sigma_true)
  data[[j]] <- tt$data
  forecast_true[[j]] <- tt$forecast_true
  b4 <- numeric((k2*k2))
  s_func <- S_matrix(data[[j]],nobs)
  s1_list[[j]] <- s_func$s1
  s2_list[[j]]<- s_func$s2
  s3_list[[j]] <- s_func$s3
  A22_func <- para_A22(s1_list[[j]],s2_list[[j]],sigma)
  gamma4_list[[j]] <- A22_func$gamma_A22
  d_A22[[j]] <- A22_func$d_A22
}
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
tau_sim <- array(0, c((p1^2), 1, (iteration-burn)))
lambda_sim <- array(0, c((p1^2),1, (iteration-burn)))

count <- 0
for(ind in 1:iteration){
  
  cov_A <- kronecker(sim_mat,as.matrix(diag(tau)))
  
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
  #=====================drawing tau,lambda=================================
  mat <- A %*% solve(sim_mat) %*% t(A)
  
  for(j in 1:(k1+k2)^2){
    lambda[j] <- sqrt(rgamma(1,((npop+1)/2+alpha),(tau[j]/2 + beta)))
    nu <- lambda[j]/sqrt(mat[j,j])
    tau[j] <- 1/rinvGauss(1, nu,lambda[j]^2 )
  }
  #================drawing similarity matrix===============================
  sim_nu <- df + (k1+k2)^2
  cc <- t(A) %*% (as.matrix(diag(1/tau))) %*% A
  sim_S <- cc + Q
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
    tau_sim [ , , count] <- tau
    lambda_sim [ , , count] <- lambda
    
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
  }
  
  forecast <- forecast_cred(pred_sim,forecast_true[[j]],horizon)
  ferr_med_low[,j] <- forecast$ferr_med_low
  
}
#PMSE values (for this misspecified case) are computed averaging the forecasting error over replicates
#and then taking square-root

#The estimation part of the code can be repeated for the first five entities and
#PMSE values can be obtained in a similar way, which can be used for normalizing the
#PMSE values obtained from the misspecified case.

