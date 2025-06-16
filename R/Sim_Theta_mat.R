setwd("C:/Users/nchak/OneDrive/Desktop/R_codes")
library(matrixcalc)
library(clusterGeneration)
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

source("function_states.R")
##################### Generating true VAR transition matrices
k1 <- 0
k2 <- 30
p1 <- k1+k2
p <- (3*k1)+k2
target <- 0.7
hetr <- 1.5
npop <- 50
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

d_tau_1 <- truth_common(k1,k2,target)*10
A_true <- list()
tt <- truth_gen_1(k1,k2,d_tau_1)
rate <- 0.02
for(i in 1:npop){
  A_true[[i]] <- tt *rate*(50-(i-1)) 
}
initial_Theta <- genPositiveDefMat(npop)$Sigma

###################### Estimation
k1 <- 0
k2 <- 30
p1 <- k1+k2
p <- (3*k1)+k2
npop <- 50
target <- 0.7
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
#===================================Gibbs sampling=====================================

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
  
  #========covariance matrix for prior of A for all pop==========
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
    sim_mat_final <- sim_mat_final + sim_mat
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

####### Estimated correlation matrix
cor_mat <- matrix(0,nrow=npop,ncol=npop)
for(i in 1:50){
  for(j in 1:50){
    cor_mat[i,j]<- sim_mat_final[i,j]/(sqrt(sim_mat_final[i,i])*sqrt(sim_mat_final[j,j]))
  }
}


