#setwd("C:/Users/user/Dropbox (UFL)/2nd project/code/Real data")
#source("truth_gen.R")
source("function_v2.R")
load("data_new_s3.dat")


start <- proc.time()[3]
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
#index=1
set.seed(index+1)
horizon <- 2
nobs <- 80 - index
npop <- 3
N <- 1000
k1 <- 10
k2 <- 5
p <- (3*k1) + k2
p1 <- (k1+k2)

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

sigma_inv <- sigma_inverse(sigma,k1,k2)$sigma_inv
sigma_11 <- sigma_inverse(sigma,k1,k2)$sigma_11
sigma_22 <- sigma_inverse(sigma,k1,k2)$sigma_22


for(i in 1:npop){
  for(j in 1: p)
  {
    data_full[[i]][j,] <- data_full[[i]][j,] - mean(data_full[[i]][j,])
   
  }
}
#======================Initialize...Gamma matrices for_all_pop==================#change to for lopp
gamma1_list <- list()
gamma2_list <- list()
gamma3_list <- list()
gamma4_list <- list()
s_list <- list();s1_list <- list();s2_list <- list();s3_list <- list()
s11_list <- list();s12_list <- list();s22_list <- list();ss11_list <- list()
ss12_list <- list();ss21_list <- list();ss22_list <- list();data <- list();forecast_true <- list()
d_A11 <- list()
d_A12 <- list()
d_A21 <- list()
d_A22 <- list()

for(j in 1:npop){
  data[[j]] <- data_full[[j]][,1:(nobs + 1)]
  forecast_true[[j]] <- data_full[[j]][,(nobs + 1 + 1):(nobs + 1 + horizon)]
  X <- as.matrix(t(data[[j]])[1:nobs,])
  Y <- as.matrix(t(data[[j]])[2:(nobs + 1),])
  LSE <- as.matrix(lm(Y ~ X - 1)$coef)
  A11_initial <- t(LSE)[1:k1,1:k1]
  b1<-  vec(A11_initial)
  b2 <- numeric((k1*k2))
  b3 <- numeric((k1*k2))
  b4 <- numeric((k2*k2))
  #A22_initial <- t(LSE)[(3*k1+1):(3*k1+k2),(3*k1+1):(3*k1+k2)]
  #b4<- as.vector(vec(A22_initial))
  theta_mat <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$theta_mat
  u <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$u
  v <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$v
  w11 <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$w11
  w12 <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$w12
  w21 <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$w21
  w22 <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$w22
  s_list[[j]] <- S_matrix(data[[j]],nobs,k1,k2)$s
  s1_list[[j]] <- S_matrix(data[[j]],nobs,k1,k2)$s1
  s2_list[[j]] <- S_matrix(data[[j]],nobs,k1,k2)$s2
  s3_list[[j]] <- S_matrix(data[[j]],nobs,k1,k2)$s3
  s11_list[[j]] <- S_matrix(data[[j]],nobs,k1,k2)$s11
  s12_list[[j]] <- S_matrix(data[[j]],nobs,k1,k2)$s12
  s22_list[[j]] <- S_matrix(data[[j]],nobs,k1,k2)$s22
  ss11_list[[j]] <- S_matrix(data[[j]],nobs,k1,k2)$ss11
  ss12_list[[j]] <- S_matrix(data[[j]],nobs,k1,k2)$ss12
  ss21_list[[j]] <- S_matrix(data[[j]],nobs,k1,k2)$ss21
  ss22_list[[j]] <- S_matrix(data[[j]],nobs,k1,k2)$ss22
  A11_func <- para_A11(k1,k2,theta_mat,s_list[[j]],sigma_inv,sigma_11,s11_list[[j]],ss11_list[[j]],s12_list[[j]],w12)
  gamma1_list[[j]] <- A11_func$gamma_A11
  d_A11[[j]] <- A11_func$d_A11
  A12_func <- para_A12(k1,k2,u,sigma_inv,sigma_11,s12_list[[j]],ss12_list[[j]],s22_list[[j]],w11)
  gamma2_list[[j]] <- A12_func$gamma_A12
  d_A12[[j]] <- A12_func$d_A12
  A21_func <- para_A21(k1,k2,v,sigma_22,s_list[[j]],s11_list[[j]],ss21_list[[j]],s12_list[[j]],w22)
  gamma3_list[[j]] <- A21_func$gamma_A21
  d_A21[[j]] <- A21_func$d_A21
  A22_func <- para_A22(v,sigma_22,s22_list[[j]],ss22_list[[j]],s12_list[[j]],w21)
  gamma4_list[[j]] <- A22_func$gamma_A22
  d_A22[[j]] <- A22_func$d_A22
}#check


#===================================Gibbs=====================================
w_sim <- list()
w_final <- list()
A_sim <- list()
A_final <- list()

for(j in 1:npop){
  w_sim[[j]] <- array(0, c(p, p, (iteration-burn)))
  w_final[[j]] <- matrix(0,p,p)
  A_sim[[j]] <- array(0, c(p1, p1, (iteration-burn)))
  A_final[[j]] <- matrix(0,p1,p1)
  
}
theta_sim <- array(0, c(npop, 1, (iteration-burn)))
sigma_sim <- array(0, c(p, p, (iteration-burn)))
sigma_final <- matrix(0,p,p)
theta_final <- rep(0,npop)

count <- 0

for(ind in 1:iteration){
 
#========covariance matrix for prior of A for all pop==========
cov_A <- kronecker(sim_mat,as.matrix(diag(tau)))

#=================Extracting the part for A11 and A21 for all pop==========
vec_A11 <- list()
vec_A21 <- list()
for(m in 1:npop)
{
  temp_A11 <- NULL
  temp_A21 <- NULL
  for(i in 1:k1){
    m1 <- c((((m-1)*((k1+k2)^2))+((i-1)*(k1+k2))+1):(((m-1)*((k1+k2)^2))+((i-1)*(k1+k2))+k1)) #for A11
    m2 <- c((((m-1)*((k1+k2)^2))+((i-1)*(k1+k2))+k1+1):(((m-1)*((k1+k2)^2))+((i-1)*(k1+k2))+k1+k2)) #for A21
    temp_A11 <- c(temp_A11,m1)
    temp_A21 <- c(temp_A21,m2)
  }
  vec_A11[[m]] <- temp_A11  
  vec_A21[[m]] <- temp_A21 
}
sig_A11 <- list()
sig_A21 <- list()
for(j in 1:npop)
{
  sig_A11[[j]] <- as.matrix(cov_A[vec_A11[[j]],])
  sig_A11[[j]] <- as.matrix(sig_A11[[j]][,vec_A11[[j]]])
  sig_A21[[j]] <- as.matrix(cov_A[vec_A21[[j]],])
  sig_A21[[j]] <- as.matrix(sig_A21[[j]][,vec_A21[[j]]])
  
}
#===========================for drawing A11 and A12 for all pop=======================
b1_all <- NULL
b3_all <- NULL

for(j in 1:npop)
{
  cov_A11 <- gamma1_list[[j]] + solve((sig_A11[[j]]))
  b1 <- rmvn(1,solve((cov_A11))%*%d_A11[[j]],solve((cov_A11)))
  b1_all <- c(b1_all,b1)
  cov_A21 <- gamma3_list[[j]] + solve((sig_A21[[j]]))
  b3 <- rmvn(1,solve((cov_A21))%*%d_A21[[j]],solve((cov_A21)))
  b3_all <- c(b3_all,b3)
}
b1 <- matrix(b1_all,nrow=k1*k1,ncol=npop)
b3 <- matrix(b3_all,nrow=k1*k2,ncol=npop)

#=================Extracting the part for A12 and A22 for all pop==========
vec_A12 <- list()
vec_A22 <- list()
for(m in 1:npop)
{
  temp_A12 <- NULL
  temp_A22 <- NULL
  for(i in 1:k2){
    m1 <- c((((m-1)*((k1+k2)^2))+ k1*(k1+k2) +((i-1)*(k1+k2))+1):(((m-1)*((k1+k2)^2))+ k1*(k1+k2) +((i-1)*(k1+k2))+k1)) #for A12
    m2 <- c((((m-1)*((k1+k2)^2))+k1*(k1+k2) +((i-1)*(k1+k2))+k1+1):(((m-1)*((k1+k2)^2))+k1*(k1+k2) +((i-1)*(k1+k2))+k1+k2)) #for A22
    temp_A12 <- c(temp_A12,m1)
    temp_A22 <- c(temp_A22,m2)
  }
  vec_A12[[m]] <- temp_A12  
  vec_A22[[m]] <- temp_A22  
}

sig_A12 <- list()
sig_A22 <- list()
for(j in 1:npop)
{
  sig_A12[[j]] <- as.matrix(cov_A[vec_A12[[j]],])
  sig_A12[[j]] <- as.matrix(sig_A12[[j]][,vec_A12[[j]]])
  sig_A22[[j]] <- as.matrix(cov_A[vec_A22[[j]],])
  sig_A22[[j]] <- as.matrix(sig_A22[[j]][,vec_A22[[j]]])
  
}

#===========================for drawing A11 and A12 for all pop=======================

b2_all <- NULL
b4_all <- NULL

for(j in 1:npop)
{
  cov_A12 <- gamma2_list[[j]] + solve((sig_A12[[j]]))
  b2 <- rmvn(1,solve((cov_A12))  %*%d_A12[[j]],solve((cov_A12)))
  b2_all <- c(b2_all,b2)
  cov_A22 <- gamma4_list[[j]] + solve((sig_A22[[j]]))
  b4 <- rmvn(1,solve((cov_A22))%*%d_A22[[j]],solve((cov_A22)))
  b4_all <- c(b4_all,b4)
}
b2 <- matrix(b2_all,nrow=k1*k2,ncol=npop)
b4 <- matrix(b4_all,nrow=k2*k2,ncol=npop)

#=================drawing tau==================
A11_list <- list()
A12_list <- list()
A21_list <- list()
A22_list <- list()
A_list <- list()
A <- matrix(0,nrow=(k1+k2)^2,ncol=npop)
#=======================A matrix for all population=======================
for(j in 1:npop){
  A11_list[[j]]<- matrix(b1[,j],nrow=k1,ncol=k1)
  A12_list[[j]]<- matrix(b2[,j],nrow=k1,ncol=k2)
  A21_list[[j]]<- matrix(b3[,j],nrow=k2,ncol=k1)
  A22_list[[j]]<- matrix(b4[,j],nrow=k2,ncol=k2)
  A_list[[j]]<- cbind(rbind(A11_list[[j]],A21_list[[j]]),rbind(A12_list[[j]],A22_list[[j]]))
  A[,j] <- vec(A_list[[j]]) #dim of A: (k1+k2)^2 by npop
}
#=====================drawing tau,lambda=================================
mat <- A %*% solve(sim_mat) %*% t(A)
for(j in 1:(k1+k2)^2){
  lambda[j] <- sqrt(rgamma(1,((npop+1)/2+alpha),(tau[j]/2 + beta)))
#  nu <- lambda[j]/norm_vec(A[j,]*sqrt(tr(solve(sim_mat))))
   nu <- lambda[j]/sqrt(mat[j,j])
 tau[j] <- 1/rinvGauss(1, nu,lambda[j]^2 )
  #lambda[j] <- 1
  #ch <- c(ch,nu)
}
#================drawing similarity matrix===============================
sim_nu <- df + (k1+k2)^2
sim_S <- t(A) %*% (as.matrix(diag(1/tau))) %*% A + Q
sim_mat <- rinvwishart(sim_nu, sim_S)
#================drawing theta and updating w,gamma_list,d=======================
##############################Parallel computing#################################
cl <- makeCluster(npop, type = "SOCK")
j <- 1:npop
clusterExport(cl,"j")

upd_dtheta <- parApply(cl, M,1, rep_dtheta,s1_list,s2_list,s3_list,
                   sigma,b1,b2,b3,b4,k1,k2,vec1)
upd_dtheta
upd_theta <- parApply(cl, M,1, rep_theta,N,upd_dtheta)
upd_theta
upd_thetamat <- parApply(cl, M,1, rep_thetamat,upd_theta,b1,b2,b3,b4,k1,k2,vec1)
upd_thetamat[[1]]
stopCluster(cl)

theta_mat <- list()
u <- list()
v <- list()
w11 <- list()
w12 <- list()
w21 <- list()
w22 <- list()
w_list <- list()

for(j in 1:npop){
  theta_mat[[j]] <- upd_thetamat[[j]]$theta_mat
  u[[j]] <- upd_thetamat[[j]]$u
  v[[j]] <- upd_thetamat[[j]]$v
  w11[[j]] <- upd_thetamat[[j]]$w11
  w12[[j]] <- upd_thetamat[[j]]$w12
  w21[[j]] <- upd_thetamat[[j]]$w21
  w22[[j]] <- upd_thetamat[[j]]$w22
  w_list[[j]] <- upd_thetamat[[j]]$w
}

cl <- makeCluster(npop, type = "SOCK")
j <- 1:npop
clusterExport(cl,"j")

upd_gamma <- parApply(cl, M,1, rep_gamma,theta_mat,s_list,sigma_inv,sigma_11,s11_list,ss11_list,s12_list,
                      w12,u,ss12_list,w11,v,sigma_22,ss21_list,ss22_list,s22_list,w22,w21,k1,k2)
upd_gamma[[1]]

stopCluster(cl)

################################################

gamma1_list <- list()
gamma2_list <- list()
gamma3_list <- list()
gamma4_list <- list()

d_A11 <- list()
d_A12 <- list()
d_A21 <- list()
d_A22 <- list()

for(j in 1:npop){
  gamma1_list[[j]] <- upd_gamma[[j]]$gamma_A11
  d_A11[[j]] <- upd_gamma[[j]]$d_A11
  
  gamma2_list[[j]] <- upd_gamma[[j]]$gamma_A12
  d_A12[[j]] <- upd_gamma[[j]]$d_A12
  
  gamma3_list[[j]] <- upd_gamma[[j]]$gamma_A21
  d_A21[[j]] <- upd_gamma[[j]]$d_A21
  
  gamma4_list[[j]] <- upd_gamma[[j]]$gamma_A22
  d_A22[[j]] <- upd_gamma[[j]]$d_A22
  
}

#gamma1 <- as.matrix(bdiag(gamma1_list))
#gamma2 <- as.matrix(bdiag(gamma2_list))
#gamma3 <- as.matrix(bdiag(gamma3_list))
#gamma4 <- as.matrix(bdiag(gamma4_list))



#================drawing error covariance:sigma==========================
para <- matrix(0,nrow=p,ncol=p)  
for(j in 1: npop){
  para <- para + (s3_list[[j]] - s2_list[[j]]%*% t( w_list[[j]]) - t(s2_list[[j]]%*% t( w_list[[j]])) + w_list[[j]]%*%(s1_list[[j]]%*%t( w_list[[j]]) ))  
  
}

sigma <- rinvwishart(nu=((nobs*npop)+df), S=(para + Q1))

sigma_inv <- sigma_inverse(sigma,k1,k2)$sigma_inv
sigma_11 <- sigma_inverse(sigma,k1,k2)$sigma_11
sigma_22 <- sigma_inverse(sigma,k1,k2)$sigma_22

if (ind > burn) 
{
  count <- count + 1
  sigma_sim [ , , count] <- as.matrix(sigma)
  theta_sim [ , , count] <- upd_theta
  theta_final <- theta_final + upd_theta
  sigma_final <- sigma_final + sigma
  
  for(j in 1: npop){
  w_sim[[j]][, , count] <- w_list[[j]]
  w_final[[j]] <- w_final[[j]] + w_list[[j]]
  A_sim[[j]][, , count] <- A_list[[j]]
  A_final[[j]] <- A_final[[j]] + A_list[[j]]
  }
 
 
}
print(ind)
print(proc.time()[3])
}

for(j in 1:npop){
  w_final[[j]] = w_final[[j]]/count  
  A_final[[j]] = A_final[[j]]/count
  
}
theta_final = theta_final/count
sigma_final = sigma_final/count

#============================== forecast ===================================
crps_mat <- list() ; logs_mat <- list()
ferr_med_low <- matrix(0,nrow=horizon,ncol=npop)
for(j in 1:npop)
{
pred <- matrix(0,nrow=((3*k1)+k2),ncol=nobs+1)
err <- matrix(0,nrow=((3*k1)+k2),ncol=nobs+1)
for ( t in 2:(nobs+1))
{
  pred[,t] <- as.vector((w_final[[j]] %*% data[[j]][,t-1]))
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
  pred_y[[1]] <- mvrnorm(1,(w_sim[[j]][,,i] %*% data[[j]][,nobs+1]),cov_mat)
  for ( t in 2:horizon)
  {
    pred_y[[t]] <- mvrnorm(1,(w_sim[[j]][,,i] %*% pred_y[[t-1]]),cov_mat)
    
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

#========================save=====================================
save(crps_mat, file=paste("/home/nchakraborty/2nd/real_data/data3/new_set3/res2",
                          "/crps_mat_test-", index, ".dat", sep=''))
save(logs_mat, file=paste("/home/nchakraborty/2nd/real_data/data3/new_set3/res2",
                          "/logs_mat_test-", index, ".dat", sep=''))
save(w_final, file=paste("/home/nchakraborty/2nd/real_data/data3/new_set3/res2",
                       "/w_est_test-", index, ".dat", sep=''))
save(ferr_med_low, file=paste("/home/nchakraborty/2nd/real_data/data3/new_set3/res2",
                                 "/ferr_med_low_test-", index, ".dat", sep=''))

save(forecast, file=paste("/home/nchakraborty/2nd/real_data/data3/new_set3/res2",
                                 "/forecast_test-", index, ".dat", sep=''))

#==============================

crps_mat <- list() ; logs_mat <- list()
ferr_med_low <- matrix(0,nrow=horizon,ncol=npop)
for(j in 1:npop)
{
  pred_sim <- array(0, c(p, horizon, (iteration-burn)))
  for(i in 1:(iteration-burn)){
    pred_y <- list()
    pred_y[[1]] <- as.vector((w_sim[[j]][,,i] %*% data[[j]][,nobs+1]))
    for ( t in 2:horizon)
    {
      pred_y[[t]] <- as.vector((w_sim[[j]][,,i] %*% pred_y[[t-1]]))
      
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

save(crps_mat, file=paste("/home/nchakraborty/2nd/real_data/data3/new_set3/res2",
                          "/crps_v1_test-", index, ".dat", sep=''))
save(logs_mat, file=paste("/home/nchakraborty/2nd/real_data/data3/new_set3/res2",
                          "/logs_v1_test-", index, ".dat", sep=''))
save(ferr_med_low, file=paste("/home/nchakraborty/2nd/real_data/data3/new_set3/res2",
                                 "/ferr_med_v1_test-", index, ".dat", sep=''))
save(data, file=paste("/home/nchakraborty/2nd/real_data/data3/new_set3/res2",
                                 "/data_test-", index, ".dat", sep=''))
save(forecast_true, file=paste("/home/nchakraborty/2nd/real_data/data3/new_set3/res2",
                                 "/forecast_true_test-", index, ".dat", sep=''))
save(sigma_final, file=paste("/home/nchakraborty/2nd/real_data/data3/new_set3/res2",
                                 "/sigma_final_test-", index, ".dat", sep=''))




