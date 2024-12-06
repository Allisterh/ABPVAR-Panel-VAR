setwd("C:/Backup/2nd project/code/simulation/Gr2/s1_new")
library(matrixcalc)
library(Matrix)

source("function_v2.R")


k1 <- 20
k2 <- 10
p1 <- k1+k2
p <- (3*k1)+k2
npop <- 10
target <- 0.7 #change
#n1 <- 25
#n2 <- 800
#sd1 <- 1.5
#sd2 <- 0.001
sd3 <- 0.01
nobs <- 300
theta_true <- 0.5
horizon <-10
N <- 1000
alpha <- 1
beta <- 1
lambda <- rep(1,(k1+k2)^2)
Q <- 0.3*diag(npop)
Q1 <- 0.3*diag(p)
df <- 1
theta <- 0.5
num <- 5

#d_tau <- truth_common(k1,k2,target,num)
#d_tau <- d_tau*0.01
#save(d_tau,file=paste("common", ".dat", sep=''))
load("common.dat")

for(i in 1:npop){
  # Create A
  p1 <- k1+k2
  z <- rnorm(p1^2)
  rand_mtx <- matrix(z,nrow=p1,ncol=p1)
  x <- matrix(sd3,nrow=p1,ncol=p1)
  mat <- hadamard.prod(x,rand_mtx)
  A <- d_tau + mat
  eig_A = max(abs(eigen(A)$values))
  eig_A
  
  temp <- max(abs(eigen(A[1:k1,1:k1])$values))
  temp
  count = 0
  while (temp > 0.61){ # should be false..no further shrinkage in A
    count = count+1
    A[1:k1,1:k1] = A[1:k1,1:k1]*0.95
    temp = max(abs(eigen(A[1:k1,1:k1])$values)) 
  }  
  
  target_A = 1 - (2*(theta^2)*0.61)
  count = 0
  eig_A = max(abs(eigen(A)$values))
  while (eig_A > target_A){# should be false..no further shrinkage in A
    count = count+1
    A[(k1+1):p1,] = A[(k1+1):p1,]*0.95
    A[1:k1,(k1+1):p1] = A[1:k1,(k1+1):p1]*0.95
    eig_A = max(abs(eigen(A)$values))
  }
  max(abs(eigen(A)$values))
  A[1:k1,1:(k1+k2)] <- A[1:k1,1:(k1+k2)]*(theta ^ 2)
  mu <- c(1/(theta* theta), 1/theta) 
  gamma <- c((theta* theta), theta)
  #eig_A = max(abs(eigen(A)$values))
  W_block_two <-kronecker(A[1:k1,], mu)
  W_block_three <-kronecker(A[,1:k1], t(gamma))
  v <- kronecker(t(gamma), mu)
  W_block_four <- kronecker(A[1:k1,1:k1], v)
  W1 <- rbind(A, W_block_two)
  W2 <- rbind(W_block_three, W_block_four)
  W <- cbind(W1,W2)
  max(abs(eigen(W)$values))
  A[1:k1,1:(k1+k2)] <- A[1:k1,1:(k1+k2)]/(theta ^ 2)
  
  
  #A_true <- list();w_true <- list()
  
  A_true[[i]] <- A
  w_true[[i]] <- W
}
# Save sequentially for all A and W.

eig <- numeric()
eig1 <- numeric()
for(j in 1:npop){
  eig <- c(eig,max(abs(eigen(w_true[[j]])$values)))
  eig1 <- c(eig1,max(abs(eigen(A_true[[j]])$values)))
}
eig
eig1


save(A_true,file=paste("A_true_0.01_10", ".dat", sep=''))
save(w_true, file=paste("w_true_0.01_10", ".dat", sep=''))


#####################
load("w_true_0.01.dat")

temp_w <- w_true[[1]]
#======================================================
sigma_true <- Gen_Sigma(k1,k2)$sigma_true
desired_SNR <- 2
k <- norm(temp_w, type="F")/norm(desired_SNR*sigma_true, type="F")
sigma_true <- k * sigma_true
save(sigma_true,file=paste("sigma_true_0.01", ".dat", sep=''))
