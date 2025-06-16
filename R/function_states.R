
norm_vec <- function(x) sqrt(sum(x^2))

truth_common <- function(k1,k2,target){
  p1 <- k1+k2
  tau <- rgamma(p1^2,shape=1.5,rate=1.5)
  tau
  d_tau <- matrix(tau,nrow=p1,ncol=p1)
  d_tau
  #target <- 0.7
  eig = max(abs(eigen(d_tau)$values))
  count = 0
  while (eig>target){
    count = count+1
    d_tau= d_tau*0.95
    eig = max(abs(eigen(d_tau)$values))
  } 
  d_tau 
}
truth_gen_1 <- function(k1,k2,A)
{
  eig_A = max(abs(eigen(A)$values))
  eig_A
  while (eig_A > 0.9){
    A = A*0.95
    eig_A = max(abs(eigen(A)$values)) 
  }  
  
  A
}



#=======================================Generate Sigma Matrix=================================
Gen_Sigma = function(k2)
{
  sigma_true <- genPositiveDefMat(k2, covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),  alphad=1, eta=1, rangeVar=c(0.1,1), lambdaLow=1, ratioLambda=2)$Sigma
  return(list(sigma_true=sigma_true))
  
}


#------permutation vectors----------

permutation <- function(k1,k2)
{
  y0<- NULL
  for(i1 in 1:k1)
  {
    x <- i1
    y0 <- c(y0,x)
  }
  y0
  
  y1<- NULL
  for(i1 in 1:k1)
  {
    x <- k1+k2+(2*i1-1)
    y1 <- c(y1,x)
  }
  y1
  
  y2<- NULL
  for(i1 in 1:k1)
  {
    x <- k1+k2+(2*i1)
    y2 <- c(y2,x)
  }
  y2
  
  z<- NULL
  for(i1 in 1:k2)
  {
    x <- k1+i1
    z <- c(z,x)
  }
  z
  vec <- c(y0,y1,y2,z)
  
  m1<- NULL
  for(i1 in 1:k1)
  {
    m <- c(k1+i1,2*k1+i1,i1)
    m1 <- c(m1,m)
  }
  m1
  
  m2<- NULL
  for(i1 in 1:k2)
  {
    m <- 3*k1+i1
    m2 <- c(m2,m)
  }
  
  vec1 <- c(m1,m2)
  
  a<- NULL
  for(i1 in 1:k1)
  {
    m <- 3*i1
    a <- c(a,m)
  }
  a
  
  b<- NULL
  for(i1 in 1:k1)
  {
    m <- 3*i1 - 2
    b <- c(b,m)
  }
  b
  
  d<- NULL
  for(i1 in 1:k1)
  {
    m <- 3*i1 - 1
    d <- c(d,m)
  }
  d
  
  m2<- NULL
  for(i1 in 1:k2)
  {
    m <- 3*k1+i1
    m2 <- c(m2,m)
  }
  m2
  vec2 <- c(a,b,d,m2)
  return(list(vec=vec,vec1=vec1,vec2=vec2))
}


data_gen <- function(W,k1,k2,horizon,sigma_true)
{
  p <- 3*k1 + k2   
  y_true <- list()
  y <- list()
  y_mod <- list()
  err <- list()
  y[[1]] <- mvrnorm(1,rep(0,p),sigma_true)
  for ( t in 2:(nobs + horizon + 1))
  {
    err[[t]] <- mvrnorm(1,rep(0,p),sigma_true)
    y[[t]] <- as.vector((W %*% y[[t-1]]) + err[[t]])
     # y_mod[[t]]<-y[[t]][vec1]
    
  }
  
  y_full <- as.matrix(do.call(cbind,y))
  #y_mod_full <- do.call(cbind,y_mod)
  data <- y_full[,1:(nobs + 1)]
  forecast_true <- y_full[,(nobs + 1 + 1):(nobs + 1 + horizon)]
  
  #data_mod <-as.matrix(y_mod_full)[,1:(nobs + 1)]
  #data_t <- t(data)
  return(list(data=data,forecast_true=forecast_true))
  
}


S_matrix <- function(data,n)
{
  s1 <- data[,1:n] %*% t(data[,1:n])
  s2 <- data[,2:(n+1)] %*% t(data[,1:n])
  s3 <- data[,2:(n+1)] %*% t(data[,2:(n+1)])
  
  return(list(s1=s1,s2=s2,s3=s3))
}

para_A22 <- function(s1,s2,sigma)
{
  
  mean <- s2 %*% solve(s1)
  m <- vec(mean)
  gamma_A22 <- solve(kronecker(solve(s1),sigma))
  d_A22 <- gamma_A22 %*% m
  return(list(gamma_A22=gamma_A22,d_A22=d_A22))
}


rep_gamma <- function(j,s1_list,s2_list,sigma)
{
  
  library(fBasics)
  
  #=========function "para_A22"======================
  
  para_A22 <- function(s1,s2,sigma)
  {
    
    mean <- s2 %*% solve(s1)
    m <- vec(mean)
    gamma_A22 <- solve(kronecker(solve(s1),sigma))
    d_A22 <- gamma_A22 %*% m
    return(list(gamma_A22=gamma_A22,d_A22=d_A22))
  }
  
  gamma_A22 <- para_A22(s1_list[[j]],s2_list[[j]],sigma)$gamma_A22
  d_A22 <- para_A22(s1_list[[j]],s2_list[[j]],sigma)$d_A22
  
  return(list(gamma_A22=gamma_A22,d_A22=d_A22))
  
  
}

forecast_cred <- function(sim,true,horizon)
{
  sim_q11 <- apply(sim,c(1,2),quantile,probs=0.025)
  sim_q1 <- apply(sim,c(1,2),quantile,probs=0.25)
  sim_med <- apply(sim,c(1,2),quantile,probs=0.50)
  sim_q3 <- apply(sim,c(1,2),quantile,probs=0.75)
  sim_q33 <- apply(sim,c(1,2),quantile,probs=0.975)
  q11 <- vec(sim_q11);q1 <- vec(sim_q1);med <- vec(sim_med);q3 <- vec(sim_q3);q33 <- vec(sim_q33)
  err <- numeric()
  for ( i in 1:horizon)
  {
    err[i] <- (sum((true[,i]-sim_med[,i])^2))/length(true[,i])
  }
  err_low <- numeric()
  for ( j in 1:horizon)
  {
    err_low[j] <- sum((true[,j][((3*k1)+1):((3*k1)+k2)]-sim_med[,j][((3*k1)+1):((3*k1)+k2)])^2)/k2
  }  
  return(list(forecast_med = sim_med,ferr_med=err,ferr_med_low=err_low ))
}

cred_A <- function(sim,est)
{
  sim_q11 <- apply(sim,c(1,2),quantile,probs=0.025)
  sim_q1 <- apply(sim,c(1,2),quantile,probs=0.25)
  sim_q2 <- apply(sim,c(1,2),quantile,probs=0.50)
  sim_q3 <- apply(sim,c(1,2),quantile,probs=0.75)
  sim_q33 <- apply(sim,c(1,2),quantile,probs=0.975)
  cred_A11 <- cbind(vec(est[1:k1,1:k1]),vec(sim_q11[1:k1,1:k1]),vec(sim_q1[1:k1,1:k1]),vec(sim_q2[1:k1,1:k1]),
                    vec(sim_q3[1:k1,1:k1]),vec(sim_q33[1:k1,1:k1]))
  colnames(cred_A11)=c("est","2.5%","25%","50%","75%","97.5%")
  cred_A12 <- cbind(vec(est[1:k1,(k1+1):(k1+k2)]),vec(sim_q11[1:k1,(k1+1):(k1+k2)]),vec(sim_q1[1:k1,(k1+1):(k1+k2)]),vec(sim_q2[1:k1,(k1+1):(k1+k2)]),
                    vec(sim_q3[1:k1,(k1+1):(k1+k2)]),vec(sim_q33[1:k1,(k1+1):(k1+k2)]))
  colnames(cred_A12)=c("est","2.5%","25%","50%","75%","97.5%")
  cred_A21 <- cbind(vec(est[(k1+1):(k1+k2),1:k1]),vec(sim_q11[(k1+1):(k1+k2),1:k1]),vec(sim_q1[(k1+1):(k1+k2),1:k1]),vec(sim_q2[(k1+1):(k1+k2),1:k1]),
                    vec(sim_q3[(k1+1):(k1+k2),1:k1]),vec(sim_q33[(k1+1):(k1+k2),1:k1]))
  colnames(cred_A21)=c("est","2.5%","25%","50%","75%","97.5%")
  cred_A22 <- cbind(vec(est[(k1+1):(k1+k2),(k1+1):(k1+k2)]),vec(sim_q11[(k1+1):(k1+k2),(k1+1):(k1+k2)]),vec(sim_q1[(k1+1):(k1+k2),(k1+1):(k1+k2)]),vec(sim_q2[(k1+1):(k1+k2),(k1+1):(k1+k2)]),
                    vec(sim_q3[(k1+1):(k1+k2),(k1+1):(k1+k2)]),vec(sim_q33[(k1+1):(k1+k2),(k1+1):(k1+k2)]))
  colnames(cred_A22)=c("est","2.5%","25%","50%","75%","97.5%")
  return(list(cred_A11=cred_A11,cred_A12=cred_A12,cred_A21=cred_A21,cred_A22=cred_A22))
}



