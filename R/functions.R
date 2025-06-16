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

truth_common_g2 <- function(k1,k2,target,num){
  p1 <- k1+k2
  tau <- runif(p1^2, 0.2, 0.9)*ifelse(rbinom(p1^2,1,0.7),-1,1)
  tau
  d_tau <- matrix(tau,nrow=p1,ncol=p1)
  d_tau
  sub <- d_tau[1:p1,1:num]
  for(i in 1:p1){
    sub[i,] <- rep(sub[i,1],num)
  }
  d_tau <- cbind(sub,d_tau[1:p1,(num+1):p1])
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


#=======================================Generate Sigma Matrix=================================
Gen_Sigma = function(k1,k2)
{
  set.seed(1233)
  c <- list()
  d <- list()
  p <- k1 + k2
  sigma_sq <- runif(p, 0.1, 1)
  ro <- rep(0.1,k1)
  sig_block1 <- diag(sigma_sq)
  for(i in 1:k1)
  {
    c[[i]] <- matrix(c(((ro[i]^2)*sigma_sq[i]),(ro[i]*sigma_sq[i])),nrow=2,ncol=1)
  }
  sig_block2 <- bdiag(c)
  null <- matrix(0,(2*k1),k2)
  sig_block2 <- cbind(sig_block2, null)
  sig_block3 <- t(sig_block2)
  for(i in 1:k1)
  {
    d[[i]] <- sigma_sq[i] * matrix(c(1,ro[i],ro[i],1), nrow=2)
  }
  sig_block4 <- bdiag(d)
  Sig1 <- rbind(sig_block1, sig_block2)
  Sig2 <- rbind(sig_block3, sig_block4)
  sigma_true <- cbind(Sig1,Sig2)
  #sigma_true
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
  y_true[[1]] <- mvrnorm(1,rep(0,p),sigma_true)
  y[[1]]<-y_true[[1]][vec]
  # y_mod[[1]]<- y[[1]][vec1]
  for ( t in 2:(nobs + horizon + 1))
  {
    err[[t]] <- mvrnorm(1,rep(0,p),sigma_true)
    y_true[[t]] <- as.vector((W %*% y_true[[t-1]]) + err[[t]])
    y[[t]]<-y_true[[t]][vec]
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

theta_matrix <- function(theta,b1,b2,b3,b4,k1,k2)
{
  theta_mat <- matrix(c(theta^2,theta^4,theta^3,1,(theta^2),theta,theta,theta^3,(theta^2)), nrow=3,ncol =3, byrow = T)
  u <- c(theta^2,1,theta)
  v <- c(1,theta^2,theta)
  A11 <- matrix(b1,nrow=k1,ncol=k1,byrow=F)
  A12 <- matrix(b2,nrow=k1,ncol=k2,byrow=F)
  A21 <- matrix(b3,nrow=k2,ncol=k1,byrow=F)
  A22 <- matrix(b4,nrow=k2,ncol=k2,byrow=F)
  w11 <- kronecker(theta_mat,A11)
  w12 <- kronecker(u,A12)
  w21 <- kronecker(t(v),A21)
  w22 <- A22
  w1 <- cbind(w11, w12)
  w2 <- cbind(w21, w22)
  w<- rbind(w1,w2)
  temp <- w[vec1,]
  w_mod <- temp[,vec1]
  return(list(theta_mat=theta_mat,u=u,v=v,w11=w11,w12=w12,w21=w21,w22=w22,w=w,w_mod=w_mod))
}


S_matrix <- function(data,n,k1,k2)
{
  s1 <- data[,1:n] %*% t(data[,1:n])
  s2 <- data[,2:(n+1)] %*% t(data[,1:n])
  s3 <- data[,2:(n+1)] %*% t(data[,2:(n+1)])
  
  s11 <- s1[1:(3*k1), 1:(3*k1)]
  s12 <- s1[1:(3*k1), ((3*k1)+1):((3*k1)+k2)]
  s21 <- s1[((3*k1)+1):((3*k1)+k2), 1:(3*k1)]
  s22 <- s1[((3*k1)+1):((3*k1)+k2), ((3*k1)+1):((3*k1)+k2)]
  
  ss11 <- s2[1:(3*k1), 1:(3*k1)]
  ss12 <- s2[1:(3*k1), ((3*k1)+1):((3*k1)+k2)]
  ss21 <- s2[((3*k1)+1):((3*k1)+k2), 1:(3*k1)]
  ss22 <- s2[((3*k1)+1):((3*k1)+k2), ((3*k1)+1):((3*k1)+k2)]
  
  s <- NULL
  for ( i in 1: 3)
  {
    s[[i]] <- list()
    
    for (j in 1: 3)
    {
      s[[i]][[j]] <- s11[(((i-1)*k1)+1):(i*k1),(((j-1)*k1)+1):(j*k1)]
    }
  }
  return(list(s1=s1,s2=s2,s3=s3,s11=s11,s12=s12,s21=s21,s22=s22,ss11=ss11,ss12=ss12,ss21=ss21,ss22=ss22,s=s))
}

sigma_inverse <- function(sigma,k1,k2)
{
  omega <- solve(sigma)
  sigma_11 <- omega[1:(3*k1),1:(3*k1)]
  sigma_22 <- omega[((3*k1)+1):((3*k1)+k2), ((3*k1)+1):((3*k1)+k2)]
  
  
  sigma_inv <- NULL
  for ( i in 1:3)
  {
    sigma_inv[[i]] <- list()
    for ( j in 1: 3)
    {
      sigma_inv[[i]][[j]] <- sigma_11[(((i-1)*k1)+1):(i*k1),(((j-1)*k1)+1):(j*k1)]
    }
  }
  return(list(sigma_11=sigma_11,sigma_22=sigma_22,sigma_inv=sigma_inv))
}



para_A11 <- function(k1,k2,theta_mat,s,sigma_inv,sigma_11,s11,ss11,s12,w12)
{
  gamma_A11 <- matrix(0,nrow=(k1*k1),ncol=(k1*k1))  
  for (i in 1:3)
  {
    for(j in 1:3)
    {
      mat1 <- matrix(0,nrow=k1,ncol=k1)
      mat2 <- matrix(0,nrow=k1,ncol=k1)
      for(k in 1:3)
      {
        mat1<- mat1+ s[[i]][[k]] * (t(theta_mat))[k,j]
        mat2 <- mat2 + sigma_inv[[j]][[k]] * theta_mat[k,i]
        
      }
      gamma1 <- kronecker(t(mat1),mat2)
      gamma_A11 <- gamma_A11 + gamma1  
    }
  }
  
  mean <- (ss11 %*% solve(s11)) - (w12 %*% t(s12) %*% solve(s11))
  
  P_mat <- sigma_11 %*% mean %*% s11
  
  P <- NULL
  for ( i in 1:3)
  {
    P[[i]] <- list()
    for ( j in 1: 3)
    {
      P[[i]][[j]] <- P_mat[(((i-1)*k1)+1):(i*k1),(((j-1)*k1)+1):(j*k1)]
    }
  }
  
  D <- matrix(0,nrow=k1,ncol=k1)  
  for (i in 1:3)
  {
    for(j in 1:3)
    {
      D <- D + (t(theta_mat))[i,j] * P[[j]][[i]]
      
    }
  }
  d_A11 <- vec(as.matrix(D))
  return(list(gamma_A11=gamma_A11,d_A11=d_A11))
}


para_A12 <- function(k1,k2,u,sigma_inv,sigma_11,s12,ss12,s22,w11)
{
  gamma <- matrix(0,nrow=k1,ncol=k1)  
  
  for(i in 1:3)
  {
    
    for(j in 1:3)
    {
      gamma <- gamma + ((u[i] * u[j]) * sigma_inv[[i]][[j]])
      
    }
  }
  gamma_A12 <- kronecker(t(s22),gamma)
  
  mean <- (ss12 %*% solve(s22)) -(w11 %*% s12 %*% solve(s22))
  
  P_mat <- sigma_11 %*% mean %*% s22
  P <- list()
  for(i in 1:3)
  {
    P[[i]] <- P_mat[(((i-1)*k1)+1):(i*k1),1:k2]
  }
  
  sum <- 0
  
  for(i in 1:3)
  {
    sum <- sum + u[i]*P[[i]] 
    
  }
  d_A12 <- vec(as.matrix(sum))
  return(list(gamma_A12=gamma_A12,d_A12=d_A12))
}

para_A21 <- function(k1,k2,v,sigma_22,s,s11,ss21,s12,w22)
{
  gamma_A21 <- matrix(0,nrow=(k1*k2),ncol=(k1*k2))  
  
  
  for(i in 1:3)
  {
    for(j in 1:3)
    {
      gamma_A21 <- gamma_A21 + kronecker((v[i] * t(s[[i]][[j]])), (sigma_22 * v[j]))
      
    }
  }
  
  mean <- (ss21 %*% solve(s11)) -(w22 %*% t(s12) %*% solve(s11))
  
  P_mat <- sigma_22 %*% mean %*% s11
  P <- list()
  for(i in 1:3)
  {
    P[[i]] <- P_mat[1:k2,(((i-1)*k1)+1):(i*k1)]
  }
  
  sum <- 0
  
  for(i in 1:3)
  {
    sum <- sum + (v[i]*P[[i]] )
    
  }
  d_A21 <- vec(as.matrix(sum))
  return(list(gamma_A21=gamma_A21,d_A21=d_A21))
}

para_A22 <- function(v,sigma_22,s22,ss22,s12,w21)
{
  
  mean <- (ss22 %*% solve(s22))-(w21 %*% s12 %*% solve(s22))
  m <- vec(mean)
  gamma_A22 <- solve(kronecker(solve(s22),solve(sigma_22)))
  d_A22 <- gamma_A22 %*% m
  return(list(gamma_A22=gamma_A22,d_A22=d_A22))
}

#-------------estimating coeff of pos dist of theta-----------------
theta_dist_coeff <- function(s1,s2,s3,sigma,b1,b2,b3,b4,k1,k2,vec1)
{
  Final_Sigma_Inv <- solve(sigma)
  theta_vec <- seq(0.1,1,length.out=20)
  state <- c(-4,-3,-2,-1,1,2,3,4)
  W <- list()
  x <- matrix(0,nrow=length(theta_vec),ncol=length(state))
  q1 <- numeric()
  d <- numeric()
  
  for(i in 1:length(theta_vec))
  {
    W[[i]]<- theta_matrix(theta_vec[i],b1,b2,b3,b4,k1,k2)[[8]] 
    q1[i] <- tr(as.matrix((Final_Sigma_Inv %*% s3) - (2*(t(W[[i]]) %*% Final_Sigma_Inv %*% s2)) + (t(W[[i]]) %*% Final_Sigma_Inv %*% W[[i]] %*% s1)  ))
    
    for (j in 1: length(state))
    {
      x[i,j]=(theta_vec[i])^state[j]
      
    }
  }
  
  d_theta <- as.vector(lm(q1 ~ x)$coef )
  d_theta <- d_theta[-1]
  return(list(d_theta=d_theta))
  
}

#---------------drawing theta--------------------

draw_theta <- function(N,d_theta)
{
  seq <- c(-4,-3,-2,-1,1,2,3,4)
  lnp <- numeric()
  x1 <- seq((1/(2*N)),((2*(N-1))/(2*N)),by=(1/N))
  for (i in 1:length(x1))
  {
    q1 <- (x1[i] ^ seq)
    lnp[i] <- sum(d_theta*q1)*(-1/2)
    
  }
  if (max(lnp) < -750 | max(lnp) > 700) lnp <- lnp + (700 - max(lnp))
  theta <- sample(x1,1,prob=exp(lnp))
  return(list(theta=theta))
  
}

data_initial <- function(j)
{
  A_true <- truth_gen(k1,k2,target,n1,n2,theta_true[j])$A
  W_true <- truth_gen(k1,k2,target,n1,n2,theta_true[j])$W
  index <- 1
  set.seed(index+j)
  data <- data_gen(W_true,k1,k2,horizon,sigma_true)$data
  X <- t(data)[1:nobs,]
  Y <- t(data)[2:(nobs + 1),]
  LSE <- as.matrix(lm(Y ~ X - 1)$coef)
  A11_initial <- t(LSE)[1:k1,1:k1]
  b1<- as.vector(vec(A11_initial))
  b2 <- numeric((k1*k2))
  b3 <- numeric((k1*k2))
  b4 <- numeric((k2*k2))
  A22_initial <- t(LSE)[(3*k1+1):(3*k1+k2),(3*k1+1):(3*k1+k2)]
  #b4<- as.vector(vec(A22_initial))
  theta_mat <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$theta_mat
  u <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$u
  v <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$v
  w11 <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$w11
  w12 <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$w12
  w21 <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$w21
  w22 <- theta_matrix(theta[j],b1,b2,b3,b4,k1,k2)$w22
  s <- S_matrix(data,nobs,k1,k2)$s
  s1 <- S_matrix(data,nobs,k1,k2)$s1
  s2 <- S_matrix(data,nobs,k1,k2)$s2
  s3 <- S_matrix(data,nobs,k1,k2)$s3
  s11 <- S_matrix(data,nobs,k1,k2)$s11
  s12 <- S_matrix(data,nobs,k1,k2)$s12
  s22 <- S_matrix(data,nobs,k1,k2)$s22
  ss11 <- S_matrix(data,nobs,k1,k2)$ss11
  ss12 <- S_matrix(data,nobs,k1,k2)$ss12
  ss21 <- S_matrix(data,nobs,k1,k2)$ss21
  ss22 <- S_matrix(data,nobs,k1,k2)$ss22
  gamma_A11 <- para_A11(k1,k2,theta_mat,s,sigma_inv,sigma_11,s11,ss11,s12,w12)$gamma_A11
  d_A11 <- para_A11(k1,k2,theta_mat,s,sigma_inv,sigma_11,s11,ss11,s12,w12)$d_A11
  gamma_A12 <- para_A12(k1,k2,u,sigma_inv,sigma_11,s12,ss12,s22,w11)$gamma_A12
  d_A12 <- para_A12(k1,k2,u,sigma_inv,sigma_11,s12,ss12,s22,w11)$d_A12
  gamma_A21 <- para_A21(k1,k2,v,sigma_22,s,s11,ss21,s12,w22)$gamma_A21
  d_A21 <- para_A21(k1,k2,v,sigma_22,s,s11,ss21,s12,w22)$d_A21
  gamma_A22 <- para_A22(v,sigma_22,s22,ss22,s12,w21)$gamma_A22
  d_A22 <- para_A22(v,sigma_22,s22,ss22,s12,w21)$d_A22
  return(list(W_true=W_true,gamma_A11=gamma_A11,d_A11=d_A11,
              gamma_A12=gamma_A12,d_A12=d_A12,
              gamma_A21=gamma_A21,d_A21=d_A21,
              gamma_A22=gamma_A22,d_A22=d_A22,s=s, s1=s1,s2=s2,s3=s3,s11=s11,s12=s12,s22=s22,ss11=ss11,
              ss12=ss12,ss21=ss21,ss22=ss22,A_true=A_true,data=data))
}


##############################################
##############################################

rep_dtheta <- function(j,s1_list,s2_list,s3_list,sigma,b1,b2,b3,b4,k1,k2,vec1)
{
  library(fBasics)
  #=========function "theta_dist_coeff"======================
  theta_dist_coeff <- function(s1,s2,s3,sigma,b1,b2,b3,b4,k1,k2,vec1)
  {
    Final_Sigma_Inv <- solve(sigma)
    theta_vec <- seq(0.1,1,length.out=20)
    state <- c(-4,-3,-2,-1,1,2,3,4)
    W <- list()
    x <- matrix(0,nrow=length(theta_vec),ncol=length(state))
    q1 <- numeric()
    d <- numeric()
    
    for(i in 1:length(theta_vec))
    {
      theta_mat <- matrix(c(theta_vec[i]^2,theta_vec[i]^4,theta_vec[i]^3,1,(theta_vec[i]^2),theta_vec[i],theta_vec[i],
                            theta_vec[i]^3,(theta_vec[i]^2)), nrow=3,ncol =3, byrow = T)
      u <- c(theta_vec[i]^2,1,theta_vec[i])
      v <- c(1,theta_vec[i]^2,theta_vec[i])
      A11 <- matrix(b1,nrow=k1,ncol=k1,byrow=F)
      A12 <- matrix(b2,nrow=k1,ncol=k2,byrow=F)
      A21 <- matrix(b3,nrow=k2,ncol=k1,byrow=F)
      A22 <- matrix(b4,nrow=k2,ncol=k2,byrow=F)
      w11 <- kronecker(theta_mat,A11)
      w12 <- kronecker(u,A12)
      w21 <- kronecker(t(v),A21)
      w22 <- A22
      w1 <- cbind(w11, w12)
      w2 <- cbind(w21, w22)
      W[[i]]<- rbind(w1,w2)
      q1[i] <- tr(as.matrix((Final_Sigma_Inv %*% s3) - (2*(t(W[[i]]) %*% Final_Sigma_Inv %*% s2)) + (t(W[[i]]) %*% Final_Sigma_Inv %*% W[[i]] %*% s1)  ))
      
      for (j in 1: length(state))
      {
        x[i,j]=(theta_vec[i])^state[j]
        
      }
    }
    
    d_theta <- as.vector(lm(q1 ~ x)$coef )
    d_theta <- d_theta[-1]
    return(list(d_theta=d_theta))
    
  }
  
  #==================== main code ==========================
  d_theta <- theta_dist_coeff(s1_list[[j]],s2_list[[j]],s3_list[[j]],sigma,b1[,j],b2[,j],b3[,j],b4[,j],k1,k2,vec1)$d_theta
  d_theta
}


#============================================================
#############################################################


rep_theta <- function(j,N,d_theta)
{
  #=========function "draw_theta"======================
  
  draw_theta <- function(N,d_theta)
  {
    seq <- c(-4,-3,-2,-1,1,2,3,4)
    lnp <- numeric()
    x1 <- seq((1/(2*N)),((2*(N-1))/(2*N)),by=(1/N))
    for (i in 1:length(x1))
    {
      q1 <- (x1[i] ^ seq)
      lnp[i] <- sum(d_theta*q1)*(-1/2)
      
    }
    if (max(lnp) < -750 | max(lnp) > 700) lnp <- lnp + (700 - max(lnp))
    theta <- sample(x1,1,prob=exp(lnp))
    return(list(theta=theta))
    
  }
  theta <- draw_theta(N,d_theta[,j])$theta
  theta
}

#============================================================
#############################################################

rep_thetamat <- function(j,theta,b1,b2,b3,b4,k1,k2,vec1)
{
  #=========function "theta_matrix"======================
  theta_matrix <- function(theta,b1,b2,b3,b4,k1,k2,vec1)
  {
    theta_mat <- matrix(c(theta^2,theta^4,theta^3,1,(theta^2),theta,theta,theta^3,(theta^2)), nrow=3,ncol =3, byrow = T)
    u <- c(theta^2,1,theta)
    v <- c(1,theta^2,theta)
    A11 <- matrix(b1,nrow=k1,ncol=k1,byrow=F)
    A12 <- matrix(b2,nrow=k1,ncol=k2,byrow=F)
    A21 <- matrix(b3,nrow=k2,ncol=k1,byrow=F)
    A22 <- matrix(b4,nrow=k2,ncol=k2,byrow=F)
    w11 <- kronecker(theta_mat,A11)
    w12 <- kronecker(u,A12)
    w21 <- kronecker(t(v),A21)
    w22 <- A22
    w1 <- cbind(w11, w12)
    w2 <- cbind(w21, w22)
    w<- rbind(w1,w2)
    temp <- w[vec1,]
    w_mod <- temp[,vec1]
    return(list(theta_mat=theta_mat,u=u,v=v,w11=w11,w12=w12,w21=w21,w22=w22,w=w,w_mod=w_mod))
  }
  theta_mat <- theta_matrix(theta[j],b1[,j],b2[,j],b3[,j],b4[,j],k1,k2,vec1)$theta_mat
  u <- theta_matrix(theta[j],b1[,j],b2[,j],b3[,j],b4[,j],k1,k2,vec1)$u
  v <- theta_matrix(theta[j],b1[,j],b2[,j],b3[,j],b4[,j],k1,k2,vec1)$v
  w11 <- theta_matrix(theta[j],b1[,j],b2[,j],b3[,j],b4[,j],k1,k2,vec1)$w11
  w12 <- theta_matrix(theta[j],b1[,j],b2[,j],b3[,j],b4[,j],k1,k2,vec1)$w12
  w21 <- theta_matrix(theta[j],b1[,j],b2[,j],b3[,j],b4[,j],k1,k2,vec1)$w21
  w22 <- theta_matrix(theta[j],b1[,j],b2[,j],b3[,j],b4[,j],k1,k2,vec1)$w22
  w <- theta_matrix(theta[j],b1[,j],b2[,j],b3[,j],b4[,j],k1,k2,vec1)$w
  return(list(theta_mat=theta_mat,u=u,v=v,w11=w11,w12=w12,w21=w21,w22=w22,w=w)) 
  
}

#=============================================================
##############################################################

rep_gamma <- function(j,theta_mat,s_list,sigma_inv,sigma_11,s11_list,ss11_list,s12_list,
                      w12,u,ss12_list,w11,v,sigma_22,ss21_list,ss22_list,s22_list,w22,w21,k1,k2)
{
  
  library(fBasics)
  #=========function "para_A11"======================
  para_A11 <- function(k1,k2,theta_mat,s,sigma_inv,sigma_11,s11,ss11,s12,w12)
  {
    gamma_A11 <- matrix(0,nrow=(k1*k1),ncol=(k1*k1))  
    for (i in 1:3)
    {
      for(j in 1:3)
      {
        mat1 <- matrix(0,nrow=k1,ncol=k1)
        mat2 <- matrix(0,nrow=k1,ncol=k1)
        for(k in 1:3)
        {
          mat1<- mat1+ s[[i]][[k]] * (t(theta_mat))[k,j]
          mat2 <- mat2 + sigma_inv[[j]][[k]] * theta_mat[k,i]
          
        }
        gamma1 <- kronecker(t(mat1),mat2)
        gamma_A11 <- gamma_A11 + gamma1  
      }
    }
    
    mean <- (ss11 %*% solve(s11)) - (w12 %*% t(s12) %*% solve(s11))
    
    P_mat <- sigma_11 %*% mean %*% s11
    
    P <- NULL
    for ( i in 1:3)
    {
      P[[i]] <- list()
      for ( j in 1: 3)
      {
        P[[i]][[j]] <- P_mat[(((i-1)*k1)+1):(i*k1),(((j-1)*k1)+1):(j*k1)]
      }
    }
    
    D <- matrix(0,nrow=k1,ncol=k1)  
    for (i in 1:3)
    {
      for(j in 1:3)
      {
        D <- D + (t(theta_mat))[i,j] * P[[j]][[i]]
        
      }
    }
    d_A11 <- vec(as.matrix(D))
    return(list(gamma_A11=gamma_A11,d_A11=d_A11))
  }
  
  #=========function "para_A12"======================
  para_A12 <- function(k1,k2,u,sigma_inv,sigma_11,s12,ss12,s22,w11)
  {
    gamma <- matrix(0,nrow=k1,ncol=k1)  
    
    for(i in 1:3)
    {
      
      for(j in 1:3)
      {
        gamma <- gamma + ((u[i] * u[j]) * sigma_inv[[i]][[j]])
        
      }
    }
    gamma_A12 <- kronecker(t(s22),gamma)
    
    mean <- (ss12 %*% solve(s22)) -(w11 %*% s12 %*% solve(s22))
    
    P_mat <- sigma_11 %*% mean %*% s22
    P <- list()
    for(i in 1:3)
    {
      P[[i]] <- P_mat[(((i-1)*k1)+1):(i*k1),1:k2]
    }
    
    sum <- 0
    
    for(i in 1:3)
    {
      sum <- sum + u[i]*P[[i]] 
      
    }
    d_A12 <- vec(as.matrix(sum))
    return(list(gamma_A12=gamma_A12,d_A12=d_A12))
  }
  
  #=========function "para_A21"======================
  para_A21 <- function(k1,k2,v,sigma_22,s,s11,ss21,s12,w22)
  {
    gamma_A21 <- matrix(0,nrow=(k1*k2),ncol=(k1*k2))  
    
    
    for(i in 1:3)
    {
      for(j in 1:3)
      {
        gamma_A21 <- gamma_A21 + kronecker((v[i] * t(s[[i]][[j]])), (sigma_22 * v[j]))
        
      }
    }
    
    mean <- (ss21 %*% solve(s11)) -(w22 %*% t(s12) %*% solve(s11))
    
    P_mat <- sigma_22 %*% mean %*% s11
    P <- list()
    for(i in 1:3)
    {
      P[[i]] <- P_mat[1:k2,(((i-1)*k1)+1):(i*k1)]
    }
    
    sum <- 0
    
    for(i in 1:3)
    {
      sum <- sum + (v[i]*P[[i]] )
      
    }
    d_A21 <- vec(as.matrix(sum))
    return(list(gamma_A21=gamma_A21,d_A21=d_A21))
  }
  
  #=========function "para_A22"======================
  para_A22 <- function(v,sigma_22,s22,ss22,s12,w21)
  {
    
    mean <- (ss22 %*% solve(s22))-(w21 %*% s12 %*% solve(s22))
    m <- vec(mean)
    gamma_A22 <- solve(kronecker(solve(s22),solve(sigma_22)))
    d_A22 <- gamma_A22 %*% m
    return(list(gamma_A22=gamma_A22,d_A22=d_A22))
  }
  
  
  gamma_A11 <- para_A11(k1,k2,theta_mat[[j]],s_list[[j]],sigma_inv,sigma_11,s11_list[[j]],ss11_list[[j]],
                        s12_list[[j]],w12[[j]])$gamma_A11
  d_A11 <- para_A11(k1,k2,theta_mat[[j]],s_list[[j]],sigma_inv,sigma_11,s11_list[[j]],ss11_list[[j]],
                    s12_list[[j]],w12[[j]])$d_A11
  
  gamma_A12 <- para_A12(k1,k2,u[[j]],sigma_inv,sigma_11,s12_list[[j]],ss12_list[[j]],s22_list[[j]],w11[[j]])$gamma_A12
  d_A12 <- para_A12(k1,k2,u[[j]],sigma_inv,sigma_11,s12_list[[j]],ss12_list[[j]],s22_list[[j]],w11[[j]])$d_A12
  
  gamma_A21 <- para_A21(k1,k2,v[[j]],sigma_22,s_list[[j]],s11_list[[j]],ss21_list[[j]],s12_list[[j]],w22[[j]])$gamma_A21
  d_A21 <- para_A21(k1,k2,v[[j]],sigma_22,s_list[[j]],s11_list[[j]],ss21_list[[j]],s12_list[[j]],w22[[j]])$d_A21
  
  gamma_A22 <- para_A22(v[[j]],sigma_22,s22_list[[j]],ss22_list[[j]],s12_list[[j]],w21[[j]])$gamma_A22
  d_A22 <- para_A22(v[[j]],sigma_22,s22_list[[j]],ss22_list[[j]],s12_list[[j]],w21[[j]])$d_A22
  
  return(list(gamma_A11=gamma_A11,d_A11=d_A11,
              gamma_A12=gamma_A12,d_A12=d_A12,
              gamma_A21=gamma_A21,d_A21=d_A21,
              gamma_A22=gamma_A22,d_A22=d_A22))
  
  
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
