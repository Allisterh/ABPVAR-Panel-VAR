library(dplyr)
library(ggplot2)
library(alfred)
library(matrixcalc)
library(clusterGeneration)
library(MASS)
library(Matrix)
index <- 1

#==================
npop <- 4

#Finland, Sweden, Denmark, Norway
var <- list()

#Finland
var[[1]] <- c("CLVMNACSCAB1GQFI","FINGDPDEFQISMEI","FINPFCEQDSMEI","FINEXPORTQDSMEI","FINGFCEQDSMEI",
"FINCPIALLMINMEI", "LFHUTTTTFIM647S", "RBFIBIS", "IRLTLT01FIM156N","IR3TIB01FIM156N",
"FINPROINDMISMEI","CCUSMA02FIM618N","FINPPDMMINMEI","FINPRMNVG01IXOBSAM","XTIMVA01FIM664S")

#Sweden
var[[2]] <- c("CLVMNACSCAB1GQSE","SWEGDPDEFQISMEI","SWEPFCEQDSMEI","SWEEXPORTQDSMEI","SWEGFCEQDSMEI",
"SWECPIALLMINMEI","LFHUTTTTSEM647S","RBSEBIS","IRLTLT01SEM156N","IR3TIB01SEM156N",
"SWEPROINDMISMEI","CCUSMA02SEM618N","SWEPPDMMINMEI","SWEPRMNVG01IXOBSAM","XTIMVA01SEM664S")

#Denmark
var[[3]] <- c("CLVMNACSCAB1GQDK","DNKGDPDEFQISMEI","DNKPFCEQDSMEI","DNKEXPORTQDSMEI","DNKGFCEQDSMEI",
"DNKCPIALLMINMEI","LFHUTTTTDKM647S","RBDKBIS","IRLTLT01DKM156N","IR3TIB01DKM156N",
"DNKPROINDMISMEI","CCUSMA02DKM618N","DNKPPDMMINMEI","DNKPRMNVG01IXOBSAM","XTIMVA01DKM664S")

#Norway
var[[4]] <- c("CLVMNACSCAB1GQNO","NORGDPDEFQISMEI","NORPFCEQDSMEI","NOREXPORTQDSMEI","NORGFCEQDSMEI",
              "NORCPIALLMINMEI", "LFHUTTTTNOM647S", "RBNOBIS", "IRLTLT01NOM156N","IR3TIB01NOM156N",
                "NORPROINDMISMEI","CCUSMA02NOM618N","NORPPDMMINMEI","NORPRMNVG01IXOBSAM","XTIMVA01NOM664S")


l <- length(var[[4]])


out <- NULL
for ( i in 1: npop)
{
  out[[i]] <- list()
  out[[i]] <- lapply(var[[i]],get_alfred_series,observation_start="2000-04-01",
                    observation_end="2021-09-01",realtime_start="2021-12-10",realtime_end="2021-12-10")#Q2 1959 to Q3 2019
  
}

out[[2]][[10]][,3][20]=out[[2]][[10]][,3][19]  #Remove NA

###################################################################################

new <- NULL
for(i in 1: npop){

new[[i]] <- list()
set1 <- c(1,3,4,5)
for(j in set1){
  temp <- as.numeric(out[[i]][[j]][,3])
  diff <- diff(log(temp))
  new[[i]][[j]] <- as.data.frame(diff[-c(1:2)])
}

set2 <- c(2)
for(j in set2){
  temp <- as.numeric(out[[i]][[j]][,3])
  diff <- diff(log(temp),diff=2)
  new[[i]][[j]] <- as.data.frame(diff[-1])
}

set3 <- c(6,13)
for(j in set3){
  temp <- as.numeric((out[[i]][[j]][,3]))
  diff <-  diff(log(temp),diff=2)
  new[[i]][[j]] <- as.data.frame(diff[-c(1:7)])
}
set4 <- c(7,11,14,15)
for(j in set4){
  temp <- as.numeric(out[[i]][[j]][,3])
  diff <- diff(log(temp))
  new[[i]][[j]] <- as.data.frame(diff[-c(1:8)])
}

set5 <- c(8:10,12)
for(j in set5){
  temp <- as.numeric((out[[i]][[j]][,3]))
  diff <- diff(temp)
  new[[i]][[j]] <- as.data.frame(diff[-c(1:8)])
}
}


k1 <- 10
k2 <- 5
new1 <- NULL
for(i in 1:npop){
new1[[i]] <- list()
for(j in 1:k1){
  new1[[i]][[j]] <- new[[i]][[k2+j]]  
}
for(j in (k1+1):(k1+k2)){
  new1[[i]][[j]] <- new[[i]][[j-k1]]  
}
}

alfred_to_ts <- function(x,freq){
  ts(x[,1],start=c(2001,1),frequency=freq)
}
#mapply(function(a,b) b(a), mf_list, list(logd,logd,logd,logd,logd,sum,diff,diff,logd,
#logd2,logd2,logd,logd,logd,logd,diff,diff,diff,diff,logd))
mf_list <- NULL
for(i in 1:npop){
mf_list[[i]] <- list()
for(j in 1:k1){
  mf_list[[i]][[j]] <- ts(new1[[i]][[j]][,1],start=c(2001,1),frequency=12) 
}
for(j in (k1+1):(k1+k2) ){
  mf_list[[i]][[j]] <- ts(new1[[i]][[j]][,1],start=c(2001,1),frequency=4) 
}
}


s1 <- NULL
for(i in 1:length(mf_list[[1]][[(k1+1)]])){
  m <- c((3*(i-1))+1)  
  s1 <-c(s1,m)
}

s2 <- NULL
for(i in 1:length(mf_list[[1]][[(k1+1)]])){
  m <- c((3*(i-1))+2)  
  s2 <-c(s2,m)
}
s3 <- NULL
for(i in 1:length(mf_list[[1]][[(k1+1)]])){
  m <- c((3*(i-1))+3)  
  s3 <-c(s3,m)
}
hf <- list()
for(i in 1:npop){
  hf1 <- NULL
  hf2 <- NULL
  hf3 <- NULL
  
for(j in 1:k1){
  hf3 <- rbind(hf3, mf_list[[i]][[j]][s3])
  hf1 <- rbind(hf1, mf_list[[i]][[j]][s1])
  hf2 <- rbind(hf2, mf_list[[i]][[j]][s2])
}
  hf[[i]] <- rbind(hf3,hf1,hf2)
}


lf <- list()

for(i in 1:npop){
  lf1 <- NULL
for(j in (k1+1):(k1+k2)){
  lf1 <- rbind(lf1, mf_list[[i]][[j]])
}
  lf[[i]] <- lf1
}

data_full <- list()
for(i in 1:npop){
  data_full[[i]] <- rbind(hf[[i]],lf[[i]])# This is the full data containing train and test data
}
for(i in 1:npop){
  data_full[[i]] <- data_full[[i]][,1:82]  
}
