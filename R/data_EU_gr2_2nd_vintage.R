setwd("C:/Users/nchak/OneDrive/Desktop/R_codes")
library(dplyr)
library(ggplot2)
library(alfred)
library(matrixcalc)
library(clusterGeneration)
library(MASS)
library(Matrix)
library(glmnet)
#==================
npop <- 3
var <- list()

#Germany
var[[1]] <- c("CLVMNACSCAB1GQDE","DEUGDPDEFQISMEI","DEUPFCEQDSMEI","DEUEXPORTQDSMEI","DEUGFCEQDSMEI",
              "RBDEBIS","IR3TIB01DEM156N","CCUSMA02DEM618N",
              "DEUCPIALLMINMEI", "LFHUTTTTDEM647S",  "IRLTLT01DEM156N",
              "DEUPROINDMISMEI","DEUPPDMMINMEI","DEUPROMANMISMEI","XTIMVA01DEM664S")

#France
var[[2]] <- c("CLVMNACSCAB1GQFR","FRAGDPDEFQISMEI","FRAPFCEQDSMEI","FRAEXPORTQDSMEI","FRAGFCEQDSMEI",
              "RBFRBIS","IR3TIB01FRM156N","CCUSMA02FRM618N",
              "FRACPIALLMINMEI","LFHUTTTTFRM647S","IRLTLT01FRM156N",
              "FRAPROINDMISMEI","FRAPPDMMINMEI","FRAPROMANMISMEI","XTIMVA01FRM664S")

#Netherlands
var[[3]] <- c("CLVMNACSCAB1GQNL","NLDGDPDEFQISMEI","NLDPFCEQDSMEI","NLDEXPORTQDSMEI","NLDGFCEQDSMEI",
              "RBNLBIS","IR3TIB01NLM156N","CCUSMA02NLM618N",
              "NLDCPIALLMINMEI","LFHUTTTTNLM647S","IRLTLT01NLM156N",
              "NLDPROINDMISMEI","NLDPPDMMINMEI","NLDPROMANMISMEI","XTIMVA01NLM664S")

out <- NULL
for ( i in 1: npop)
{
  out[[i]] <- list()
  out[[i]] <- lapply(var[[i]],get_alfred_series,observation_start="2000-04-01",
                     observation_end="2020-06-01",
                     realtime_start="2020-10-09",realtime_end="2020-10-09"
  )
}

for(i in 1:npop){
  for(j in 1:length(var[[3]])){
    id=which(is.na(out[[i]][[j]][,3])==T) 
    for(k in id){
      out[[i]][[j]][,3][k]=out[[i]][[j]][,3][k-1]
    }   
  }
}  
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
  
  set3 <- c(9,13)
  for(j in set3){
    temp <- as.numeric((out[[i]][[j]][,3]))
    diff <-  diff(log(temp),diff=2)
    new[[i]][[j]] <- as.data.frame(diff[-c(1:7)])
  }
  set4 <- c(10,12,14,15)
  for(j in set4){
    temp <- as.numeric(out[[i]][[j]][,3])
    diff <- diff(log(temp))
    new[[i]][[j]] <- as.data.frame(diff[-c(1:8)])
  }
  set5 <- c(6,7,8,11)
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
  data_full[[i]] <- rbind(hf[[i]],lf[[i]])[,1:78]  
  }

save(data_full,file=paste("data_EU_g2_2nd_vintage", ".dat", sep=''))
