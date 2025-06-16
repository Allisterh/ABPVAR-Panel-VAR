#Data, using 2nd vintage, for analysis using ABPVAR with grouping 2
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
npop <- 5
var <- list()

#Minnesota
var[[1]] <- c("MNMFG","MNCONS","SMU27000002000000003","SMU27000003000000007","SMU27000002000000002",
              "LBSSA27","MNUR","MNNA","SMU27000005552410001SA","SMS27000004200000001","SMU27000000500000003SA",
              "SMU27000003000000003","SMU27000005500000003","SMU27000000500000002","SMU27000005500000002")

#Ohio
var[[2]] <- c("OHMFG","OHCONS","SMU39000002000000003","SMU39000003000000007","SMU39000002000000002",
              "LBSSA39","OHUR","OHNA","SMU39000005552410001SA","SMS39000004200000001","SMU39000000500000003SA",
              "SMU39000003000000003","SMU39000005500000003","SMU39000000500000002","SMU39000005500000002")

#Michigan
var[[3]] <- c("MIMFG","MICONS","SMU26000002000000003","SMU26000003000000007","SMU26000002000000002",
              "LBSSA26","MIUR","MINA","SMU26000005552400001SA","SMS26000004200000001","SMU26000000500000003SA",
              "SMU26000003000000003","SMU26000005500000003","SMU26000000500000002","SMU26000005500000002")

#Wisconsin
var[[4]] <- c("WIMFG","WICONS","SMU55000002000000003","SMU55000003000000007","SMU55000002000000002",
              "LBSSA55","WIUR","WINA","SMU55000005552400001SA","SMS55000004200000001","SMU55000000500000003SA",
              "SMU55000003000000003","SMU55000005500000003","SMU55000000500000002","SMU55000005500000002")
#Illinois
var[[5]] <- c("ILMFG","ILCONS","SMU17000002000000003","SMU17000003000000007","SMU17000002000000002",
              "LBSSA17","ILUR","ILNA","SMU17000005552400001SA","SMS17000004200000001","SMU17000000500000003SA",
              "SMU17000003000000003","SMU17000005500000003","SMU17000000500000002","SMU17000005500000002")

out <- NULL
for ( i in 1: npop)
{
  out[[i]] <- list()
  out[[i]] <- lapply(var[[i]],get_alfred_series,observation_start="2007-04-01",
                     observation_end="2020-04-01",realtime_start="2020-06-19",realtime_end="2020-06-19")#Q2 1959 to Q3 2019
  
}

###################################################################################

new <- NULL
for(i in 1: npop){

new[[i]] <- list()
set1 <- c(6,7)
for(j in set1){
  temp <- as.numeric((out[[i]][[j]][,3]))
  diff <- diff(temp)
  new[[i]][[j]] <- as.data.frame(diff)
}

set2 <- c(8, 1, 2, 9, 10, 11, 3, 12, 13)
for(j in set2){
  temp <- as.numeric(out[[i]][[j]][,3])
  diff <- diff(log(temp))
  new[[i]][[j]] <- as.data.frame(diff)
}


set3 <- c(4,14,5,15)
for(j in set3){
  temp <- as.numeric((out[[i]][[j]][,3]))
  new[[i]][[j]] <- as.data.frame(temp[-1])
}
}


k1 <- 15
alfred_to_ts <- function(x,freq){
  ts(x[,1],start=c(2001,1),frequency=freq)
}
mf_list <- NULL
for(i in 1:npop){
mf_list[[i]] <- list()
for(j in 1:k1){
  mf_list[[i]][[j]] <- ts(new[[i]][[j]][,1],start=c(2007,5),frequency=12) 
}
}

lf <- list()

for(i in 1:npop){
  lf1 <- NULL
for(j in 1:15){
  lf1 <- rbind(lf1, mf_list[[i]][[j]])
}
  lf[[i]] <- lf1
}

save(lf,file=paste("data_US_states_gr2_2nd", ".dat", sep=''))

  