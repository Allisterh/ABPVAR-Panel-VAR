#Data, using 2nd vintage, for analysis using ABPVAR with grouping 1
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
npop <- 7
var <- list()

#Minnesota
var[[1]] <- c("LBSSA27","MNUR","MNNA","MNMFG","MNCONS","SMU27000005552410001SA","SMS27000004200000001",
              "SMU27000000500000003SA","SMU27000002000000003","SMU27000003000000003",
              "SMU27000005500000003","SMU27000003000000007","SMU27000000500000002","SMU27000002000000002","SMU27000005500000002")

#Ohio
var[[2]] <- c("LBSSA39","OHUR","OHNA","OHMFG","OHCONS","SMU39000005552410001SA","SMS39000004200000001",
              "SMU39000000500000003SA","SMU39000002000000003","SMU39000003000000003",
              "SMU39000005500000003","SMU39000003000000007","SMU39000000500000002","SMU39000002000000002","SMU39000005500000002")

#Iowa
var[[3]] <- c("LBSSA19","IAUR","IANA","IAMFG","IACONS","SMU19000005552400001SA","SMS19000004200000001",
              "SMU19000000500000003SA","SMU19000002000000003","SMU19000003000000003",
              "SMU19000005500000003","SMU19000003000000007","SMU19000000500000002","SMU19000002000000002","SMU19000005500000002")

#Michigan
var[[4]] <- c("LBSSA26","MIUR","MINA","MIMFG","MICONS","SMU26000005552400001SA","SMS26000004200000001",
              "SMU26000000500000003SA","SMU26000002000000003","SMU26000003000000003",
              "SMU26000005500000003","SMU26000003000000007","SMU26000000500000002","SMU26000002000000002","SMU26000005500000002")

#Wisconsin
var[[5]] <- c("LBSSA55","WIUR","WINA","WIMFG","WICONS","SMU55000005552400001SA","SMS55000004200000001",
              "SMU55000000500000003SA","SMU55000002000000003","SMU55000003000000003",
              "SMU55000005500000003","SMU55000003000000007","SMU55000000500000002","SMU55000002000000002","SMU55000005500000002")
#Illinois
var[[6]] <- c("LBSSA17","ILUR","ILNA","ILMFG","ILCONS","SMU17000005552400001SA","SMS17000004200000001",
              "SMU17000000500000003SA","SMU17000002000000003","SMU17000003000000003",
              "SMU17000005500000003","SMU17000003000000007","SMU17000000500000002","SMU17000002000000002","SMU17000005500000002")


#Georgia
var[[7]] <- c("LBSSA13","GAUR","GANA","GAMFG","GACONS","SMU13000005552400001SA","SMS13000004200000001",
              "SMU13000000500000003SA","SMU13000002000000003","SMU13000003000000003",
              "SMU13000005500000003","SMU13000003000000007","SMU13000000500000002","SMU13000002000000002","SMU13000005500000002")


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
  set1 <- c(1,2)
  for(j in set1){
    temp <- as.numeric((out[[i]][[j]][,3]))
    diff <- diff(temp)
    new[[i]][[j]] <- as.data.frame(diff)
  }
  
  set2 <- c(3:11)
  for(j in set2){
    temp <- as.numeric(out[[i]][[j]][,3])
    diff <- diff(log(temp))
    new[[i]][[j]] <- as.data.frame(diff)
  }
  
  
  set3 <- c(12:15)
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
save(lf,file=paste("data_US_states_gr1_2nd", ".dat", sep=''))

