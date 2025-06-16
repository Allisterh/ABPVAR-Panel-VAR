setwd("C:/Users/nchak/OneDrive/Desktop/R_codes")
library(dplyr)
library(ggplot2)
library(alfred)
library(matrixcalc)
library(clusterGeneration)
library(MASS)
library(Matrix)
library(glmnet)
#================== Set 1
npop <- 4
var <- list()

#Finland
var[[1]] <- c("CLVMNACSCAB1GQFI","FINGDPDEFQISMEI","FINPFCEQDSMEI","FINEXPORTQDSMEI","FINGFCEQDSMEI",
              "RBFIBIS","IR3TIB01FIM156N","CCUSMA02FIM618N",
              "FINCPIALLMINMEI", "LFHUTTTTFIM647S",  "IRLTLT01FIM156N",
"FINPROINDMISMEI","FINPPDMMINMEI","FINPRMNVG01IXOBSAM","XTIMVA01FIM664S")

#Sweden
var[[2]] <- c("CLVMNACSCAB1GQSE","SWEGDPDEFQISMEI","SWEPFCEQDSMEI","SWEEXPORTQDSMEI","SWEGFCEQDSMEI",
              "RBSEBIS","IR3TIB01SEM156N","CCUSMA02SEM618N",
              "SWECPIALLMINMEI","LFHUTTTTSEM647S","IRLTLT01SEM156N",
"SWEPROINDMISMEI","SWEPPDMMINMEI","SWEPRMNVG01IXOBSAM","XTIMVA01SEM664S")

#Denmark
var[[3]] <- c("CLVMNACSCAB1GQDK","DNKGDPDEFQISMEI","DNKPFCEQDSMEI","DNKEXPORTQDSMEI","DNKGFCEQDSMEI",
              "RBDKBIS","IR3TIB01DKM156N","CCUSMA02DKM618N",
"DNKCPIALLMINMEI","LFHUTTTTDKM647S","IRLTLT01DKM156N",
"DNKPROINDMISMEI","DNKPPDMMINMEI","DNKPRMNVG01IXOBSAM","XTIMVA01DKM664S")

#Norway
var[[4]] <- c("CLVMNACSCAB1GQNO","NORGDPDEFQISMEI","NORPFCEQDSMEI","NOREXPORTQDSMEI","NORGFCEQDSMEI",
              "RBNOBIS", "IR3TIB01NOM156N","CCUSMA02NOM618N",
              "NORCPIALLMINMEI", "LFHUTTTTNOM647S", "IRLTLT01NOM156N",
                "NORPROINDMISMEI","NORPPDMMINMEI","NORPRMNVG01IXOBSAM","XTIMVA01NOM664S")


out <- NULL
for ( i in 1: npop)
{
  out[[i]] <- list()
  out[[i]] <- lapply(var[[i]],get_alfred_series,observation_start="2000-04-01",
                    observation_end="2021-09-01",realtime_start="2021-12-10",realtime_end="2021-12-10")#Q2 1959 to Q3 2019
  
}

out[[2]][[7]][,3][20]=out[[2]][[7]][,3][19]  #Remove NA

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
  data_full[[i]] <- rbind(hf[[i]],lf[[i]])[,1:82]  
  }


save(data_full,file=paste("data_EU_g2_set1", ".dat", sep=''))

############### Set 2

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

l <- length(var[[3]])


out <- NULL
for ( i in 1: npop)
{
  out[[i]] <- list()
  out[[i]] <- lapply(var[[i]],get_alfred_series,observation_start="2000-04-01",
                     observation_end="2021-09-01",realtime_start="2022-02-01",realtime_end="2022-02-01")#Q2 1959 to Q3 2019
  
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
  data_full[[i]] <- rbind(hf[[i]],lf[[i]])[,1:82]  
}

save(data_full,file=paste("data_EU_g2_set2", ".dat", sep=''))

################# Set 3

npop <- 2
var <- list()

#ITALY
var[[1]] <- c("CLVMNACSCAB1GQIT","ITAGDPDEFQISMEI","ITAPFCEQDSMEI","ITAEXPORTQDSMEI","ITAGFCEQDSMEI",
              "RBITBIS","IR3TIB01ITM156N","CCUSMA02ITM618N",
              "ITACPIALLMINMEI", "LFHUTTTTITM647S", "IRLTLT01ITM156N",
              "ITAPROINDMISMEI","ITAPPDMMINMEI","ITAPROMANMISMEI","XTIMVA01ITM664S")

#SPAIN
var[[2]] <- c("CLVMNACSCAB1GQES","ESPGDPDEFQISMEI","ESPPFCEQDSMEI","ESPEXPORTQDSMEI","ESPGFCEQDSMEI",
              "RBESBIS", "IR3TIB01ESM156N","CCUSMA02ESM618N",
              "ESPCPIALLMINMEI", "LFHUTTTTESM647S", "IRLTLT01ESM156N",
              "ESPPROINDMISMEI","ESPPPDMMINMEI","ESPPROMANMISMEI","XTIMVA01ESM664S")

out <- NULL
for ( i in 1: npop)
{
  out[[i]] <- list()
  out[[i]] <- lapply(var[[i]],get_alfred_series,observation_start="2000-04-01",
                     observation_end="2021-09-01",realtime_start="2022-02-01",realtime_end="2022-02-01")#Q2 1959 to Q3 2019
  
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
  data_full[[i]] <- rbind(hf[[i]],lf[[i]])[,1:82]  
}

save(data_full,file=paste("data_EU_g2_set3", ".dat", sep=''))

######################## Set 4
npop <- 2
var <- list()

#Austria
var[[1]] <- c("CLVMNACSCAB1GQAT","AUTGDPDEFQISMEI","AUTPFCEQDSMEI","AUTEXPORTQDSMEI","AUTGFCEQDSMEI",
              "RBATBIS","IR3TIB01ATM156N", "CCUSMA02ATM618N",
              "AUTCPIALLMINMEI","LFHUTTTTATM647S","IRLTLT01ATM156N",
              "AUTPROINDMISMEI","AUTPPDMMINMEI","AUTPROMANMISMEI","XTIMVA01ATM664S")

#BELGIUM
var[[2]] <- c("CLVMNACSCAB1GQBE","BELGDPDEFQISMEI","BELPFCEQDSMEI","BELEXPORTQDSMEI","BELGFCEQDSMEI",
              "RBBEBIS", "IR3TIB01BEM156N","CCUSMA02BEM618N","BELCPIALLMINMEI", 
              "LFHUTTTTBEM647S","IRLTLT01BEM156N","BELPROINDMISMEI","BELPPDMMINMEI","BELPROMANMISMEI","XTIMVA01BEM664S")

l <- length(var[[2]])


out <- NULL
for ( i in 1: npop)
{
  out[[i]] <- list()
  out[[i]] <- lapply(var[[i]],get_alfred_series,observation_start="2000-04-01",
                     observation_end="2021-09-01",realtime_start="2022-02-01",realtime_end="2022-02-01")#Q2 1959 to Q3 2019
  
}

####################################################

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
  data_full[[i]] <- rbind(hf[[i]],lf[[i]])[,1:82]  
}

save(data_full,file=paste("data_EU_g2_set4", ".dat", sep=''))

