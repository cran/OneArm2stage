## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=F, warning=F----------------------------------------------
library(OneArm2stage)
library(survival)
library(IPDfromKM)

## ---- echo=T, eval=F----------------------------------------------------------
#  Design1 <- Optimal.rKJ(dist="WB", shape=1, S0=0.62, x0=2, hr=0.467, x=2, rate=5,
#                          alpha=0.05, beta=0.2)

## ----eval=F-------------------------------------------------------------------
#  Design1$Two_stage
#  #  n1     c1    n      c     t1    MTSL      ES     PS
#  # 1 26 0.1349  46   1.6171 5.0946  11.2   34.6342 0.5537

## ----echo=F-------------------------------------------------------------------
dat1 <- read.csv(system.file("extdata", "rkj1_interim.csv", package = "OneArm2stage"))

#order by entry
dat1[order(dat1$Entry),] 


## ----echo=F-------------------------------------------------------------------
KM1 <- survfit(Surv(time, status) ~ 1, data = dat1, conf.type = "log-log")
plot(KM1, mark.time = FALSE, main = "Kaplan-Meier survival curve", xlab = "Time", ylab = "Survival probability")

## ----echo=F-------------------------------------------------------------------
summary(KM1, times=KM1$time)

## -----------------------------------------------------------------------------
LRT(dist = "WB", shape = 1, S0=0.62, x0=2, data=dat1)

## ----echo=F-------------------------------------------------------------------
dat2.1 <- read.csv(system.file("extdata", "rkj2_interim.csv", package = "OneArm2stage"))

# order by entry
dat2.1[order(dat2.1$Entry),]

## ----echo=F-------------------------------------------------------------------
KM2 <- survfit(Surv(time, status) ~ 1, data = dat2.1, conf.type = "log-log")
plot(KM2, mark.time = FALSE, main = "Kaplan-Meier survival curve", xlab = "Time", ylab = "Survival probability")

## ----echo=F-------------------------------------------------------------------
summary(KM2, times=KM2$time)

## -----------------------------------------------------------------------------
LRT(dist = "WB", shape = 1, S0=0.62, x0=2, data=dat2.1)

## ----echo=F-------------------------------------------------------------------
dat2.2 <- read.csv(system.file("extdata", "rkj2_final.csv", package = "OneArm2stage"))

# order by entry
dat2.2[order(dat2.2$Entry),]

## ----echo=F-------------------------------------------------------------------
KM3 <- survfit(Surv(time, status) ~ 1, data = dat2.2, conf.type = "log-log")
plot(KM2, mark.time = FALSE, main = "Kaplan-Meier survival curve", xlab = "Time", ylab = "Survival probability")

## ----echo=F-------------------------------------------------------------------
summary(KM3, times=KM3$time)

## -----------------------------------------------------------------------------
LRT(dist="WB", shape=1, S0=0.62, x0=2, data=dat2.2)

## ---- echo=T, eval=F----------------------------------------------------------
#  Design2 <- Optimal.KJ(dist="WB", shape=1, S0=0.62, x0=2, hr=0.467, tf=2, rate=5,
#                        alpha=0.05, beta=0.2)

## ----eval=F-------------------------------------------------------------------
#  Design2$Two_stage
#  #   n1     c1   n      c     t1    MTSL      ES     PS
#  # 1 16 -0.302  26  1.6135 3.0593   7.2    21.9187 0.3813

## ---- echo=T, eval=F----------------------------------------------------------
#  Design3 <- Optimal.rKJ(dist="WB", shape=0.5, S0=0.3, x0=1, hr=0.65, x=1, rate=10,
#         alpha=0.05, beta=0.2)

## ----eval=F-------------------------------------------------------------------
#  Design3$Two_stage
#  #   n1     c1  n      c     t1  MTSL      ES    PS
#  # 1 38 0.1688 63 1.6306 3.7084  7.3    48.3058 0.567

## -----------------------------------------------------------------------------
# calculate empirical power
Sim_rKJ(dist="WB", shape=0.5, S0=0.3, S1=0.3^(0.65), x0=1, x=1, rate=10, t1=3.7084,
     c1=0.1688, c=1.6306, n1=38, n=63, N=10000, seed=5868)


## ----echo=T, eval=F-----------------------------------------------------------
#  Design6 <- Optimal.rKJ(dist="LN", shape=0.5, S0=0.3, x0=1, hr=0.65, x=1, rate=10,
#                          alpha=0.05, beta=0.2)

## ----eval=F-------------------------------------------------------------------
#  Design6$Two_stage
#  #>   n1     c1  n      c     t1  MTSL      ES     PS
#  #> 1 39 0.1097 63 1.6281 3.8408  7.3    49.6291 0.5437

## -----------------------------------------------------------------------------
# calculate empirical power
Sim_rKJ(dist="LN", shape=0.5, S0=0.3, S1=0.3^(0.65), x0=1, x=1, rate=10, t1=3.8408, c1=0.1097,  c=1.6281, n1=39, n=63, N=10000, seed=20198)


## ----echo=T, eval=F-----------------------------------------------------------
#  Design7 <- Optimal.rKJ(dist="LG", shape=0.5, S0=0.3, x0=1, hr=0.65, x=1, rate=10,
#                          alpha=0.05, beta=0.2)
#  

## ----eval=F-------------------------------------------------------------------
#  Design7$Two_stage
#  #   n1     c1  n      c     t1  MTSL     ES     PS
#  # 1 38 0.1958 63 1.6306 3.7089  7.3   48.034 0.5776

## -----------------------------------------------------------------------------
# calculate empirical power
Sim_rKJ(dist="LG", shape=0.5, S0=0.3, S1=0.3^(0.65), x0=1, x=1, rate=10, t1=3.7089, c1=0.1958,  c=1.6306, n1=38, n=63, N=10000, seed=20198)


## ----echo=T, eval=F-----------------------------------------------------------
#  Design8 <- Optimal.rKJ(dist="GM", shape=0.5, S0=0.3, x0=1, hr=0.65, x=1, rate=10,
#                          alpha=0.05, beta=0.2)

## ----eval=F-------------------------------------------------------------------
#  Design8$Two_stage
#  
#  #   n1     c1  n      c     t1  MTSL      ES   PS
#  # 1 38 0.1511 63 1.6306 3.7079  7.3   48.4841 0.56

## -----------------------------------------------------------------------------
# calculate empirical power
Sim_rKJ(dist="GM", shape=0.5, S0=0.3, S1=0.3^(0.65), x0=1, x=1, rate=10, t1=3.7079, c1=0.1511,  c=1.6306, n1=38, n=63, N=10000, seed=20198)


## ----eval=F-------------------------------------------------------------------
#  ## First, we take a screenshot of the plot and save it as a .png file. In this example it is named as "wang2018.png".
#  
#  ## Next, create a path directing to the .png file.
#  ## e.g. path <-"./vignettes/wang2018.png"
#  
#  ## Last, type the following code in console, and follow the prompts to extract data points by clicking your mouse on the plot.
#  ## df <- getpoints(path, 0, 24, 0, 100)  # 0,24,0,100 are the ranges of x and y-axes in the image
#  

## -----------------------------------------------------------------------------
# a sample dataset that we already extracted from Wang et al, 2018.
df<- read.csv(system.file("extdata", "df.csv", package = "OneArm2stage"))

# risk time points
trisk <- c(0,2,4,6,8,10,12,14,16,18,20,22,24) 

# number of patients at risk at each risk time point
nrisk.radio <- c(124,120,115,110,107,104,103,95,46,18,11,8,0) 

# Preprocess the raw coordinates into an proper format for reco""nstruct IPD
pre_radio <- preprocess(dat=df, trisk=trisk, nrisk=nrisk.radio,totalpts=NULL,maxy=100)

#Reconstruct IPD
est_radio <- getIPD(prep=pre_radio,armID=0,tot.events=NULL)

## -----------------------------------------------------------------------------
plot_ex5 <- survreport(ipd1=est_radio$IPD,ipd2=NULL,arms=1,interval=2,s=c(0.75,0.5,0.25),showplots=F)
plot_ex5$arm1$survprob

## -----------------------------------------------------------------------------
ipd <- est_radio$IPD
dat3 <- as.data.frame(cbind(rep(0, nrow(ipd)),ipd$time, ipd$status))
colnames(dat3) <- c("Entry", "Time", "Cens")

## -----------------------------------------------------------------------------
modelWB<- FitDat(dat3)

#AIC
modelWB$AIC

modelWB$parameter.estimates

## ----eval=F-------------------------------------------------------------------
#  Optimal.rKJ(dist="WB", shape=1, S0=0.62, x0=2, hr=0.467, x=2, rate=5,
#                    alpha=0.05, beta=0.2, prStop=0.65)
#  
#  # $Two_stage
#  #   n1    c1  n      c     t1 MTSL      ES     PS
#  # 1 30 0.386 46 1.6116 5.8247 11.2 35.0258 0.6503

## ----eval=F-------------------------------------------------------------------
#  # early stopping probability under H0
#  P <- t(seq(0.6, 0.85, l=6))
#  
#  # a design without pre-defined stopping probability
#  o1 <- Optimal.rKJ(dist="WB", shape=1, S0=0.62, x0=2, hr=0.467, x=2, rate=5,
#                    alpha=0.05, beta=0.2)
#  
#  # designs with varying pre-defined stopping probabilities
#  res1 <- t(apply(P, 2, Optimal.rKJ,
#                  dist="WB", shape=1, S0=0.62, x0=2, hr=0.467, x=2, rate=5, alpha=0.05, beta=0.2))
#  
#  # compare results
#  r1 <- rbind(o1$Two_stage,
#              res1[[1]]$Two_stage, res1[[2]]$Two_stage,
#              res1[[3]]$Two_stage, res1[[4]]$Two_stage,
#              res1[[5]]$Two_stage, res1[[6]]$Two_stage)
#  
#  # first row shows the design without pre-defined stopping probability (prStop=NULL)
#  # > r1
#  #   n1     c1  n      c     t1 MTSL      ES     PS
#  # 1 26 0.1349 46 1.6171 5.0946 11.2 34.6350 0.5537
#  # 2 28 0.2544 46 1.6147 5.4413 11.2 34.7161 0.6004
#  # 3 30 0.3860 46 1.6116 5.8247 11.2 35.0258 0.6503
#  # 4 31 0.5247 47 1.6013 6.1331 11.4 35.5642 0.7001
#  # 5 33 0.6760 47 1.5952 6.5867 11.4 36.4433 0.7505
#  # 6 36 0.8436 46 1.5927 7.1199 11.2 37.6740 0.8005
#  # 7 39 1.0382 46 1.5768 7.6639 11.2 39.4686 0.8504

