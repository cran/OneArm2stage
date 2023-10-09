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
#  Design1 <- phase2.TTE(shape=1,S0=0.72,x0=3,hr=0.459,tf=3,rate=20,alpha=0.05,
#  					 beta=0.2, restricted=0)

## ----eval=F-------------------------------------------------------------------
#  Design1
#  # $Two_stage_Optimal
#  #   n1      c1  n      c     t1 MTSL      ES     PS
#  # 1 22 -0.7215 46 1.5952 1.0619  5.3 40.1736 0.2353
#  #
#  # $Two_stage_minmax
#  #   n1     c1  n      c     t1 MTSL      ES     PS
#  # 1 34 -0.731 44 1.6257 1.6997  5.2 41.6745 0.2324
#  #
#  # $Two_stage_Admissible
#  #     n1      c1  n      c     t1 MTSL      ES     PS      Rho
#  # 88  26 -0.7274 45 1.6061 1.2512 5.25 40.3361 0.2335 42.66805
#  # 98  29 -1.0069 44 1.6293 1.4269 5.20 41.5723 0.1570 42.78615
#  # 108 21 -0.7393 46 1.5952 1.0326 5.30 40.1736 0.2298 43.08680
#  

## ---- echo=T, eval=F----------------------------------------------------------
#  Design2 <- phase2.TTE(shape=1,S0=0.72,x0=3,hr=0.459,tf=3,rate=20,alpha=0.05,
#  					 beta=0.2, prStop=0.35, restricted=0)

## ----eval=F-------------------------------------------------------------------
#  Design2
#  # $Two_stage_Optimal
#  #   n1      c1  n      c     t1 MTSL      ES     PS
#  # 1 29 -0.3844 47 1.5799 1.4364 5.35 40.5985 0.3503
#  #
#  # $Two_stage_minmax
#  #   n1      c1  n      c     t1 MTSL      ES     PS
#  # 1 40 -0.3829 44 1.6171 1.9808  5.2 42.4619 0.3509
#  #
#  # $Two_stage_Admissible
#  #     n1      c1  n      c     t1 MTSL      ES     PS      Rho
#  # 92  34 -0.3850 45 1.5988 1.6732 5.25 40.9611 0.3501 42.98055
#  # 102 40 -0.3846 44 1.6171 1.9791 5.20 42.4522 0.3503 43.22610
#  # 921 31 -0.3851 46 1.5884 1.5339 5.30 40.6365 0.3501 43.31825
#  # 91  29 -0.3850 47 1.5799 1.4353 5.35 40.5949 0.3501 43.79745
#  
#  

## ----eval=F-------------------------------------------------------------------
#  Design3 <- phase2.TTE(shape=1,S0=0.72,x0=3,hr=0.459,tf=3,rate=20,alpha=0.05,
#  					 beta=0.2, q_value=0.1, restricted=0)
#  
#  Design3_2 <- phase2.TTE(shape=1,S0=0.72,x0=3,hr=0.459,tf=3,rate=20,alpha=0.05,
#  					 beta=0.2, q_value=0.7, restricted=0)

## ----eval=F-------------------------------------------------------------------
#  
#  # q_value = 0.1
#  Design3
#  #
#  # $Two_stage_Admissible
#  #     n1      c1  n      c     t1 MTSL      ES     PS      Rho
#  # 108 21 -0.7393 46 1.5952 1.0326 5.30 40.1736 0.2298 40.75624
#  # 88  26 -0.7274 45 1.6061 1.2512 5.25 40.3361 0.2335 40.80249
#  # 98  29 -1.0069 44 1.6293 1.4269 5.20 41.5723 0.1570 41.81507
#  
#  # q_value = 0.7
#  Design3_2
#  # $Two_stage_Admissible
#  #     n1      c1  n      c     t1 MTSL      ES     PS      Rho
#  # 98  29 -1.0069 44 1.6293 1.4269 5.20 41.5723 0.1570 43.27169
#  # 88  26 -0.7274 45 1.6061 1.2512 5.25 40.3361 0.2335 43.60083
#  # 108 21 -0.7393 46 1.5952 1.0326 5.30 40.1736 0.2298 44.25208

## ----echo=F-------------------------------------------------------------------
dat1 <- read.csv(system.file("extdata", "d4_t1.csv", package = "OneArm2stage"))
dat1 <- dat1[order(dat1$Entry), c("Entry", "time", "status", "Total")] #order by entry


print(dat1)

## ----echo=F-------------------------------------------------------------------
KM1 <- survfit(Surv(time, status) ~ 1, data = dat1, conf.type = "log-log")
plot(KM1, mark.time = FALSE, main = "Kaplan-Meier survival curve", xlab = "Time", ylab = "Survival probability")

## ----echo=F-------------------------------------------------------------------
summary(KM1, times=KM1$time)

## -----------------------------------------------------------------------------
LRT(shape = 1, S0=0.72, x0=3, data=dat1)

## ----echo=F-------------------------------------------------------------------
dat2 <- read.csv(system.file("extdata", "d4.csv", package = "OneArm2stage"))
dat2 <- dat2[order(dat2$Entry), c("Entry", "time", "status", "Total")]
# order by entry
print(dat2)

## ----echo=F-------------------------------------------------------------------
KM2 <- survfit(Surv(time, status) ~ 1, data = dat2, conf.type = "log-log")
plot(KM2, mark.time = FALSE, main = "Kaplan-Meier survival curve", xlab = "Time", ylab = "Survival probability")

## ----echo=F-------------------------------------------------------------------
summary(KM2, times=KM2$time)

## -----------------------------------------------------------------------------
LRT(shape=1, S0=0.72, x0=3, data=dat2)

## ---- echo=T, eval=F----------------------------------------------------------
#  Design5 <- phase2.TTE(shape=1,S0=0.72,x0=3,hr=0.459,tf=3,rate=20,alpha=0.05,
#  					 beta=0.2, restricted=1)

## ----eval=F-------------------------------------------------------------------
#  Design5
#  # $Two_stage_Optimal
#  #   n1      c1  n   c     t1 MTSL      ES    PS
#  # 1 36 -0.3665 60 1.6 1.7523    6 51.0914 0.357
#  #
#  # $Two_stage_minmax
#  #   n1      c1  n      c     t1 MTSL      ES     PS
#  # 1 37 -0.5517 58 1.6193 1.8222  5.9 51.7359 0.2906
#  #
#  # $Two_stage_Admissible
#  #     n1      c1  n      c     t1 MTSL      ES     PS      Rho
#  # 120 37 -0.5625 58 1.6196 1.8066 5.90 51.7260 0.2869 54.86300
#  # 88  36 -0.4531 59 1.6086 1.7515 5.95 51.2044 0.3252 55.10220
#  # 232 36 -0.3643 60 1.6000 1.7556 6.00 51.0949 0.3578 55.54745
#  

## ----eval=F-------------------------------------------------------------------
#  Design1$param
#  # shape   S0    hr alpha beta rate x0 tf q_value prStop restricted
#  #     1 0.72 0.459  0.05  0.2   20  3  3     0.5      0          0
#  
#  Design1$Two_stage_Optimal
#  #   n1      c1  n      c     t1 MTSL      ES     PS
#  # 1 22 -0.7215 46 1.5952 1.0619  5.3 40.1736 0.2353
#  
#  Design1$Two_stage_minmax
#  #   n1     c1  n      c     t1 MTSL      ES     PS
#  # 1 34 -0.731 44 1.6257 1.6997  5.2 41.6745 0.2324
#  
#  Design1$Two_stage_Admissible
#  #     n1      c1  n      c     t1 MTSL      ES     PS      Rho
#  # 88  26 -0.7274 45 1.6061 1.2512 5.25 40.3361 0.2335 42.66805
#  # 98  29 -1.0069 44 1.6293 1.4269 5.20 41.5723 0.1570 42.78615
#  # 108 21 -0.7393 46 1.5952 1.0326 5.30 40.1736 0.2298 43.08680
#  

## -----------------------------------------------------------------------------
# calculate empirical power based on optimal method, power=0.778
Sim(shape=1, S0=0.72, S1=0.72^(0.459), x0=3, tf=3, rate=20, t1=1.0619,
     c1=-0.7215, c=1.5952, n1=22, n=46, N=10000, seed=5868)

# calculate empirical power based on minmax method, power=0.803
Sim(shape=1, S0=0.72, S1=0.72^(0.459), x0=3, tf=3, rate=20, t1=1.6997,
     c1=-0.731, c=1.6257, n1=34, n=44, N=10000, seed=5868)

# calculate empirical power based on admissible method, power=0.789
Sim(shape=1, S0=0.72, S1=0.72^(0.459), x0=3, tf=3, rate=20, t1=1.2512,
     c1=-0.7274, c=1.6061, n1=26, n=45, N=10000, seed=5868)


## ----eval=F-------------------------------------------------------------------
#  Design5
#  # $Two_stage_Optimal
#  #   n1      c1  n   c     t1 MTSL      ES    PS
#  # 1 36 -0.3665 60 1.6 1.7523    6 51.0914 0.357
#  #
#  # $Two_stage_minmax
#  #   n1      c1  n      c     t1 MTSL      ES     PS
#  # 1 37 -0.5517 58 1.6193 1.8222  5.9 51.7359 0.2906
#  #
#  # $Two_stage_Admissible
#  #     n1      c1  n      c     t1 MTSL      ES     PS      Rho
#  # 120 37 -0.5625 58 1.6196 1.8066 5.90 51.7260 0.2869 54.86300
#  # 88  36 -0.4531 59 1.6086 1.7515 5.95 51.2044 0.3252 55.10220
#  # 232 36 -0.3643 60 1.6000 1.7556 6.00 51.0949 0.3578 55.54745
#  

## -----------------------------------------------------------------------------
# calculate empirical power based on optimal method, power=0.799
Sim(shape=1, S0=0.72, S1=0.72^(0.459), x0=3, tf=3, rate=20, t1=1.7523,
     c1=-0.3665, c=1.6, n1=36, n=60, N=10000, restricted=1, seed=5868)

# calculate empirical power based on minmax method, power=0.796
Sim(shape=1, S0=0.72, S1=0.72^(0.459), x0=3, tf=3, rate=20, t1=1.8222,
     c1=-0.5517, c=1.6193, n1=37, n=58, N=10000, restricted=1, seed=5868)

# calculate empirical power based on admissible method, power=0.796
Sim(shape=1, S0=0.72, S1=0.72^(0.459), x0=3, tf=3, rate=20, t1=1.8066,
     c1=-0.5625, c=1.6196, n1=37, n=58, N=10000, restricted=1, seed=5868)


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

