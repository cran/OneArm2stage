##' Performs the one-sample log-rank test (OSLR) for the time-to-event data from two-stage Phase II clinical trials, assuming the failure time follows one of the four distributions: Weibull, Gamma, log-normal or log-logistic.\cr
##' This can be used for both unrestricted and restricted follow-up designs.
##'
##' @title Conduct Interim or Final Analyses Using One-Sample Log-Rank Test for the Optimal Two-Stage Trials
##' @param dist distribution options with 'WB' as Weibull, 'GM' as Gamma, 'LN' as log-normal, 'LG' as log-logistic.
##' @param shape shape parameter for one of the four parametric distributions ('WB', 'GM', 'LN' and 'LG').
##' @param S0 the survival probability at a fixed time point x0 under the null hypothesis.
##' @param x0 a fixed time point when the survival probability is S0 under null.
##' @param data the time-to-event data for either the interim or final analysis from a two-stage survival trial, contains 2 variables:\cr
##'             \emph{time} time period under observation before the time of interim analysis (for interim analysis) or during entire trial (for final analysis) for each patient.\cr
##'             \emph{status} status indicator of patients (event = 1, censored = 0).
##' @return  \emph{z} the OSLR test statistic for the interim or final analysis, depending on data used.\cr
##'          \emph{O} the observed number of events.\cr
##'          \emph{E} the expected number of events.
##'
##' @import utils
##' @export
##' @examples
##' dat<- read.csv(system.file("extdata", "kj1_final.csv", package = "OneArm2stage"))
##' LRT(dist="WB", shape=1, S0=0.62, x0=2, data=dat)
##'# O       E       Z
##'# 18.0000 16.3598 -0.4055

##'@references
##' Wu, J, Chen L, Wei J, Weiss H, Chauhan A. (2020). Two-stage phase II survival trial design. Pharmaceutical Statistics. 2020;19:214-229. https://doi.org/10.1002/pst.1983


LRT <- function(dist, shape, S0, x0, data) {

  if (dist == "WB") {
    scale0 = x0/(-log(S0))^(1/shape)
    s0 = function(u) 1 - pweibull(u, shape, scale0)
  }


  if (dist == "LN") {
    scale0 = log(x0) - shape * qnorm(1 - S0)
    s0 = function(u) 1 - plnorm(u, scale0, shape)
  }

  if (dist == "LG") {
    scale0 = x0/(1/S0 - 1)^(1/shape)
    s0 = function(u) 1/(1 + (u/scale0)^shape)
  }

  if (dist == "GM") {
    root0 = function(t) 1 - pgamma(x0, shape, t) - S0
    rate0 = uniroot(root0, c(0, 10))$root
    s0 = function(u) 1 - pgamma(u, shape, rate0)
  }

  H = function(u) -log(s0(u)) # Hazard function for null distribution s0

  Time = data$time
  Cens = data$status
  O = sum(Cens)
  E = sum(H(Time))
  Z = (E - O)/sqrt(E)  # OSLR test statistic
  ans = c(O = O, E = round(E, 4), Z = round(Z, 4))
  return(ans)
}


