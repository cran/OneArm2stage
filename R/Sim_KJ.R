##' Sim_KJ() can be used to calculate empirical power and type-I error by simulation
##' given the design parameters (e.g., n1, n, c1, c) obtained from the optimal two-stage design with
##' unrestricted follow-up.
##' @title Calculate Empirical Power by Simulation for the Optimal Two-Stage Design Using
##'        One-Sample Log-Rank Test with Unrestricted Follow-Up
##' @param dist distribution options with 'WB' as Weibull, 'GM' as Gamma,
##' 'LN' as log-normal, 'LG' as log-logistic for the baseline hazard function.
##' @param shape shape parameter of the baseline hazard function assuming one of the four possible parametric distributions ('WB', 'GM',
##'  'LN' and 'LG').
##' @param S0 survival probability at a fixed time point x0 under the null hypothesis.
##' @param S1 survival probability at a fixed time point x0 under the alternative hypothesis.
##' @param x0 a fixed time point where the survival probabilities are known for both
##'        null and alternative hypotheses.
##' @param tf unrestricted follow-up time, the time period from the entry of the last patient to the end of the trial.
##' @param rate a constant accrual rate.
##' @param t1 the interim analysis time as given by the two-stage design.
##' @param c1 the critical value given by the two-stage design for the interim analysis.
##' @param c the critical value given by the two-stage design for the final analysis.
##' @param n1 the required sample size given by the two-stage design for the interim analysis.
##' @param n the required sample size given by the two-stage design for the final stage.
##' @param N number of trials in the simulation.
##' @param seed seed for random number generation.
##'
##' @return a numeric value that is either the empirical power (when S1=S0^hr) or the type I-error (when S1=S0).
##'
##' @export
##'
##' @examples
##' \donttest{Design2 <- Optimal.KJ(dist="WB", shape=1, S0=0.62, x0=2, hr=0.467, tf=2, rate=5,
##' alpha=0.05, beta=0.2)}
##' # ##' # $Two_stage
##' #   n1     c1  n      c     t1 MTSL      ES     PS
##' # 1 16 -0.302 26 1.6135 3.0593  7.2 21.9187 0.3813
##'
##' # calculate empirical power and type I error given the above design parameters
##' Sim_KJ(dist="WB", shape=1, S0=0.62, S1=0.62^(0.467), x0=2, tf=2, rate=5, t1=3.0593,
##'     c1=-0.302, c=1.6135, n1=16, n=26, N=10000, seed=5868)
##' # empirical power
##' # 0.813
##'
##' Sim_KJ(dist="WB", shape=1, S0=0.62, S1=0.62, x0=2, tf=2, rate=5, t1=3.0593,
##'     c1=-0.302, c=1.6135, n1=16, n=26, N=10000, seed=5868)
##' # empirical type-I error
##' # 0.037


Sim_KJ <- function(dist, shape, S0, S1, x0, tf, rate, t1, c1, c, n1, n, N, seed=123)
{
  ta=n/rate # accrual time
  MTSL=ta+tf

  # survival function based on distribution, for calculating expected events under Null
  if (dist == "WB") {
    scale0=x0/(-log(S0))^(1/shape) # scale0 is b, shape is a, in the paper Wu 2020
    scale=x0/(-log(S1))^(1/shape) # WB under proportional hazard assumption is still WB
    s=function(a,b,u){1-pweibull(u,a,b)} # u is quantile
  }

  if (dist == "LN") {
    inverse_hr=log(S0)/log(S1)
    scale0=log(x0)-shape*qnorm(1-S0)
    s=function(a,b,u){1-plnorm(u,b,a)} # u, q; meanlog, b(scale); sdlog, a(shape)
  }

  if (dist == "LG") {
    inverse_hr=log(S0)/log(S1)
    scale0=x0/(1/S0-1)^(1/shape)
    s=function(a,b,u){1/(1+(u/b)^a)}
  }

  if (dist == "GM") {
    inverse_hr=log(S0)/log(S1)
    root0=function(t){1-pgamma(x0, shape, t)-S0}
    ### solution t is rate not scale; scale=1/rate
    rate0=uniroot(root0,c(0,10))$root
    # we use rate0 to calculate H but change its name to scale0
    # the true scale is the inverse of rate
    scale0=rate0
    s=function(a,b,u){1-pgamma(u,a,b)} ## shape=a; rate=b
  }

  H=function(a,b,u){-log(s(a,b,u))}

  set.seed(seed)
  k=N1=0
  for (i in 1:N)
  {
    #entry time, unordered
    u=runif(n, 0, ta)

    # Generate observed failure time w
    if (dist == "WB") {
      w=rweibull(n, shape, scale)
    }

    # inverse cdf method
    if (dist == "LN") {
      U=runif(n, 0, 1)
      w=exp(scale0+shape*qnorm(1-U^inverse_hr))
    }

    if (dist == "LG") {
      U=runif(n, 0, 1)
      w=scale0*(1/U^inverse_hr-1)^(1/shape)
    }

    if (dist == "GM") {
      U=runif(n, 0, 1)
      w=(1/scale0)*qgamma(1-U^inverse_hr, shape, rate=1) # 1/scale0 is the true scale0
    }

    # interim analysis
    w1=w[1:n1]
    u1=u[1:n1]

    # if t1-u1 >0, observed time xt1 is the smaller btw entry to interim t1-u1 and survival time w1;
    # if t1-u1 <= 0, patient enrolled after or at t1, xt1=0
    xt1=pmax(0, pmin(w1, t1-u1))

    delta1 = as.numeric(w1 < t1-u1)  # death indicator
    O1=sum(delta1)                          # observed death events
    M1=H(shape,scale0,xt1)                  # H(a,b,u) under Null assumption
    E1=sum(M1)                              # expected death events
    Z1=(E1-O1)/E1^(1/2)

    # final analysis
    xt2=pmax(0, pmin(w, MTSL-u))
    delta2 = as.numeric(w < MTSL-u)

    O2=sum(delta2)
    M2=H(shape,scale0,xt2)
    E2=sum(M2)
    Z=(E2-O2)/E2^(1/2)

    if (E1!=0 & E2!=0)
    {
      # if(Z1 >c1)(k1=k1+1) # go at interim

      if (Z1>c1 & Z>c) (k=k+1) #success at final
      N1=N1+1


    }
  }
  ans=c(round(k/N1,3))
  return(ans)
}
