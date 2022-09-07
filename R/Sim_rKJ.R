##' Sim_rKJ() can be used to calculate empirical power and type-I error by simulation
##' given the design parameters (e.g., n1, n, c1, c) obtained from the optimal two-stage design with
##' restricted follow-up.
##' @title Calculate Empirical Power by Simulation for the Optimal Two-Stage Design Using
##'        One-Sample Log-Rank Test with Restricted Follow-Up
##' @param dist distribution options with 'WB' as Weibull, 'GM' as Gamma,
##' 'LN' as log-normal, 'LG' as log-logistic for the baseline hazard function.
##' @param shape shape parameter of the baseline hazard function assuming one of the four possible parametric distributions ('WB', 'GM',
##'  'LN' and 'LG').
##' @param S0 survival probability at a fixed time point x0 under the null hypothesis.
##' @param S1 survival probability at a fixed time point x0 under the alternative hypothesis.
##' @param x0 a fixed time point where the survival probabilities are known for both
##'        null and alternative hypotheses.
##' @param x the restricted follow-up time period.
##' @param rate a constant accrual rate.
##' @param t1 the interim analysis time in the two-stage design.\cr
##' @param c1 the critical value in two-stage designs for the interim analysis.
##' @param c the critical value in two-stage designs for the final analysis.\cr
##' @param n1 the required sample size in the two-stage design for the interim analysis.
##' @param n the required sample size in the two-stage design for the final stage.\cr
##' @param N number of trials in the simulation.
##' @param seed seed for random number generation.
##'
##' @return a numeric value that is either the empirical power (when S1=S0^hr) or the type I-error (when S1=S0).
##'
##' @export
##'
##' @examples
##' \donttest{Design3 <- Optimal.rKJ(dist="WB", shape=0.5, S0=0.3, x0=1, hr=0.65, x=1, rate=10,
##' alpha=0.05, beta=0.2)}
##' # Design3$Two_stage
##' # n1     c1      n    c      t1      MTSL   ES      PS
##' # 38     0.1688  63   1.6306 3.7084  7.3    48.3058 0.567
##'
##' # calculate empirical power and type I error given the above design parameters
##' Sim_rKJ(dist="WB", shape=0.5, S0=0.3, S1=0.3^(0.65), x0=1, x=1, rate=10, t1=3.7084,
##'     c1=0.1688, c=1.6306, n1=38, n=63, N=10000, seed=5868)
##' # empirical power
##' # 0.796
##'
##' Sim_rKJ(dist="WB", shape=0.5, S0=0.3, S1=0.3, x0=1, x=1, rate=10, t1=3.7084,
##'     c1=0.1688, c=1.6306, n1=38, n=63, N=10000, seed=5868)
##' # empirical type-I error
##' # 0.041


Sim_rKJ <- function(dist, shape, S0, S1, x0, x, rate, t1, c1, c, n1, n, N, seed=123)
{
  ta=n/rate # accrual time

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
    #entry time
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

    # if t1-u1 >0, observed time xt1 is the smallest among t1-u1, x, and w1;
    # if t1-u1 <= 0, patient enrolled after or at t1, xt1=0
    xt1=pmax(0, pmin(w1, pmin(x, t1-u1)))

    delta1 = as.numeric(w1 < pmin(x, t1-u1))  # status indicator
    O1=sum(delta1)                          # observed death events
    M1=H(shape,scale0,xt1)                  # H(a,b,u) under Null assumption
    E1=sum(M1)                              # expected death events
    Z1=(E1-O1)/E1^(1/2)

    # final analysis
    xt2=pmax(0, pmin(w, x)) # MTSL-u= ta+x-u > x
    delta2 = as.numeric(w < x)
    O2=sum(delta2)
    M2=H(shape,scale0,xt2)
    E2=sum(M2)
    Z=(E2-O2)/E2^(1/2)

    if (E1!=0 & E2!=0)
    {
      if (Z1>c1 & Z>c) (k=k+1)
      N1=N1+1

    }
  }
  ans=round(k/N1,3)
  return(ans)
}
