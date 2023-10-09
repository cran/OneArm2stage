##' Sim() can be used to calculate the empirical power and type-I error by
##' simulation given the design parameters obtained from the two-stage
##' designs using phase2.TTE().
##' @title Calculate Empirical Power by Simulation for Two-Stage Designs Using
##'        One-Sample Log-Rank Test
##' @param shape shape parameter of the baseline hazard function assuming
##' Weibull distribution.
##' @param S0 survival probability at a fixed time point x0 under the null
##' hypothesis.
##' @param S1 survival probability at a fixed time point x0 under the
##' alternative hypothesis.
##' @param x0 a fixed time point where the survival probabilities are known for
##' both null and alternative hypotheses.
##' @param tf the follow-up time, the time period from the entry of the last
##' patient to the end of the trial.
##' @param rate a constant accrual rate.
##' @param t1 the interim analysis time as given by the two-stage design.
##' @param c1 the critical value for the interim analysis as given by the
##' two-stage design.
##' @param c the critical value for the final analysis as given by the two-stage
##' design.
##' @param n1 the required sample size for the interim analysis as given by the
##' two-stage design.
##' @param n the required total sample size for the entire trial as given by the
##' two-stage design.
##' @param N number of trials in the simulation.
##' @param restricted if restricted = 0 (default value), unrestricted follow-up
##' is used. If restricted = 1 or any other values, restricted follow-up is used.
##' @param seed seed for random number generation.
##'
##' @return a numeric value that is either the empirical power (when S1=S0^hr)
##' or the type I-error (when S1=S0).
##'
##' @export
##'
##' @examples
##' # calculate empirical power and type I error given design parameters
##' Sim(shape=1, S0=0.62, S1=0.62^(0.467), x0=2, tf=2, rate=5, t1=3.0593,
##'     c1=-0.302, c=1.6135, n1=16, n=26, N=10000, seed=5868)
##' # empirical power
##' # 0.813
##'
##' Sim(shape=1, S0=0.62, S1=0.62, x0=2, tf=2, rate=5, t1=3.0593,
##'     c1=-0.302, c=1.6135, n1=16, n=26, N=10000, seed=5868)
##' # empirical type-I error
##' # 0.037

Sim <- function(shape, S0, S1, x0, tf, rate, t1, c1, c, n1, n, N, restricted=0, seed=123)
{
	ta=n/rate # accrual time
	MTSL=ta+tf

	# survival function based on distribution, for calculating expected events under Null
	scale0=x0/(-log(S0))^(1/shape) # scale0 is b, shape is a, in the paper Wu 2020
	scale=x0/(-log(S1))^(1/shape) # WB under proportional hazard assumption is still WB
	s=function(a,b,u){1-pweibull(u,a,b)} # u is quantile
	H=function(a,b,u){-log(s(a,b,u))}

	set.seed(seed)
	k=N1=0
	for (i in 1:N)
	{
		#entry time, unordered
		u=runif(n, 0, ta)

		# Generate observed failure time w
		w=rweibull(n, shape, scale)

		# interim analysis
		w1=w[1:n1]
		u1=u[1:n1]

		if(restricted==0){
			# if t1-u1 >0, observed time xt1 is the smaller btw entry to interim t1-u1 and survival time w1;
			# if t1-u1 <= 0, patient enrolled after or at t1, xt1=0
			xt1=pmax(0, pmin(w1, t1-u1))
			delta1 = as.numeric(w1 < t1-u1)  # death indicator
		}
		else {
			# if t1-u1 >0, observed time xt1 is the smallest among t1-u1, tf, and w1;
			# if t1-u1 <= 0, patient enrolled after or at t1, xt1=0
			xt1=pmax(0, pmin(w1, pmin(tf, t1-u1)))
			delta1 = as.numeric(w1 < pmin(tf, t1-u1))  # status indicator
		}

		O1=sum(delta1)                          # observed death events
		M1=H(shape,scale0,xt1)                  # H(a,b,u) under Null assumption
		E1=sum(M1)                              # expected death events
		Z1=(E1-O1)/E1^(1/2)

		# final analysis
		if(restricted==0){
			xt2=pmax(0, pmin(w, MTSL-u))
			delta2 = as.numeric(w < MTSL-u)
		} else{
			xt2=pmax(0, pmin(w, tf)) # MTSL-u= ta+tf-u > tf
			delta2 = as.numeric(w < tf)
		}

		O2=sum(delta2)
		M2=H(shape,scale0,xt2)
		E2=sum(M2)
		Z=(E2-O2)/E2^(1/2)

		if (E1!=0 & E2!=0)
		{
			if (Z1>c1 & Z>c) (k=k+1) #success at final
			N1=N1+1


		}
	}
	ans=c(round(k/N1,3))
	return(ans)
}
