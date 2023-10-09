##' phase2.TTE() provides the clinical trial design solutions for two-stage
##' trials with time-to-event outcomes based on the one-sample log-rank (OSLR)
##' test. It calculates the design parameters (e.g., t1, n1, n, c1, c) using
##' optimal, minmax and admissible methods.
##'
##' @title Two-Stage Designs with TTE Outcomes Using the One-Sample Log-Rank
##' Test
##'
##' @param shape shape parameter for the baseline hazard function assuming
##'  that the failure time follows a Weibull distribution.
##' @param S0 survival probability at the fixed time point x0 under the null
##'  hypothesis (i.e, H0).
##' @param x0 a fixed time point where the survival probability is S0 under
##' the null.
##' @param hr hazard ratio, hr < 1. s1=s0^hr, where s1 is the survival
##'  probability under the alternative hypothesis (i.e., HA) and s0 is that
##'  under H0.
##' @param tf the follow-up time (restricted or unrestricted), the time period
##'  from the entry of the last patient to the end of the trial.
##' @param rate a constant accrual rate. Please consider use a reasonable rate
##'  value. If the rate is too small, the function might throw an error.
##' @param alpha type I error.
##' @param beta type II error.
##' @param prStop the lower limit of the early stopping probability under H0,
##'  with default=0.
##' @param q_value the relative importance between the maximum sample size (n)
##'  and the expected sample size under H0 (ES) when deriving the design
##'  based on the admissible method. The default is 0.5 and the range is
##'  (0, 1). The smaller q_value is, the more importance is given to ES. The
##'  greater it is, the more importance is given to n.
##' @param dfc1 the value defines the stopping criterion of the optimization
##'  process in the minmax method with smaller values lead to more iterations,
##'  default=0.001. Change is not recommended.
##' @param dfc2 the value defines the stopping criterion of the optimization
##'  process in the optimal method, smaller values lead to more iterations,
##'  default=0.001. Change is not recommended.
##' @param dfc3 the value defines the stopping criterion of the optimization
##'  process in the admissible method, smaller values lead to more iterations,
##'  default=0.001. Change is not recommended.
##' @param maxEn the maximum of the expected sample size under null, default=
##'  10000. Change is not recommended.
##' @param range the value defines how far the parameters can deviate from the
##'  last iteration in the computation of the second-stage design parameters,
##'  default=1. Change is not recommended.
##' @param t1_p1 this value defines the lower limit of the possible range of t1
##'  depending on the single-stage accrual time, default=0.2. Change is not
##'  recommended.
##' @param t1_p2 this value defines the upper limit of the possible range of t1
##'  depending on the single-stage accrual time, default=0.2. Change is not
##'  recommended.
##' @param c1_p this value defines the initial center in the possible range of
##' c1, default=0.25. Change is not recommended.
##' @param nbpt_p this value defines the initial number of points checked
##' within possible ranges for n, t1 and c1, default=11. Change is not
##' recommended.
##' @param pascote_p this value defines how fast the possible ranges of the
##'  two-stage design parameters shrink on each iteration, default=1.26. Change
##'  is not recommended.
##' @param restricted whether using restricted (1) or unrestricted (0)
##' follow-up, default = 0.
##'
##' @return The function returns a list that includes \bold{Single_stage},
##'  \bold{Two_stage_Optimal}, \bold{Two_stage_minmax}, and
##'  \bold{Two_stage_Admissible}, etc.
##' @return \bold{Single_stage} contains the design parameters for the
##' single-stage design:
##' \itemize{
##'           \item \emph{nsignle} the required sample size for the
##'           single-stage design.\cr
##'           \item \emph{tasingle} the estimated accrual time for the single-stage
##'           design.\cr
##'          \item \emph{csingle} the critical value for the
##'          single-stage design.\cr
##'          }
##' @return \bold{Two_stage_Optimal} contains the design parameters for the
##' two-stage design based on the optimal method (i.e., minimizing ES):\cr
##' \itemize{
##'          \item \emph{n1} and \emph{n} required sample sizes in the two-stage
##'           design by the interim and final stage, respectively.\cr
##'          \item \emph{c1} and \emph{c} critical values in the two-stage
##'          design for interim and final analysis, respectively.\cr
##'          \item \emph{t1} the interim analysis time in the two-stage design.\cr
##'          \item \emph{MTSL} the maximum total study length (the sum of the accrual
##'          time and the follow-up time).\cr
##'          \item \emph{ES} the expected sample size under null in the two-stage
##'           design.\cr
##'          \item \emph{PS} the probability of early stopping under null in the
##'           two-stage design.
##'           }
##' @return \bold{Two_stage_minmax} contains the design parameters for the two-
##' stage design based on the minmax method (i.e., minimizing the total sample
##'  size, n), including the same parameters as for the optimal method.
##' @return \bold{Two_stage_Admissible} contains the design parameters for the
##'  two-stage design based on the admissible method (i.e., a "compromise"
##'  between the optimal and the minmax method), including the same parameters
##'  as for the optimal method, as well as:\cr
##'  \itemize{
##'          \item \emph{Rho} The expected loss. Between the total sample size n
##'          derived from the minmax and the optimal method, the admissible
##'          method calculates a design for each possible value of n. The design
##'          with the lowest Rho value (i.e., first row in the output) is the
##'          recommended design based on the admissible method with the
##'          specified q-value.\cr
##'          }
##'  \bold{Other outputs}:\cr
##'  \itemize{
##'          \item \emph{param} The input values to the arguments.\cr
##'          \item \emph{difn_opSg} The difference in n between the single-stage
##'           design and the optimal two-stage designs.\cr
##'          \item \emph{difn_opminmax} The difference in n between the optimal
##'          and the minmax two-stage designs.\cr
##'          \item \emph{minmax.err} 0 or 1. If minmax.err=1, optimization
##'           for the minmax method is incomplete.\cr
##'          \item \emph{optimal.err} 0 or 1. If optimal.err=1, optimization for
##'           the optimal method is incomplete.\cr
##'          \item \emph{admiss.err} 0 or 1. If admiss.err=1, optimization for
##'          a given n value in the admissible method is incomplete.\cr
##'          \item \emph{admiss.null1} 0 or 1. If admiss.null1=1, the
##'          admissible result is unavailable due to incomplete optimization
##'          with either the minmax or the optimal method.\cr
##'          \item \emph{admiss.null2} 0 or 1. If admiss.null2=1, the
##'          admissible result is unavailable as either the minmax or the
##'          optimal result is unavailable.\cr
##'          \item \emph{admiss.null3} 0 or 1. If admiss.null3=1, the admissible
##'           result is unavailable as n in the minmax result and n in the
##'           optimal result are equal.
##'           }
##'
##' @import survival
##' @export
##' @examples
##' # 1. An example when q_value=0.1, i.e, more importance is given to ES.
##' # phase2.TTE(shape=0.5, S0=0.6, x0=3, hr=0.5, tf=1, rate=5,
##' # 					 alpha=0.05, beta=0.15, q_value=0.1, prStop=0, restricted=0)
##' # $param
##' #   shape  S0  hr alpha beta rate x0 tf q_value prStop restricted
##' # 1   0.5 0.6 0.5  0.05 0.15    5  3  1     0.1      0          0
##' #
##' # $Single_stage
##' #   nsingle tasingle  csingle
##' # 1      45        9 1.644854
##' #
##' # $Two_stage_Optimal
##' #   n1     c1  n      c     t1 MTSL    ES     PS
##' # 1 29 0.1389 48 1.6159 5.7421 10.6 37.29 0.5552
##' #
##' # $Two_stage_minmax
##' #   n1     c1  n      c     t1 MTSL      ES     PS
##' # 1 34 0.1151 45 1.6391 6.7952   10 38.9831 0.5458
##' #
##' # $Two_stage_Admissible
##' #      n1      c1  n      c     t1 MTSL      ES     PS      Rho
##' # 123  29  0.0705 47 1.6232 5.7261 10.4 37.2993 0.5281 38.26937
##' # 285  28  0.0792 48 1.6171 5.5663 10.6 37.2790 0.5316 38.35110
##' # 1701 31  0.0733 46 1.6293 6.0191 10.2 37.5828 0.5292 38.42452
##' # 170  33 -0.0405 45 1.6391 6.4245 10.0 38.7692 0.4839 39.39228
##' #
##' # $difn_opSg
##' # [1] 3
##' #
##' # $difn_opminmax
##' # [1] 3
##' #
##' # $minmax.err
##' # [1] 0
##' #
##' # $optimal.err
##' # [1] 0
##' #
##' # $admiss.err
##' # [1] 0
##' #
##' # $admiss.null1
##' # [1] 0
##' #
##' # $admiss.null2
##' # [1] 0
##' #
##' # $admiss.null3
##' # [1] 0
##'
##' # 2. An example when q_value=0.75, i.e., more importance is given to n.
##' # phase2.TTE(shape=0.5, S0=0.6, x0=3, hr=0.5, tf=1, rate=5,
##' # alpha=0.05, beta=0.15, q_value=0.75, prStop=0, restricted=0)
##' # $param
##' #   shape  S0  hr alpha beta rate x0 tf q_value prStop restricted
##' # 1   0.5 0.6 0.5  0.05 0.15    5  3  1    0.75      0          0
##' #
##' # $Single_stage
##' #   nsingle tasingle  csingle
##' # 1      45        9 1.644854
##' #
##' # $Two_stage_Optimal
##' #   n1     c1  n      c     t1 MTSL    ES     PS
##' # 1 29 0.1389 48 1.6159 5.7421 10.6 37.29 0.5552
##' #
##' # $Two_stage_minmax
##' #   n1     c1  n      c     t1 MTSL      ES     PS
##' # 1 34 0.1151 45 1.6391 6.7952   10 38.9831 0.5458
##' #
##' # $Two_stage_Admissible
##' #      n1      c1  n      c     t1 MTSL      ES     PS      Rho
##' # 170  33 -0.0405 45 1.6391 6.4245 10.0 38.7692 0.4839 43.44230
##' # 1701 31  0.0733 46 1.6293 6.0191 10.2 37.5828 0.5292 43.89570
##' # 123  29  0.0705 47 1.6232 5.7261 10.4 37.2993 0.5281 44.57483
##' # 285  28  0.0792 48 1.6171 5.5663 10.6 37.2790 0.5316 45.31975
##' #
##' # $difn_opSg
##' # [1] 3
##' #
##' # $difn_opminmax
##' # [1] 3
##' #
##' # $minmax.err
##' # [1] 0
##' #
##' # $optimal.err
##' # [1] 0
##' #
##' # $admiss.err
##' # [1] 0
##' #
##' # $admiss.null1
##' # [1] 0
##' #
##' # $admiss.null2
##' # [1] 0
##' #
##' # $admiss.null3
##' # [1] 0
##' @references
##' Wu, J, Chen L, Wei J, Weiss H, Chauhan A. (2020). Two-stage phase II survival trial design. Pharmaceutical Statistics. 2020;19:214-229. https://doi.org/10.1002/pst.1983


phase2.TTE <- function(shape, S0, x0, hr, tf, rate, alpha, beta, prStop=0,
											q_value=0.5, dfc1=0.001, dfc2=0.001, dfc3=0.001,
											maxEn=10000, range=1, t1_p1=0.2, t1_p2=1.2,
											c1_p=0.25, nbpt_p=11, pascote_p=1.26, restricted=0) {

 if (restricted==0) {
	# error flags
	minmax.err <- optimal.err <- admiss.err <- admiss.null1 <- admiss.null2 <- admiss.null3 <- 0

	difn_opSg <- difn_opminmax <- NA

	calculate_alpha <- function(c2, c1, rho0) {
		fun1 <- function(z, c1, rho0) {
			f <- dnorm(z) * pnorm((rho0 * z - c1)/sqrt(1 - rho0^2))
			return(f)
		}
		alpha <- integrate(fun1, lower = c2, upper = Inf, c1, rho0)$value
		return(alpha)
	}

	calculate_power <- function(cb, cb1, rho1) {
		fun2 <- function(z, cb1, rho1) {
			f <- dnorm(z) * pnorm((rho1 * z - cb1)/sqrt(1 - rho1^2))
			return(f)
		}
		pwr <- integrate(fun2, lower = cb, upper = Inf, cb1 = cb1, rho1 = rho1)$value
		return(pwr)
	}

	###########################################################################
	fct <- function(zi, ceps = 1e-04, alphaeps = 1e-04, nbmaxiter = 100) {
		ta <- as.numeric(zi[1])
		t1 <- as.numeric(zi[2])
		c1 <- as.numeric(zi[3])

		f0 = function(u) {
			(shape/scale0) * (u/scale0)^(shape - 1) * exp(-(u/scale0)^shape)
		}
		s0 = function(u) {
			exp(-(u/scale0)^shape)
		}
		h0 = function(u) {
			(shape/scale0) * (u/scale0)^(shape - 1)
		}
		H0 = function(u) {
			(u/scale0)^shape
		}
		s = function(b, u) {
			exp(-b * (u/scale0)^shape)
		}
		h = function(b, u) {
			b * h0(u)
		}
		H = function(b, u) {
			b * H0(u)
		}
		scale0 = x0/(-log(S0))^(1/shape)
		scale1 = hr
		# above code differs from rKJ

		G = function(t) {
			1 - punif(t, tf, ta + tf)
		} # 214-231 differs from restricted follow-up

		g0 = function(t) {
			s(scale1, t) * h0(t) * G(t)
		}
		g1 = function(t) {
			s(scale1, t) * h(scale1, t) * G(t)
		}
		g00 = function(t) {
			s(scale1, t) * H0(t) * h0(t) * G(t)
		}
		g01 = function(t) {
			s(scale1, t) * H0(t) * h(scale1, t) * G(t)
		}
		p0 = integrate(g0, 0, ta + tf)$value
		p1 = integrate(g1, 0, ta + tf)$value
		p00 = integrate(g00, 0, ta + tf)$value
		p01 = integrate(g01, 0, ta + tf)$value
		sigma2.1 = p1 - p1^2 + 2 * p00 - p0^2 - 2 * p01 + 2 * p0 * p1
		sigma2.0 = p0
		om = p0 - p1

		G1 = function(t) {
			1 - punif(t, t1 - ta, t1)
		}
		g0 = function(t) {
			s(scale1, t) * h0(t) * G1(t)
		}
		g1 = function(t) {
			s(scale1, t) * h(scale1, t) * G1(t)
		}
		g00 = function(t) {
			s(scale1, t) * H0(t) * h0(t) * G1(t)
		}
		g01 = function(t) {
			s(scale1, t) * H0(t) * h(scale1, t) * G1(t)
		}
		p0 = integrate(g0, 0, ta + tf)$value # differs from restricted fu
		p1 = integrate(g1, 0, ta + tf)$value
		p00 = integrate(g00, 0, ta + tf)$value
		p01 = integrate(g01, 0, ta + tf)$value
		sigma2.11 = p1 - p1^2 + 2 * p00 - p0^2 - 2 * p01 + 2 * p0 * p1
		sigma2.01 = p0
		om1 = p0 - p1

		q1 = function(t) {
			s0(t) * h0(t) * G1(t)
		}
		q = function(t) {
			s0(t) * h0(t) * G(t)
		}
		v1 = integrate(q1, 0, ta + tf)$value # differs from restricted fu
		v = integrate(q, 0, ta + tf)$value
		rho0 = sqrt(v1/v)
		rho1 <- sqrt(sigma2.11/sigma2.1)


		########## Calculate Type I error alpha ###############################

		cL <- (-10)
		cU <- (10)
		alphac <- 100
		iter <- 0
		while ((abs(alphac - alpha) > alphaeps | cU - cL > ceps) & iter < nbmaxiter) {
			iter <- iter + 1
			c <- (cL + cU)/2
			alphac <- calculate_alpha(c, c1, rho0)
			if (alphac > alpha) {
				cL <- c
			} else {
				cU <- c
			}
		}

		################ Calculate power ###########################################

		cb1 <- sqrt((sigma2.01/sigma2.11)) * (c1 - (om1 * sqrt(rate * t1)/sqrt(sigma2.01)))
		cb <- sqrt((sigma2.0/sigma2.1)) * (c - (om * sqrt(rate * ta))/sqrt(sigma2.0))
		pwrc <- calculate_power(cb = cb, cb1 = cb1, rho1 = rho1)
		res <- c(cL, cU, alphac, 1 - pwrc, rho0, rho1, cb1, cb)
		return(res)
	}

	##############################################################################
	##############################################################################
	c1 <- 0
	rho0 <- 0
	cb1 <- 0
	rho1 <- 0
	hz <- c(0, 0)
	ceps <- 0.001
	alphaeps <- 0.001
	nbmaxiter <- 100

	Duration = function(shape, S0, x0, hr, tf, rate, alpha, beta) {

		f0 = function(u) {
			(shape/scale0) * (u/scale0)^(shape - 1) * exp(-(u/scale0)^shape)
		}
		s0 = function(u) {
			exp(-(u/scale0)^shape)
		}
		h0 = function(u) {
			(shape/scale0) * (u/scale0)^(shape - 1)
		}
		H0 = function(u) {
			(u/scale0)^shape
		}
		s = function(b, u) {
			exp(-b * (u/scale0)^shape)
		}
		h = function(b, u) {
			b * h0(u)
		}
		H = function(b, u) {
			b * H0(u)
		}
		scale0 = x0/(-log(S0))^(1/shape)
		scale1 = hr
		#}

		root = function(ta) {
			tau = ta + tf
			G = function(t) {
				1 - punif(t, tf, tau)
			}
			g0 = function(t) {
				s(scale1, t) * h0(t) * G(t)
				# if tau is too big, h0(t) can be NA, maybe try use a smaller parameter before tf
			}
			g1 = function(t) {
				s(scale1, t) * h(scale1, t) * G(t)
			}
			g00 = function(t) {
				s(scale1, t) * H0(t) * h0(t) * G(t)
			}
			g01 = function(t) {
				s(scale1, t) * H0(t) * h(scale1, t) * G(t)
			}

			p0 = integrate(g0, 0, tau)$value
			p1 = integrate(g1, 0, tau)$value
			p00 = integrate(g00, 0, tau)$value
			p01 = integrate(g01, 0, tau)$value
			s1 = sqrt(p1 - p1^2 + 2 * p00 - p0^2 - 2 * p01 + 2 * p0 * p1)
			s0 = sqrt(p0)
			om = p0 - p1
			rate * ta - (s0 * qnorm(1 - alpha) + s1 * qnorm(1 - beta))^2/om^2
		}

		tasingle = tryCatch(
			expr={uniroot(root, lower = 0, upper = 10 * tf)$root},

			error=function(e){
				stop("Solution for ta in single stage cannot be found.
             You may consider using a longer tf or a faster rate.")}
		)

		nsingle <- ceiling(tasingle * rate)
		tasingle = ceiling(tasingle)
		ans = list(nsingle = nsingle, tasingle = tasingle)
		return(ans)
	}

	################################# Single Stage ###############################
	ans = Duration(shape, S0, x0, hr, tf, rate, alpha, beta)
	tasingle <- ans$tasingle
	nsingle = ans$nsingle

	Single_stage <- data.frame(nsingle = nsingle, tasingle = tasingle, csingle = qnorm(1 - alpha))

	################################# Two-Stage Optimal ##########################

	atc0 <- data.frame(n = nsingle, t1 = tasingle, c1 = c1_p)
	nbpt <- nbpt_p
	pascote <- pascote_p
	cote <- range * pascote
	c1.lim <- nbpt * c(-1, 1)
	EnH0 <- maxEn
	iter <- 0

	while (iter < nbmaxiter & diff(c1.lim)/nbpt > dfc2) {
		iter <- iter + 1
		cat("Calculating optimal results, iter=", iter, "&EnH0=", round(EnH0, 2),
				"& Dc1/nbpt=", round(diff(c1.lim)/nbpt, 5), "\n", sep = "")

		if (iter%%2 == 0)
			nbpt <- nbpt + 1
		cote <- cote/pascote
		n.lim <- atc0$n + c(-1, 1) * nsingle * cote
		t1.lim <- atc0$t1 + c(-1, 1) * tasingle * cote
		c1.lim <- atc0$c1 + c(-1, 1) * cote
		ta.lim <- n.lim/rate
		t1.lim <- pmax(0, t1.lim)
		n <- seq(n.lim[1], n.lim[2], l = nbpt)
		n <- ceiling(n)
		n <- unique(n)
		ta <- n/rate
		t1 <- seq(t1.lim[1], t1.lim[2], l = nbpt)
		c1 <- seq(c1.lim[1], c1.lim[2], l = nbpt)

		## set c1 > qnorm(prStop)
		c1 <- c1[c1 >= qnorm(prStop)]

		if(length(c1)==0) {
			warning("length(c1)==0 in iter ", iter, ". Check if prStop is too large. ",
							"Optimization halted, the output (Two_stage_Optimal) might not be optimal.")
			optimal.err <- 1
			break
		}

		ta <- ta[ta > 0]
		t1 <- t1[t1 >= t1_p1 * tasingle & t1 <= t1_p2 * tasingle]
		z <- expand.grid(list(ta = ta, t1 = t1, c1 = c1))
		z <- z[z$ta > z$t1, ]
		nz <- dim(z)[1]

		z$pap <- pnorm(z$c1)
		z$eta <- z$ta - pmax(0, z$ta - z$t1) * z$pap
		z$enh0 <- z$eta * rate   # estimated sample size under null
		z <- z[z$enh0 <= EnH0, ]
		nz1 <- dim(z)[1]
		if (nz1==0) next

		resz <- t(apply(z, 1, fct, ceps = ceps, alphaeps = alphaeps, nbmaxiter = nbmaxiter))
		resz <- as.data.frame(resz)
		names(resz) <- c("cL", "cU", "alphac", "betac", "rho0", "rho1", "cb1", "cb")

		r <- cbind(z, resz)
		r$n <- r$ta * rate
		r$c <- (r$cL + r$cU)/2
		r$diffc <- r$cU - r$cL
		r$diffc <- ifelse(r$diffc <= ceps, 1, 0)
		r <- r[1 - r$betac >= 1 - beta, ]             # select ones reaching required power
		r <- r[order(r$enh0), ]


		if (dim(r)[1] > 0) {
			atc <- r[1, ]

			if (atc$enh0 < EnH0) {
				EnH0 <- atc$enh0
				atc0 <- atc
			}
		} else {
			atc <- data.frame(ta = NA, t1 = NA, c1 = NA, n = NA,
												pap = NA, eta = NA, enh0 = NA,cL = NA,
												cU = NA,  alphac= NA, betac = NA, rho0 = NA,
												rho1 = NA, cb1 = NA, cb  = NA, c = NA, diffc=NA)
		}

		if (iter == 1) {
			atcs <- atc
		} else {
			atcs <- rbind(atcs, atc)
		}
	}

	a <- atcs[!is.na(atcs$n), ]

	a <- a[!(a$n==ceiling(a$t1*rate)), ] ## remove designs with n=n1

	res <- a[order(a$enh0), ]
	des <- round(res[1, ], 4)

	Two_stage_Optimal <- data.frame(n1 = ceiling(des$t1 * rate), c1 = des$c1,
																	n = ceiling(des$n), c = des$c,
																	t1 = des$t1, MTSL = des$ta + tf,
																	ES = des$enh0, PS = des$pap)

	###################################### Two stage Minmax ##############################################
	if(anyNA(Two_stage_Optimal$n)==1){
		maxn <- maxEn
	} else {maxn <- Two_stage_Optimal$n}

	atc0 <- data.frame(n = nsingle, t1 = tasingle, c1 = c1_p)
	nbpt <- nbpt_p
	pascote <- pascote_p
	cote <- range* pascote
	c1.lim <- nbpt * c(-1, 1)
	iter <- 0

	while (iter < nbmaxiter & diff(c1.lim)/nbpt > dfc1) {
		iter <- iter + 1
		cat("Calculating minmax results, iter=", iter, " & n=", maxn, "& Dc1/nbpt=", round(diff(c1.lim)/nbpt, 5), "\n", sep = "")

		if (iter%%2 == 0) nbpt <- nbpt + 1             # increase number of points to be checked within range

		cote <- cote/pascote                           # decrease the parameter cote to decrease search range
		n.lim <- atc0$n + c(-1, 1) * nsingle * cote    # atc0 is optimal design from last iteration
		t1.lim <- atc0$t1 + c(-1, 1) * tasingle * cote

		c1.lim <- atc0$c1 + c(-1, 1) * cote            # range of c1

		ta.lim <- n.lim/rate
		t1.lim <- pmax(0, t1.lim)

		n <- seq(n.lim[1], n.lim[2], l = nbpt)
		n <- ceiling(n)
		n <- unique(n)
		ta <- n/rate
		t1 <- seq(t1.lim[1], t1.lim[2], l = nbpt)
		c1 <- seq(c1.lim[1], c1.lim[2], l = nbpt)

		## set c1 > qnorm(prStop)

		c1 <- c1[c1 >= qnorm(prStop)]

		# if length(c1)==0, no need to go to next iter b/c a narrower c1.lim
		# centered at same c1 will also fail, but it might be useful to keep results
		# from previous iterations
		if(length(c1)==0) {
			warning("length(c1)==0 in iter ", iter, ". Check if prStop is too large. "
							, "Optimization halted, the output (Two_stage_minmax) might not be
							optimal.")
			minmax.err <- 1
			break
		}



		ta <- ta[ta > 0]
		t1 <- t1[t1 >= t1_p1 * tasingle & t1 <= t1_p2 * tasingle] # set to parameter

		z <- expand.grid(list(ta = ta, t1 = t1, c1 = c1))
		z <- z[z$ta > z$t1, ]
		nz <- dim(z)[1]
		z$pap <- pnorm(z$c1)                          ## stopping prob under null
		z$eta <- z$ta - pmax(0, z$ta - z$t1) * z$pap
		z$enh0 <- z$eta * rate                        ## expected sample size under null

		z$n <- z$ta*rate                              ## keep design with n < maxn
		z <- z[z$n <= maxn, ]

		nz1 <- dim(z)[1]
		if (nz1 ==0) break # the vector z$n will always contain maxn, due to ceiling()

		### For each row of z, apply fct
		resz <- t(apply(z, 1, fct, ceps = ceps, alphaeps = alphaeps, nbmaxiter = nbmaxiter))

		resz <- as.data.frame(resz)
		names(resz) <- c("cL", "cU", "alphac", "betac", "rho0", "rho1", "cb1", "cb")
		r <- cbind(z, resz)
		r$c <- r$cL + (r$cU - resz$cL)/2             # c, criterion for final analysis
		r$diffc <- r$cU - r$cL
		r$diffc <- ifelse(r$diffc <= ceps, 1, 0)
		r <- r[1 - r$betac >= 1 - beta, ]            ## determine if reaches required power

		r <- r[order(r$n, r$enh0), ]                 # order by max sample size and enh0

		if (dim(r)[1] > 0) {
			atc <- r[1, ]                              # pick the row with smallest n
			if (atc$n < maxn) {
				maxn <- atc$n                            # update upper limit of n, maxn to the smallest n in r
				atc0 <- atc                              # update atc0, the (ta, t1, c1, n) w/ smallest n in this iteration
			}
		} else {
			atc <- data.frame(ta = NA, t1 = NA, c1 = NA, n = NA,
												pap = NA, eta = NA, enh0 = NA,cL = NA,
												cU = NA,  alphac= NA, betac = NA, rho0 = NA,
												rho1 = NA, cb1 = NA, cb  = NA, c = NA, diffc= NA)
		}

		atc$iter <- iter
		atc$tai <- ta.lim[1]
		atc$tas <- ta.lim[2]
		atc$ti <- t1.lim[1]
		atc$ts <- t1.lim[2]
		atc$ci <- c1.lim[1]
		atc$cs <- c1.lim[2]
		atc$cote <- cote
		if (iter == 1) {
			atcs <- atc
		} else {
			atcs <- rbind(atcs, atc)
		}
	}

	a <- atcs[!is.na(atcs$n), ]

	a <- a[!(a$n==ceiling(a$t1*rate)), ] ## remove designs with n=n1

	res <- a[order(a$n, a$enh0), ]     ## find the design minimizes n and secondarily enh0

	des <- round(res[1, ], 4)

	Two_stage_minmax <- data.frame(n1 = ceiling(des$t1 * rate), c1 = des$c1,
																 n = des$n, c = des$c,
																 t1 = des$t1, MTSL = des$ta + tf,
																 ES = des$enh0, PS = des$pap)

	####################################### Two stage Admissible ##############################################

	if (minmax.err ==1|optimal.err==1){
		message("Admissible design unavailabe: Incomplete optimization in minmax
						or optimal results")
		admiss.null1 <- 1
		Two_stage_Admissible=NULL

	}else if (anyNA(c(Two_stage_minmax$n, Two_stage_Optimal$n))){
		message("Admissible design unavailabe: n = NA in minmax or optimal results")
		admiss.null2 <- 1
		Two_stage_Admissible=NULL
	}
	else if (Two_stage_Optimal$n-Two_stage_minmax$n == 0){
		message("Admissible design unavailabe: Two_stage_Optimal$n - Two_stage_minmax$n = 0.")
		admiss.null3 <- 1
		Two_stage_Admissible=NULL

	}else{
		nvec <- seq(Two_stage_minmax$n, Two_stage_Optimal$n, by = 1)

		for (i in 1:length(nvec)) {

			atc0 <- data.frame(n=nvec[i], t1 = tasingle, c1 = c1_p)
			nbpt <- nbpt_p
			pascote <- pascote_p
			cote <- range* pascote
			c1.lim <- nbpt * c(-1, 1)
			EnH0 <- maxEn
			iter <- 0

			while (iter < nbmaxiter & diff(c1.lim)/nbpt > dfc3) {
				iter <- iter + 1
				cat("Calculating admissible results, n=", nvec[i], " & EnH0=", round(EnH0, 2), "& Dc1/nbpt=", round(diff(c1.lim)/nbpt, 5), "\n", sep = "")

				if (iter%%2 == 0) nbpt <- nbpt + 1             # increase number of points to be checked within range

				cote <- cote/pascote                           # decrease the parameter cote to decrease search range

				n <- atc0$n

				t1.lim <- atc0$t1 + c(-1, 1) * tasingle * cote

				c1.lim <- atc0$c1 + c(-1, 1) * cote            # range of c1

				t1.lim <- pmax(0, t1.lim)

				ta <- n/rate
				t1 <- seq(t1.lim[1], t1.lim[2], l = nbpt)
				c1 <- seq(c1.lim[1], c1.lim[2], l = nbpt)


				## set c1 > qnorm(prStop)
				c1 <- c1[c1 >= qnorm(prStop)]

				# if length(c1)==0, no need to go to next iter b/c a narrower c1.lim
				# will also fail, but it might be useful to keep results from previous
				# iterations
				if(length(c1)==0) {
					warning("length(c1)==0 in iter ", iter,
									". The admissible result might not be optimized for n=",
									nvec[i])
					admiss.err <- 1
					break

				}
				t1 <- t1[t1 >= t1_p1 * tasingle & t1 <= t1_p2 * tasingle]

				z <- expand.grid(list(ta = ta, t1 = t1, c1 = c1))
				z <- z[z$ta > z$t1, ]
				z$pap <- pnorm(z$c1)                          ## stopping prob under null
				z$eta <- z$ta - pmax(0, z$ta - z$t1) * z$pap
				z$enh0 <- z$eta * rate                        ## expected sample size under null
				z <- z[z$enh0 <= EnH0, ]                      ## exclude combos with enh0 higher than EnH0

				### For each row of z, apply fct
				resz <- t(apply(z, 1, fct, ceps = ceps, alphaeps = alphaeps, nbmaxiter = nbmaxiter))

				resz <- as.data.frame(resz)
				names(resz) <- c("cL", "cU", "alphac", "betac", "rho0", "rho1", "cb1", "cb")
				r <- cbind(z, resz)
				r$n <- n
				r$c <- (r$cU + resz$cL)/2             # c, criterion for final analysis
				r <- r[1 - r$betac >= 1 - beta, ]            ## determine if reaches required power
				r <- r[order(r$enh0), ]                      # order by expected sample size under null

				if (dim(r)[1] > 0) {
					atc <- r[1, ]                          # pick the row with smallest enh0

					if (atc$enh0 < EnH0) {
						EnH0 <- atc$enh0                         # update upper limit of enh0, EnH0 to the smallest enh0 in r
						atc0 <- atc                             # update atc0, the (ta, t1, c1, n) w/ smallest enh0 in this iteration
					}
				} else {
					atc <- data.frame(ta = NA, t1 = NA, c1 = NA, n = NA,
														pap = NA, eta = NA, enh0 = NA,cL = NA,
														cU = NA,  alphac= NA, betac = NA, rho0 = NA,
														rho1 = NA, cb1 = NA, cb  = NA, c = NA)
				}

				atc$iter <- iter
				atc$EnH0 <- EnH0
				atc$ti <- t1.lim[1]
				atc$ts <- t1.lim[2]
				atc$ci <- c1.lim[1]
				atc$cs <- c1.lim[2]
				atc$cote <- cote

				if (iter == 1) {
					atcs <- atc
				} else {
					atcs <- rbind(atcs, atc)
				}
			}

			a <- atcs[!is.na(atcs$n), ]

			a <- a[!(a$n==ceiling(a$t1*rate)), ] ## remove designs with n=n1

			res <- a[order(a$enh0), ] ## find the optimal design minimizes enh0

			if (i==1) {des3 <- round(res[1, ], 4)
			}else des3 <- rbind(des3, round(res[1, ], 4))

		}

		# find design minimize rho
		des3$Rho <- q_value*des3$n + (1 - q_value)*des3$enh0
		des3$n1 <- ceiling(des3$t1 * rate)
		des3$MTSL <- des3$ta + tf
		des3$ES <-  des3$enh0
		des3$PS <- des3$pap
		des3 <- des3[, c("n1", "c1", "n", "c", "t1", "MTSL", "ES", "PS", "Rho")]
		des3 <- des3[order(des3$Rho), ]

		Two_stage_Admissible <- des3

		# difn_opminmax <- Two_stage_Optimal$n-Two_stage_minmax$n
		# difn_opSg<- Two_stage_Optimal$n-Single_stage$n
	}
###############################################################################
#                           restricted follow-up                              #
###############################################################################
 } else if (restricted==1) {
	# error flags
	minmax.err <- optimal.err <- admiss.err <- admiss.null1 <- admiss.null2 <- admiss.null3 <- 0
	difn_opSg <- difn_opminmax <- NA

	calculate_alpha <- function(c2, c1, rho0) {
		fun1 <- function(z, c1, rho0) {
			f <- dnorm(z) * pnorm((rho0 * z - c1)/sqrt(1 - rho0^2))
			return(f)
		}
		alpha <- integrate(fun1, lower = c2, upper = Inf, c1, rho0)$value
		return(alpha)
	}

	calculate_power <- function(cb, cb1, rho1) {
		fun2 <- function(z, cb1, rho1) {
			f <- dnorm(z) * pnorm((rho1 * z - cb1)/sqrt(1 - rho1^2))
			return(f)
		}
		pwr <- integrate(fun2, lower = cb, upper = Inf, cb1 = cb1, rho1 = rho1)$value
		return(pwr)
	}

	##############################################################################
	fct <- function(zi, ceps = 1e-04, alphaeps = 1e-04, nbmaxiter = 100) {
		# 49-238
		ta <- as.numeric(zi[1])
		t1 <- as.numeric(zi[2])
		c1 <- as.numeric(zi[3])

		s0 = function(u) {
			1 - pweibull(u, shape, scale0)
		}
		f0 = function(u) {
			dweibull(u, shape, scale0)
		}
		h0 = function(u) {
			f0(u)/s0(u)
		}
		H0 = function(u) {
			-log(s0(u))
		}
		s = function(b, u) {
			s0(u)^b
		}
		h = function(b, u) {
			b * f0(u)/s0(u)
		}
		H = function(b, u) {
			-b * log(s0(u))
		}
		scale0 = x0/(-log(S0))^(1/shape)
		scale1 = hr
		## when b=scale1(=hr) give the survival, hazard and cumulative hazard at
		## alternative

		############################################################################
		g0 = function(t) {
			s(scale1, t) * h0(t)
		}
		g1 = function(t) {
			s(scale1, t) * h(scale1, t)
		}
		g00 = function(t) {
			s(scale1, t) * H0(t) * h0(t)
		}
		g01 = function(t) {
			s(scale1, t) * H0(t) * h(scale1, t)
		}
		p0 = integrate(g0, 0, tf)$value
		p1 = integrate(g1, 0, tf)$value
		p00 = integrate(g00, 0, tf)$value
		p01 = integrate(g01, 0, tf)$value
		sigma2.1 = p1 - p1^2 + 2 * p00 - p0^2 - 2 * p01 + 2 * p0 * p1
		# exact var(W) under H1
		sigma2.0 = p0
		om = p0 - p1

		############################################################################
		G1 = function(t) {
			1 - punif(t, t1 - ta, t1)
		}
		g0 = function(t) {
			s(scale1, t) * h0(t) * G1(t)
		}
		g1 = function(t) {
			s(scale1, t) * h(scale1, t) * G1(t)
		}
		g00 = function(t) {
			s(scale1, t) * H0(t) * h0(t) * G1(t)
		}
		g01 = function(t) {
			s(scale1, t) * H0(t) * h(scale1, t) * G1(t)
		}
		p0 = integrate(g0, 0, tf)$value
		p1 = integrate(g1, 0, tf)$value
		p00 = integrate(g00, 0, tf)$value
		p01 = integrate(g01, 0, tf)$value
		sigma2.11 = p1 - p1^2 + 2 * p00 - p0^2 - 2 * p01 + 2 * p0 * p1
		# exact var(W1) under H1
		sigma2.01 = p0
		om1 = p0 - p1

		###########################################################################
		q1 = function(t) {
			s0(t) * h0(t) * G1(t)
		}
		q = function(t) {
			s0(t) * h0(t)
		}
		v1 = integrate(q1, 0, tf)$value
		v = integrate(q, 0, tf)$value
		rho0 = sqrt(v1/v)                 # correlation between W1 and W under H0?
		rho1 <- sqrt(sigma2.11/sigma2.1) # correlation between W1 and W under H1?

		################ Calculate Type I error alpha ##############################
		cL <- (-10)
		cU <- (10)
		alphac <- 100
		iter <- 0
		while ((abs(alphac - alpha) > alphaeps | cU - cL > ceps) & iter < nbmaxiter) {
			iter <- iter + 1
			c <- (cL + cU)/2                        # c2 in calculate_alpha()
			alphac <- calculate_alpha(c, c1, rho0)
			if (alphac > alpha) {
				cL <- c # if alphac is too big, increase the lower limit of alpha
			} else {
				cU <- c
			}
		}

		################################### Calculate power ########################

		cb1 <- sqrt((sigma2.01/sigma2.11)) * (c1 - (om1 * sqrt(rate * t1)/sqrt(sigma2.01))) # rate*t1=n1
		cb <- sqrt((sigma2.0/sigma2.1)) * (c - (om * sqrt(rate * ta))/sqrt(sigma2.0)) # rate*ta=n
		pwrc <- calculate_power(cb = cb, cb1 = cb1, rho1 = rho1)

		res <- c(cL, cU, alphac, 1 - pwrc, rho0, rho1, cb1, cb)
		return(res)
	}

	##############################################################################
	##############################################################################

	c1 <- 0
	rho0 <- 0
	cb1 <- 0
	rho1 <- 0
	hz <- c(0, 0)
	ceps <- 0.001
	alphaeps <- 0.001
	nbmaxiter <- 100

	s0 = function(u) {
		1 - pweibull(u, shape, scale0)
	}
	f0 = function(u) {
		dweibull(u, shape, scale0)
	}
	h0 = function(u) {
		f0(u)/s0(u)
	}
	H0 = function(u) {
		-log(s0(u))
	}
	s = function(b, u) {
		s0(u)^b
	}
	h = function(b, u) {
		b * f0(u)/s0(u)
	}
	H = function(b, u) {
		-b * log(s0(u))
	}
	scale0 = x0/(-log(S0))^(1/shape)
	scale1 = hr

	##################################### Single stage ###########################
	g0 = function(t) {
		s(scale1, t) * h0(t)
	}
	g1 = function(t) {
		s(scale1, t) * h(scale1, t)
	}
	g00 = function(t) {
		s(scale1, t) * H0(t) * h0(t)
	}
	g01 = function(t) {
		s(scale1, t) * H0(t) * h(scale1, t)
	}
	p0 = integrate(g0, 0, tf)$value    # eq(3)
	p1 = integrate(g1, 0, tf)$value    # eq(4)
	p00 = integrate(g00, 0, tf)$value  # eq(5)
	p01 = integrate(g01, 0, tf)$value  # eq(6)
	s1 = sqrt(p1 - p1^2 + 2 * p00 - p0^2 - 2 * p01 + 2 * p0 * p1)
	s0 = sqrt(p0)
	om = p0 - p1 # omega in eq(2)

	#############################

	nsingle <- (s0 * qnorm(1 - alpha) + s1 * qnorm(1 - beta))^2/om^2 # eq(2)
	nsingle <- ceiling(nsingle)
	tasingle <- nsingle/rate

	Single_stage <- data.frame(nsingle = nsingle, tasingle = tasingle,
														 csingle = qnorm(1 - alpha))

	######################### Two stage Optimal ##################################
	atc0 <- data.frame(n = nsingle, t1 = tasingle, c1 = 0.25)
	atcs <- NULL
	nbpt <- nbpt_p
	pascote <- pascote_p
	cote <- range * pascote
	c1.lim <- nbpt * c(-1, 1)
	EnH0 <- maxEn
	iter <- 0

	while (iter < nbmaxiter & diff(c1.lim)/nbpt > dfc2 ) {
		iter <- iter + 1
		cat("Calculating optimal results, iter=", iter, "&EnH0=", round(EnH0, 2),
				"& Dc1/nbpt=", round(diff(c1.lim)/nbpt, 5), "\n", sep = "")
		if (iter%%2 == 0) nbpt <- nbpt + 1
		# increase number of points to be checked within range
		cote <- cote/pascote
		# decrease the parameter cote to decrease search range

		# set the search range of ta, t1 and c1
		# atc0 is optimal design from last iteration
		n.lim <- atc0$n + c(-1, 1) * nsingle * cote
		ta.lim <- n.lim/rate
		t1.lim <- atc0$t1 + c(-1, 1) * tasingle * cote
		t1.lim <- pmax(0, t1.lim)
		c1.lim <- atc0$c1 + c(-1, 1) * cote  # range of c1

		# search points of t1, ta, c1
		n <- seq(n.lim[1], n.lim[2], l = nbpt)
		n <- ceiling(n)
		n <- unique(n)
		ta <- n/rate
		ta <- ta[ta > 0]
		t1 <- seq(t1.lim[1], t1.lim[2], l = nbpt)
		t1 <- t1[t1 >=  t1_p1 * tasingle & t1 <= t1_p2 * tasingle]
		c1 <- seq(c1.lim[1], c1.lim[2], l = nbpt)
		## set c1 > qnorm(prStop)
		c1 <- c1[c1 >= qnorm(prStop)]
		# if length(c1)==0, no need to go to next iter b/c a narrower c1.lim
		# will also fail, but it might be useful to keep results from previous
		# iterations
		if(length(c1)==0) {
			warning("length(c1)==0 in iter ", iter, ". Check if prStop is too large.
        				", "Optimization halted, the output (Two_stage_Optimal) might
        				not be optimal.")
			optimal.err <- 1
			break
		}

		# create a grid
		z <- expand.grid(list(ta = ta, t1 = t1, c1 = c1))
		z <- z[z$ta > z$t1, ]
		z$pap <- pnorm(z$c1)                          ## stopping prob under null
		z$eta <- z$ta - pmax(0, z$ta - z$t1) * z$pap

		z$enh0 <- z$eta * rate                 ## expected sample size under null
		z <- z[z$enh0 <= EnH0, ]     ## exclude combos with enh0 higher than EnH0
		nz1 <- dim(z)[1]
		if(nz1==0) next

		### For each row of z, apply fct to get type I and II errors and c
		resz <- t(apply(z, 1, fct, ceps = ceps, alphaeps = alphaeps, nbmaxiter =
											nbmaxiter))
		resz <- as.data.frame(resz)
		names(resz) <- c("cL", "cU", "alphac", "betac", "rho0", "rho1", "cb1", "cb")
		r <- cbind(z, resz)
		r$n <- r$ta * rate
		r$c <- (r$cL + r$cU)/2   # c, criterion for final analysis
		# r$diffc <- r$cU - r$cL
		# r$diffc <- ifelse(r$diffc <= ceps, 1, 0)
		r <- r[1 - r$betac >= 1 - beta, ]  # if reach required power

		# for current inter, find the combo w/ smallest enh0 that < EnH0
		r <- r[order(r$enh0), ]  # order by expected sample size under nul
		if (dim(r)[1] > 0) {
			atc <- r[1, ]          # pick the row with smallest enh0
			if (atc$enh0 < EnH0) {
				EnH0 <- atc$enh0
				# update upper limit of enh0, EnH0 to the smallest enh0 in r
				atc0 <- atc
				# update atc0 (ta, t1, c1, n) w/ smallest enh0 in this iteration
			}
		} else {
			atc <- data.frame(ta = NA, t1 = NA, c1 = NA, n = NA,
												pap = NA, eta = NA, enh0 = NA,cL = NA,
												cU = NA,  alphac= NA, betac = NA, rho0 = NA,
												rho1 = NA, cb1 = NA, cb  = NA, c = NA)
		}
		#  atc$iter <- iter
		#  atc$enh0 <- r$enh0[1] # ???
		#  atc$EnH0 <- EnH0
		#  atc$tai <- ta.lim[1]
		#  atc$tas <- ta.lim[2]
		#  atc$ti <- t1.lim[1]
		#  atc$ts <- t1.lim[2]
		#  atc$ci <- c1.lim[1]
		#  atc$cs <- c1.lim[2]
		#  atc$cote <- cote

		# combine results from each iter
		if (is.null(atcs)) {
			atcs <- atc
		} else {
			atcs <- rbind(atcs, atc)}
	}

	###########################
	a <- atcs[!is.na(atcs$n), ] # remove designs with n=NA
	a <- a[!(a$n==ceiling(a$t1*rate)), ] # remove designs with n=n1
	res <- a[order(a$enh0), ] ## find the optimal design that minimizes enh0
	des <- round(res[1, ], 4)

	Two_stage_Optimal <- data.frame(n1 = ceiling(des$t1 * rate), c1 = des$c1,
																	n = des$n, c = des$c,
																	t1 = des$t1, MTSL = des$ta + tf,
																	ES = des$enh0, PS = des$pap)

	###################################### Two stage Minmax ######################
	# if availabe, make n from optimal the maximum of n
	if(anyNA(Two_stage_Optimal$n)==1){
		maxn <- maxEn
	} else {maxn <- Two_stage_Optimal$n}

	atc0 <- data.frame(n = nsingle, t1 = tasingle, c1 = c1_p)
	atcs <- NULL
	nbpt <- nbpt_p
	pascote <- pascote_p
	cote <- range * pascote
	c1.lim <- nbpt * c(-1, 1)
	iter <- 0

	while (iter < nbmaxiter & diff(c1.lim)/nbpt > dfc1) {
		iter <- iter + 1
		cat("Calculating minmax results, iter=", iter, " & n=", round(maxn, 2),
				"& Dc1/nbpt=", round(diff(c1.lim)/nbpt, 5), "\n", sep = "")

		if (iter%%2 == 0) nbpt <- nbpt + 1
		# increase number of points to be checked within range
		cote <- cote/pascote # decrease the parameter cote to decrease search range

		# set the search range of ta, t1 and c1
		# atc0 is optimal design from last iteration
		n.lim <- atc0$n + c(-1, 1) * nsingle * cote
		ta.lim <- n.lim/rate
		t1.lim <- atc0$t1 + c(-1, 1) * tasingle * cote
		t1.lim <- pmax(0, t1.lim)
		c1.lim <- atc0$c1 + c(-1, 1) * cote

		# set the search points
		n <- seq(n.lim[1], n.lim[2], l = nbpt)
		n <- ceiling(n)
		n <- unique(n)
		ta <- n/rate
		ta <- ta[ta > 0]
		t1 <- seq(t1.lim[1], t1.lim[2], l = nbpt)
		t1 <- t1[t1 >= t1_p1 * tasingle & t1 <= t1_p2 * tasingle]
		c1 <- seq(c1.lim[1], c1.lim[2], l = nbpt)
		## set c1 > qnorm(prStop)
		if (is.null(prStop)==0){
			c1 <- c1[c1 >= qnorm(prStop)]
			# if length(c1)==0, no need to go to next iter b/c a narrower c1.lim
			# will also fail, but it might be use full to keep results from previous
			# iterations
			if(length(c1)==0) {
				warning("length(c1)==0 in iter ", iter, ". Check if prStop is too large.
  							 ", "Optimization halted, the output (Two_stage_minmax) might
  							 not be optimal.")
				minmax.err <- 1
				break
			}
		}

		# create a grid of ta, t1, and c1
		z <- expand.grid(list(ta = ta, t1 = t1, c1 = c1))
		z <- z[z$ta > z$t1, ]
		#nz <- dim(z)[1]
		z$pap <- pnorm(z$c1) # stopping prob under null
		z$eta <- z$ta - pmax(0, z$ta - z$t1) * z$pap
		z$enh0 <- z$eta * rate # expected sample size under null

		# keep design with n < maxn
		z$n <- z$ta*rate
		z <- z[z$n <= maxn, ]
		nz1 <- dim(z)[1]
		if(nz1==0) next

		### For each row of z, apply fct to get type I and II errors and c
		resz <- t(apply(z, 1, fct, ceps = ceps, alphaeps = alphaeps, nbmaxiter = nbmaxiter))
		resz <- as.data.frame(resz)
		names(resz) <- c("cL", "cU", "alphac", "betac", "rho0", "rho1", "cb1", "cb")
		r <- cbind(z, resz)
		r$c <- (r$cL + r$cL)/2             # c, criterion for final analysis
		# r$diffc <- r$cU - r$cL
		# r$diffc <- ifelse(r$diffc <= ceps, 1, 0)
		r <- r[1 - r$betac >= 1 - beta, ]  # determine if reach required power
		r <- r[order(r$n, r$enh0), ]       # order by max sample size and enh0

		# for current inter, find the combo w/ smallest n that < maxn
		if (dim(r)[1] > 0) {
			atc <- r[1, ]   # pick the row with smallest n
			if (atc$n < maxn) {
				maxn <- atc$n # update upper limit of n, maxn to the smallest n in r
				atc0 <- atc   # update atc0 to the design w/ smallest n in this iter
			}
		} else {
			atc <- data.frame(ta = NA, t1 = NA, c1 = NA, n = NA,
												pap = NA, eta = NA, enh0 = NA,cL = NA,
												cU = NA,  alphac= NA, betac = NA, rho0 = NA,
												rho1 = NA, cb1 = NA, cb  = NA, c = NA)
		}
		# atc$iter <- iter
		# atc$tai <- ta.lim[1]
		# atc$tas <- ta.lim[2]
		# atc$ti <- t1.lim[1]
		# atc$ts <- t1.lim[2]
		# atc$ci <- c1.lim[1]
		# atc$cs <- c1.lim[2]
		# atc$cote <- cote

		# combin results from each iter
		if (is.null(atcs)) {
			atcs <- atc
		} else {
			atcs <- rbind(atcs, atc)
		}
	}

	############################
	a <- atcs[!is.na(atcs$n), ] # remove designs with n=NA
	a <- a[!(a$n==ceiling(a$t1*rate)), ] # remove designs with n=n1
	res <- a[order(a$n, a$enh0), ] # find the design minimizes n and then enh0
	des <- round(res[1, ], 4)

	Two_stage_minmax <- data.frame(n1 = ceiling(des$t1 * rate), c1 = des$c1,
																 n = des$n, c = des$c,
																 t1 = des$t1, MTSL = des$ta + tf,
																 ES = des$enh0, PS = des$pap)

	############################# Two stage Admissible ###########################

	if (minmax.err ==1|optimal.err==1){
		message("Admissible design unavailabe: Incomplete optimization in minmax
						or optimal results")
		admiss.null1 <- 1
		Two_stage_Admissible=NULL

	}else if (anyNA(c(Two_stage_minmax$n, Two_stage_Optimal$n))){
		message("Admissible design unavailabe: n = NA in minmax or optimal results
  					")
		admiss.null2 <- 1
		Two_stage_Admissible=NULL

	}else if (Two_stage_Optimal$n-Two_stage_minmax$n == 0){
		message("Admissible design unavailabe: Two_stage_Optimal$n - Two_stage_min
  					max$n = 0.")
		admiss.null3 <- 1
		Two_stage_Admissible=NULL

	}else{
		nvec <- seq(Two_stage_minmax$n, Two_stage_Optimal$n, by = 1)
		for (i in 1:length(nvec)) {
			atc0 <- data.frame(n=nvec[i], t1 = tasingle, c1 = c1_p)
			atcs <- NULL
			nbpt <- nbpt_p
			pascote <- pascote_p
			cote <- range* pascote
			c1.lim <- nbpt * c(-1, 1)
			EnH0 <- maxEn
			iter <- 0

			while (iter < nbmaxiter & diff(c1.lim)/nbpt > dfc3) {
				iter <- iter + 1
				cat("Calculating admissible results, n=", nvec[i], " & EnH0=",
						round(EnH0,2), "& Dc1/nbpt=", round(diff(c1.lim)/nbpt, 5),
						"\n", sep = "")

				if (iter%%2 == 0) nbpt <- nbpt + 1
				# increase number of points to be checked within range
				cote <- cote/pascote
				# decrease the parameter cote to decrease search range

				# set search range for t1, c1
				n <- atc0$n
				t1.lim <- atc0$t1 + c(-1, 1) * tasingle * cote
				t1.lim <- pmax(0, t1.lim)
				c1.lim <- atc0$c1 + c(-1, 1) * cote
				ta <- n/rate

				# search points
				t1 <- seq(t1.lim[1], t1.lim[2], l = nbpt)
				t1 <- t1[t1 >= t1_p1 * tasingle & t1 <= t1_p2 * tasingle]
				c1 <- seq(c1.lim[1], c1.lim[2], l = nbpt)
				## set c1 > qnorm(prStop)
				if (is.null(prStop)==0){
					c1 <- c1[c1 >= qnorm(prStop)]
					# if length(c1)==0, no need to go to next iter b/c a narrower c1.lim
					# will also fail, but it might be useful to keep results from previous
					# iterations
					if(length(c1)==0) {
						warning("length(c1)==0 in iter ", iter,
										". The admissible result might not be optimized for n=",
										nvec[i])
						admiss.err <- 1
						break
					}
				}

				# Create a grid for ta, t1 and c1
				z <- expand.grid(list(ta = ta, t1 = t1, c1 = c1))
				z <- z[z$ta > z$t1, ]
				z$pap <- pnorm(z$c1)    # stopping prob under null
				z$eta <- z$ta - pmax(0, z$ta - z$t1) * z$pap
				z$enh0 <- z$eta * rate  # expected sample size under null
				z <- z[z$enh0 <= EnH0, ] # exclude combos with enh0 higher than EnH0
				nz1 <- dim(z)[1]
				if(nz1==0) next

				# For each row of z, apply fct to get type I and II errors and c
				resz <- t(apply(z, 1, fct, ceps = ceps, alphaeps = alphaeps, nbmaxiter =
													nbmaxiter))
				resz <- as.data.frame(resz)
				names(resz) <- c("cL", "cU", "alphac", "betac", "rho0", "rho1", "cb1",
												 "cb")
				r <- cbind(z, resz)
				r$n <- n
				r$c <- (r$cU + r$cL)/2             # c, criterion for final analysis
				r <- r[1 - r$betac >= 1 - beta, ]  # determine if reach required power
				r <- r[order(r$enh0), ]       # order by expected sample size under null

				# for current inter, find the combo w/ smallest enh0 that < EnH0
				if (dim(r)[1] > 0) {
					atc <- r[1, ]   # pick the row with smallest enh0
					if (atc$enh0 < EnH0) {
						EnH0 <- atc$enh0
						# update upper limit of enh0, EnH0 to the smallest enh0 in r
						atc0 <- atc  # update atc0 (ta, t1, c1, n) w/ smallest enh0
					}
				} else {
					atc <- data.frame(ta = NA, t1 = NA, c1 = NA, n = NA,
														pap = NA, eta = NA, enh0 = NA,cL = NA,
														cU = NA,  alphac= NA, betac = NA, rho0 = NA,
														rho1 = NA, cb1 = NA, cb  = NA, c = NA)
				}

				# atc$iter <- iter
				# atc$EnH0 <- EnH0
				# atc$tai <- ta.lim[1]
				# atc$tas <- ta.lim[2]
				# atc$ti <- t1.lim[1]
				# atc$ts <- t1.lim[2]
				# atc$ci <- c1.lim[1]
				# atc$cs <- c1.lim[2]
				# atc$cote <- cote
				# combine results from all inter
				if (is.null(atcs)) {
					atcs <- atc
				} else {
					atcs <- rbind(atcs, atc)
				}
			}

			a <- atcs[!is.na(atcs$n), ]
			# find the optimal design minimizes enh0 for a specific n
			res <- a[order(a$enh0), ]
			if (i==1) {des3 <- round(res[1, ], 4)
			}else des3 <- rbind(des3, round(res[1, ], 4))
		}

		# find design minimize rho
		des3$Rho <- q_value*des3$n + (1 - q_value)*des3$enh0
		des3$n1 <- ceiling(des3$t1 * rate)
		des3$MTSL <- des3$ta + tf
		des3$ES <-  des3$enh0
		des3$PS <- des3$pap
		des3 <- des3[, c("n1", "c1", "n", "c", "t1", "MTSL", "ES", "PS", "Rho")]
		des3 <- des3[order(des3$Rho), ]

		Two_stage_Admissible <-des3
	}
}



	########################## Final Results #####################################
	difn_opSg<- Two_stage_Optimal$n-Single_stage$n
	difn_opminmax <- Two_stage_Optimal$n-Two_stage_minmax$n

	param <- data.frame(shape = shape, S0 = S0, hr = hr, alpha = alpha,
											beta = beta, rate = rate, x0 = x0, tf = tf,  q_value=q_value,
											prStop=prStop, restricted=restricted)

	DESIGN <- list(param = param,
								 Single_stage = Single_stage,
								 Two_stage_Optimal = Two_stage_Optimal,
								 Two_stage_minmax=Two_stage_minmax,
								 Two_stage_Admissible=Two_stage_Admissible,
								 difn_opSg=difn_opSg,
								 difn_opminmax=difn_opminmax,
								 minmax.err = minmax.err,
								 optimal.err = optimal.err,
								 admiss.err = admiss.err,
								 admiss.null1 = admiss.null1,
								 admiss.null2 = admiss.null2,
								 admiss.null3 = admiss.null3)
	return(DESIGN)
}
