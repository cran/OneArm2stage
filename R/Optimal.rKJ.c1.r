##' Optimal.rKJ() calculates the design parameters (e.g., t1, n1, n, c1, c) in the optimal two-stage design with restricted follow-up based on the one-sample log-rank (OSLR) test.
##' @title Optimal Two-Stage Design Using One-Sample Log-Rank Test with Restricted Follow-Up
##' @param dist distribution options with 'WB' as Weibull, 'GM' as Gamma, 'LN' as log-normal, 'LG' as log-logistic for the baseline hazard function. Default="WB".
##' @param shape shape parameter for the baseline hazard function assuming one of the four parametric distributions ('WB', 'GM', 'LN' and 'LG').
##' @param S0 survival probability at a fixed time point x0 under the null hypothesis.
##' @param x0 a fixed time point where the survival probability is S0 under null.
##' @param hr the hazard ratio, s1=s0^hr where s1 is the survival probability under HA and s0 is that under H0.
##' @param x the restricted follow-up time period.
##' @param rate a constant accrual rate.
##' @param alpha type I error.
##' @param beta type II error.
##' @param prStop the user-defined early stopping probability under H0, default to be NULL.
##' @return  \emph{nsignle} the required sample size for the single-stage design.\cr
##'          \emph{tasingle} the estimated accrual time for the single-stage design.\cr
##'          \emph{csingle} the critical value for the single-stage design.\cr
##'          \emph{n1} and \emph{n} required sample sizes in the two-stage design for interim and final stage, respectively.\cr
##'          \emph{c1} and \emph{c} critical values in two-stage designs for the interim and final analysis, respectively.\cr
##'          \emph{t1} the interim analysis time in the two-stage design.\cr
##'          \emph{MTSL} the maximum total study length (the sum of accrual time and restricted follow-up time).\cr
##'          \emph{ES} the expected sample size in the two-stage design.\cr
##'          \emph{PS} the probability of early stopping under null in the two-stage design.
##' @import survival
##' @export
##'
##' @examples
##' # when early stopping probability is NOT pre-defined
##' # Optimal.rKJ(dist="WB", shape=1,S0=0.20,x0=2,hr=0.569,x=2,rate=5,alpha=0.05,
##' # beta=0.2)
##' # $param
##' #   dist shape  S0    hr alpha beta rate x0 x
##' # 1   WB     1 0.2 0.569  0.05  0.2    5  2 2
##' #
##' # $Single_stage
##' #   nsingle tasingle  csingle
##' # 1      32      6.4 1.644854
##' #
##' # $Two_stage
##' #   n1     c1  n     c     t1 MTSL      ES     PS
##' # 1 21 0.0355 33 1.633 4.1882  8.6 26.7993 0.5142
##'
##' # when early stopping probability is pre-defined as 0.65
##' # Optimal.rKJ(dist="WB", shape=1,S0=0.20,x0=2,hr=0.569,x=2,rate=5,alpha=0.05,
##' # beta=0.2, prStop=0.65)
##' # $Single_stage
##' #  nsingle tasingle  csingle
##' # 1     32      6.4 1.644854
##' #
##' # $Two_stage
##' #   n1    c1  n     c     t1 MTSL      ES     PS
##' # 1 24 0.387 34 1.622 4.6959  8.8 27.1552 0.6506

##'
##'@references
##' Wu, J, Chen L, Wei J, Weiss H, Chauhan A. (2020). Two-stage phase II survival trial design. Pharmaceutical Statistics. 2020;19:214-229. https://doi.org/10.1002/pst.1983

Optimal.rKJ <- function(dist="WB", shape, S0, x0, hr, x, rate, alpha, beta, prStop=NULL){
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

  ##################################################################################
  fct <- function(zi, ceps = 1e-04, alphaeps = 1e-04, nbmaxiter = 100, dist) {
    # 49-238
    ta <- as.numeric(zi[1])
    t1 <- as.numeric(zi[2])
    c1 <- as.numeric(zi[3])

    if (dist == "WB") {
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
      ## when b=scale1(=hr) give the survival, hazard and cumulative hazard at alternative
    }

    if (dist == "LN") {
      s0 = function(u) {
        1 - plnorm(u, scale0, shape)
      }
      f0 = function(u) {
        dlnorm(u, scale0, shape)
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
      scale0 = log(x0) - shape * qnorm(1 - S0)
      scale1 = hr
    }

    if (dist == "LG") {
      s0 = function(u) {
        1/(1 + (u/scale0)^shape)
      }
      f0 = function(u) {
        (shape/scale0) * (u/scale0)^(shape - 1)/(1 + (u/scale0)^shape)^2
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
      scale0 = x0/(1/S0 - 1)^(1/shape)
      scale1 = hr
    }

    if (dist == "GM") {
      s0 = function(u) {
        1 - pgamma(u, shape, scale0)
      }
      f0 = function(u) {
        dgamma(u, shape, scale0)
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
      root0 = function(t) {
        1 - pgamma(x0, shape, t) - S0
      }
      scale0 = uniroot(root0, c(0, 10))$root
      scale1 = hr
    }

    ####################################################################################
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
    p0 = integrate(g0, 0, x)$value
    p1 = integrate(g1, 0, x)$value
    p00 = integrate(g00, 0, x)$value
    p01 = integrate(g01, 0, x)$value
    sigma2.1 = p1 - p1^2 + 2 * p00 - p0^2 - 2 * p01 + 2 * p0 * p1 # exact var(W) under H1
    sigma2.0 = p0
    om = p0 - p1

    ####################################################################################
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
    p0 = integrate(g0, 0, x)$value
    p1 = integrate(g1, 0, x)$value
    p00 = integrate(g00, 0, x)$value
    p01 = integrate(g01, 0, x)$value
    sigma2.11 = p1 - p1^2 + 2 * p00 - p0^2 - 2 * p01 + 2 * p0 * p1  # exact var(W1) under H1
    sigma2.01 = p0
    om1 = p0 - p1

    ####################################################################################
    q1 = function(t) {
      s0(t) * h0(t) * G1(t)
    }
    q = function(t) {
      s0(t) * h0(t)
    }
    v1 = integrate(q1, 0, x)$value
    v = integrate(q, 0, x)$value
    rho0 = sqrt(v1/v)                 # correlation between W1 and W under H0?
    rho1 <- sqrt(sigma2.11/sigma2.1) # correlation between W1 and W under H1?

    ################################### Calculate Type I error alpha #################################################
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

    ################################### Calculate power #################################################

    cb1 <- sqrt((sigma2.01/sigma2.11)) * (c1 - (om1 * sqrt(rate * t1)/sqrt(sigma2.01))) # rate*t1=n1
    cb <- sqrt((sigma2.0/sigma2.1)) * (c - (om * sqrt(rate * ta))/sqrt(sigma2.0)) # rate*ta=n
    pwrc <- calculate_power(cb = cb, cb1 = cb1, rho1 = rho1)

    res <- c(cL, cU, alphac, 1 - pwrc, rho0, rho1, cb1, cb)
    return(res)
  }

  ####################################################################################
  ####################################################################################

  c1 <- 0
  rho0 <- 0
  cb1 <- 0
  rho1 <- 0
  hz <- c(0, 0)
  ceps <- 0.001
  alphaeps <- 0.001
  nbmaxiter <- 100

  if (dist == "WB") {
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
  }

  if (dist == "LN") {
    s0 = function(u) {
      1 - plnorm(u, scale0, shape)
    }
    f0 = function(u) {
      dlnorm(u, scale0, shape)
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
    scale0 = log(x0) - shape * qnorm(1 - S0)
    scale1 = hr
  }

  if (dist == "LG") {
    s0 = function(u) {
      1/(1 + (u/scale0)^shape)
    }
    f0 = function(u) {
      (shape/scale0) * (u/scale0)^(shape - 1)/(1 + (u/scale0)^shape)^2
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
    scale0 = x0/(1/S0 - 1)^(1/shape)
    scale1 = hr
  }

  if (dist == "GM") {
    s0 = function(u) {
      1 - pgamma(u, shape, scale0)
    }
    f0 = function(u) {
      dgamma(u, shape, scale0)
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
    root0 = function(t) {
      1 - pgamma(x0, shape, t) - S0
    }
    scale0 = uniroot(root0, c(0, 10))$root
    scale1 = hr
  }

  ##################################### Single stage ###############################################
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
  p0 = integrate(g0, 0, x)$value    # eq(3)
  p1 = integrate(g1, 0, x)$value    # eq(4)
  p00 = integrate(g00, 0, x)$value  # eq(5)
  p01 = integrate(g01, 0, x)$value  # eq(6)
  s1 = sqrt(p1 - p1^2 + 2 * p00 - p0^2 - 2 * p01 + 2 * p0 * p1)
  s0 = sqrt(p0)
  om = p0 - p1 # omega in eq(2)

  #############################

  nsingle <- (s0 * qnorm(1 - alpha) + s1 * qnorm(1 - beta))^2/om^2 # eq(2)
  nsingle <- ceiling(nsingle)
  tasingle <- nsingle/rate

  Single_stage <- data.frame(nsingle = nsingle, tasingle = tasingle, csingle = qnorm(1 - alpha))
  atc0 <- data.frame(n = nsingle, t1 = tasingle, c1 = 0.25)
  nbpt <- 11
  pascote <- 1.26
  cote <- 1 * pascote
  c1.lim <- nbpt * c(-1, 1)
  EnH0 <- 10000
  iter <- 0

  ###################################### Two stage ##############################################

  while (iter < nbmaxiter & diff(c1.lim)/nbpt > 0.001) {
    iter <- iter + 1
    cat("iter=", iter, "&EnH0=", round(EnH0, 2), "& Dc1/nbpt=", round(diff(c1.lim)/nbpt, 5), "\n", sep = "")

    if (iter%%2 == 0)
      nbpt <- nbpt + 1                             # increase number of points to be checked within range

    cote <- cote/pascote                           # decrease the parameter cote to decrease search range
    n.lim <- atc0$n + c(-1, 1) * nsingle * cote    # atc0 is optimal design from last iteration
    t1.lim <- atc0$t1 + c(-1, 1) * tasingle * cote

    #if(iter==1){c1.lim <- atc0$c1 + c(-2.1, 2.1) * cote} else
    c1.lim <- atc0$c1 + c(-1, 1) * cote       ## range of c1

    ta.lim <- n.lim/rate
    t1.lim <- pmax(0, t1.lim)

    n <- seq(n.lim[1], n.lim[2], l = nbpt)
    n <- ceiling(n)
    n <- unique(n)
    ta <- n/rate
    t1 <- seq(t1.lim[1], t1.lim[2], l = nbpt)
    c1 <- seq(c1.lim[1], c1.lim[2], l = nbpt)

    ## set c1 > qnorm(prStop)
    if (is.null(prStop)==0){

        c1 <- c1[c1 >= qnorm(prStop)]

        if(length(c1)==0) {

        stop("length(c1)==0, use prStop <= 0.89 ")

        }
    }


    ta <- ta[ta > 0]
    t1 <- t1[t1 >= 0.2 * tasingle & t1 <= 1.2 * tasingle]

    z <- expand.grid(list(ta = ta, t1 = t1, c1 = c1))
    z <- z[z$ta > z$t1, ]
    nz <- dim(z)[1]
    z$pap <- pnorm(z$c1)                          ## stopping prob under null
    z$eta <- z$ta - pmax(0, z$ta - z$t1) * z$pap
    z$enh0 <- z$eta * rate                        ## expected sample size under null
    z <- z[z$enh0 <= EnH0, ]
    nz1 <- dim(z)[1]

    ### For each row of z, apply fct
    resz <- t(apply(z, 1, fct, ceps = ceps, alphaeps = alphaeps, nbmaxiter = nbmaxiter, dist = dist))

    resz <- as.data.frame(resz)
    names(resz) <- c("cL", "cU", "alphac", "betac", "rho0", "rho1", "cb1", "cb")
    r <- cbind(z, resz)
    r$n <- r$ta * rate
    r$c <- r$cL + (r$cU - resz$cL)/2             # c, criterion for final analysis
    r$diffc <- r$cU - r$cL
    r$diffc <- ifelse(r$diffc <= ceps, 1, 0)     #
    r <- r[1 - r$betac >= 1 - beta, ]            ## determine if reaches required power
    r <- r[order(r$enh0), ]                      # order by expected sample size under null

    if (dim(r)[1] > 0) {
      #atc <- r[, c("ta", "t1", "c1", "n")][1, ]  # pick the row with smallest enh0
      atc <- r[, ][1, ]

      r1 <- r[1, ]
      if (r1$enh0 < EnH0) {
        EnH0 <- r1$enh0                         # update upper limit of enh0, EnH0 to the smallest enh0 in r
        atc0 <- atc                             # update atc0, the (ta, t1, c1, n) w/ smallest enh0 in this iteration
      }
    } else {
      atc <- data.frame(ta = NA, t1 = NA, c1 = NA, n = NA,
                        pap = NA, eta = NA, enh0 = NA,cL = NA,
                        cU = NA,  alphac= NA, betac = NA, rho0 = NA,
                        rho1 = NA, cb1 = NA, cb  = NA, c = NA, diffc= NA)
    }

    atc$iter <- iter
    atc$enh0 <- r$enh0[1]
    atc$EnH0 <- EnH0
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

  ####################################################################################
  a <- atcs[!is.na(atcs$n), ]

  res <- a[order(a$enh0), ] ## find the optimal design minimizes enh0
  des <- round(res[1, ], 4)

  if (is.null(prStop)==0){
  param <- data.frame(dist= dist, shape = shape, S0 = S0, hr = hr, alpha = alpha,
                      beta = beta, rate = rate, x0 = x0, x = x,  prStop=prStop)
  } else{
  param <- data.frame(dist= dist, shape = shape, S0 = S0, hr = hr, alpha = alpha,
                        beta = beta, rate = rate, x0 = x0, x = x)
   }

  Two_stage <- data.frame(n1 = ceiling(des$t1 * param$rate), c1 = des$c1,
                          n = des$n, c = des$c,
                          t1 = des$t1, MTSL = des$ta + param$x,
                          ES = des$enh0, PS = des$pap)
  DESIGN <- list(param = param, Single_stage = Single_stage, Two_stage = Two_stage)
  return(DESIGN)

}


