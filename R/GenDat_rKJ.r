##' Generate survival data with restricted follow-up from Weibull distribution.
##'
##' @title Generate Survival Data with Restricted Follow-Up
##' @param shape shape parameter for WB distribution.
##' @param S_fixed survival probability at the fixed time point x_fixed.
##' @param x_fixed a fixed time point where the survival probability is known.
##' @param rate a constant accrual rate.
##' @param t1 the time when interim analysis happens.
##' @param x the restricted follow-up time period, during which patients are monitored.
##' @param n sample size.
##' @param seed seed used for the random sample generation.
##' @return Two dataframes: data (the full data for final analysis) and data_t1 (the data for interim analysis), both containing the following 4 variables:\cr
##'             \emph{Entry} time when each patient enters the trial, assumed to be uniform from 0 to ta with ta as the time when enrollment ends.\cr
##'             \emph{time} the time under observation during entire trial (in data) or before t1 (in data_t1) for each patient.\cr
##'             \emph{status} the status indicator of patients (event = 1, censored = 0), which can be different between data and data_t1.\cr
##'             \emph{Total} the sum of \emph{Entry} and \emph{time}.
##' @import stats
##' @keywords internal
GenDat_rKJ <- function(shape, S_fixed, x_fixed, rate, t1, x, n, seed) {
    set.seed(seed)
    ta = n/rate
    MTSL = ta + x

    #dist=WB
    scale0 = x_fixed/(-log(S_fixed))^(1/shape)
    time = rweibull(n, shape, scale = scale0)

    # unordered
    A = runif(n, 0, ta)

    Cens = as.numeric(time < pmin(x, MTSL - A))
    # MTSL=ta + x, min(MTSL-A) = MTSL-ta=x

    Time = pmin(time, pmin(x, MTSL - A))
    # Time, the observed time capped by the smaller of restricted follow-up and MTSL-A

    #### generate Full data #######
    data=data.frame(Entry=A, time=Time, status=Cens, Total=A+Time)

    ### generate stage I data ######
    data_t1=data[(data$Entry<=t1),]
    ind_cens=which(data_t1$Total>t1)
    data_t1$status<-replace(data_t1$status, ind_cens, 0)
    data_t1$time<-replace(data_t1$time, ind_cens, pmin(t1-data_t1$Entry[ind_cens], x))

    return(list(data=data, data_t1=data_t1))
}

