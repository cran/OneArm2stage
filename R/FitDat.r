##' The function fits parametric models with the underlying distributions assumed to be Weibull.
##' @title Fit Historical Survival Data Assuming the Failure Time Follows the Weibull
##'        Distribution
##' @param data a historical survival data sample, has to contain two variables 'Time' and 'Cens': \cr
##'             \emph{Time}, time under observation during trial for each patient.\cr
##'             \emph{Cens}, the status indicator of patients (event = 1, censored = 0).
##' @return \emph{fit.Weibull} Fitted models assuming Weibull distributions.\cr
##'         \emph{AIC} AIC values from the fitted model. \cr
##'         \emph{parameter.estimates} the estimated parameters from the fitted model.
##' @export
##' @import flexsurv IPDfromKM
##'
##' @examples
##' library(IPDfromKM)
##' # a sample dataset that we already extracted from Wang et al, 2018.
##' df<- read.csv(system.file("extdata", "df.csv", package = "OneArm2stage"))
##'
##' # risk time points
##' trisk <- c(0,2,4,6,8,10,12,14,16,18,20,22,24)
##'
##' # number of patients at risk at each risk time point
##' nrisk.radio <- c(124,120,115,110,107,104,103,95,46,18,11,8,0)
##'
##' # Preprocess the raw coordinates into an proper format for reconstruct IPD
##' pre_radio <- preprocess(dat=df, trisk=trisk,
##'                      nrisk=nrisk.radio,totalpts=NULL,maxy=100)
##'
##' #Reconstruct IPD
##' est_radio <- getIPD(prep=pre_radio,armID=0,tot.events=NULL)
##'
##' # shift the IPD data into the proper format for 'FitDat()'
##' ipd <- est_radio$IPD
##' dat3 <- as.data.frame(cbind(rep(0, nrow(ipd)),ipd$time, ipd$status))
##' colnames(dat3) <- c("Entry", "Time", "Cens")
##'
##' # use FitDat function to fit the historical dat
##' modelSelect <- FitDat(dat3)
##' modelSelect$AIC
##' # Weibull
##' # 301.7776
##'
##' # check the estimated parameters from the modeling results
##' modelSelect$parameter.estimates
##' # $Weibull
##' # shape     scale
##' # 0.1133671 3.9939753
##'
##'@references
##' Wang, M., Rule, S., Zinzani, P. L., Goy, A., Casasnovas, O., Smith, S. D.,..., Robak, T. (2018).
##' Acalabrutinib in relapsed or refractory mantle cell lymphoma (ACE-LY-004):
##' a single-arm, multicentre, phase 2 trial. The Lancet, 391(10121), 659–667.
##' https://doi.org/10.1016/s0140-6736(17)33108-2

FitDat <- function(data) {
    fit.w <- flexsurvreg(Surv(Time, Cens) ~ 1, data = data, dist = "weibull")  # use indep var similar to target population?

    AICvec <- c(fit.w$AIC)
    names(AICvec) <- c("Weibull")
    parameter <- list(fit.w$coefficients)
    names(parameter) <- c("Weibull")

    res <- list(fit.Weibull = fit.w, AIC = AICvec, parameter.estimates = parameter)

    return(res)
}
