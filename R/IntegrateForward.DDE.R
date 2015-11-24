
dderhs <- function(t, r, parms){
    tau = parms$tau
    beta = parms$beta
    proc = parms$proc
    pars = parms$pars
    ndelay <- parms$ndelay
    fdobj0 <- parms$fdobj0
    fdobj.d <- parms$fdobj.d
    t0 <- parms$t0
    tau.max <- parms$tau.max
    if(t - t0 < tau.max){
        delayObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = t, tau = tau, beta= beta, ndelay = ndelay, lik = TRUE)
        proc$more$y.d <- delayObj$y.d
    }
    else{
        y.d <- rep(0, length(ndelay))
        i <- 0
        for(idelay in ndelay){
            i <- i + 1
            for(j in 1:length(tau[[i]])){
                y.d[i] <- y.d[i] + beta[j] * deSolve::lagvalue(t - tau[[i]][j], idelay)
            }
        }
        proc$more$y.d <- matrix(y.d, length(ndelay),1)
    }
    r = matrix(r, 1, length(r))
    if (!is.null(proc$names))
        colnames(r) = proc$names
    y = as.vector(proc$fn(t, r, pars, proc$more))
    return(list(y))
}

##' Solves a delay differential equation going forward based on a proc object.
##' @title Intergrate Forward a DDE Model
##' @param ts A time series at which the state of the process is predicted
##' @param pars Estimated parameters.
##' @param beta Estimated contributions of lags of delay.
##' @param proc The \code{proc} object returned from estimation functions.
##' @param more An object specifying additional arguments to fn.
##' @param tau A list of delay lags.
##' @param ndelay A vector inidicating which state process has a delay term.
##' @param fdobj0  A functional data object fitted with the initial history part of the data.
##' @param fdobj.d A functional data object fitted by generalized profiling.
##' @return A list of objects
##' \describe{
##' \item{times}{The time points at where the predictions are made.}
##' \item{states}{The predicted states of the process.}
##' }
##' @export
##' @author Ziqian Zhou
IntegrateForward.dde <- function(ts, pars, beta,  proc, more, tau, ndelay, fdobj0, fdobj.d){
    proc = proc$more
    tau.max <- 0
    for(i in 1:length(tau)){
        tau[[i]] = tau[[i]][beta > 0]
        tau.max <- max(tau.max, max(tau[[i]]))
    }
    beta = beta[beta > 0]
    t0 <- ts[1]
    y0 <- eval.fd(t0, fdobj.d)
    parms = list(pars = pars, proc = proc, more = more, fdobj0 = fdobj0, fdobj.d = fdobj.d, beta = beta, tau = tau, ndelay = ndelay, tau.max = tau.max, t0 = t0)
    out = deSolve::dede(y = y0, times = ts, func = dderhs, parms = parms)
    return(list(times = out[, 1], states = out[, 2:ncol(out)]))

}
