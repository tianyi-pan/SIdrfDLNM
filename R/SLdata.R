#' Prepare Super Learner quantities for a single-index DRF-DLNM
#'
#' @param object An object of class \code{sidDLNM_fit}.
#' @param dat.c Required counterfactual data, including the full exposure
#'   histories, on the predictor scale expected by \code{predict()}. Retain
#'   the fitted reference records and model structure.
#' @param verbose Whether to print progress.
#' @param AIC.only Prepare singleton AIC quantities for \code{mDLNM::AICdata}
#'   instead of the Super Learner predicted means, default \code{FALSE}.
#' @param nthreads Number of threads for deletion calculations, default \code{1}.
#'   Uses OpenMP for parallel calculation.
#' @param ... Unused.
#' @exportS3Method mDLNM::SLdata
SLdata.sidDLNM_fit <- function(object, dat.c, verbose = FALSE,
                              AIC.only = FALSE, nthreads = 1, ...) {
  if(identical(object$family, "gaussian")) stop("SLdata() requires a negative-binomial fit.")
  if(missing(dat.c) || is.null(dat.c)) stop("dat.c is required.")
  if(length(list(...))) stop("Unused arguments.")

  dat <- object$data
  point <- object$point
  if(is.null(dat$Bindex)) stop("Refit this model to store Bindex before using SLdata().")

  coef.samples <- NULL
  if(AIC.only) {
    samples <- object$CI.sample
    if(is.null(samples$alpha_f_sample) || nrow(samples$alpha_f_sample) < 2) {
      stop("AICdata() requires at least two stored coefficient samples.")
    }
    index.coef <- t(solve(dat$Bindex, t(samples$index_sample)))
    index.coef <- index.coef[, -1, drop = FALSE] / index.coef[, 1]
    coef.samples <- cbind(samples$alpha_f_sample, index.coef,
                          samples$betaR_sample, samples$betaF_sample)
  }

  mod.build <- sidDLNMbuild(dat$y, dat$B_inner, dat$knots_f, dat$knots_w,
                            dat$Sw, dat$SwR, dat$Sf, dat$SfR,
                            dat$Xrand, dat$Xfix, dat$Zf, dat$offset$Xoffset,
                            dat$r, dat$Bindex, point$alpha_f, point$delta,
                            point$log_theta, point$log_smoothing_f, point$log_smoothing_w,
                            point$betaR, point$betaF, point$log_smoothing)
  mod.c <- stats::predict(object, dat = dat.c, .build.only = TRUE)
  result <- SLdatasidDLNM(mod.build$address.eigen, mod.c$address.eigen, verbose,
                          AIC.only, nthreads, coef.samples)
  theta <- exp(point$log_theta)

  if(AIC.only) {
    aic.names <- c("AIC_conventional", "AIC_cAIC")
    colnames(result$aic.loo) <- aic.names
    return(list(y = dat$y, row.id = rownames(object$modeldata), theta = theta,
                aic = stats::setNames(as.vector(result$aic), aic.names),
                aic.loo = result$aic.loo,
                delta = result$delta, delta.loo = as.vector(result$delta.loo),
                variance = stats::var(as.vector(result$delta.sample)),
                estimand = "ratio"))
  }

  return(list(y = dat$y,
              row.id = rownames(object$modeldata),
              theta = theta,
              mu.loo = as.vector(result$mu.loo),
              mu.lto = result$mu.lto,
              delta = result$delta,
              delta.loo = as.vector(result$delta.loo),
              estimand = "ratio"))
}
