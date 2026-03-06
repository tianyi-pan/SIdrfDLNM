#' Compute Rate Ratio
#'
#' @param object An object.
#' @param ... Passed to methods.
#' @export
RateRatio <- function(object, ...) {
  UseMethod("RateRatio")
}


#' Title Rate Ratio for sidDLNM_fit objects
#'
#' @param object object of class \code{sidDLNM_fit}.
#' @param verbose whether to print messages during the process, default \code{FALSE}.
#'
#' @return
#' @importFrom mgcv s
#' @export
RateRatio.sidDLNM_fit <- function(object, x0, x1, verbose = FALSE, ...) {


  pc <- object$data$pc
  alpha_f <- object$point$alpha_f
  kw <- object$data$kw
  kE <- object$data$kE

  maxL <- object$data$maxL

  Blag <- sapply(seq(0, maxL), function(l0) Bsplinevec2(l0, object$data$knots_w, 4))

  cen <- object$data$knots_f[as.integer(kE/2) + 2] # does not affect the contrasts
  # cen <- 0

  ## point estimate

  B_inner_sum_x0 <- Reduce("+", mapply("*", x0, object$point$index, SIMPLIFY = FALSE))
  B_inner_sum_x1 <- Reduce("+", mapply("*", x1, object$point$index, SIMPLIFY = FALSE))


  gridEl.x0 <- data.frame(E = B_inner_sum_x0,
                           l = 0:maxL)
  gridEl.x1 <- data.frame(E = B_inner_sum_x1,
                          l = 0:maxL)

  eta.E.0 <- sum(apply(gridEl.x0, 1, function(row.) {
    SurfaceEval(row.[1], cen, row.[2], alpha_f,
                          object$data$knots_f, object$data$Zf, Blag)
  }))

  eta.E.1 <- sum(apply(gridEl.x1, 1, function(row.) {
    SurfaceEval(row.[1], cen, row.[2], alpha_f,
                          object$data$knots_f, object$data$Zf, Blag)
  }))

  RR.est <- data.frame(eta.E.0 = eta.E.0,
                       eta.E.1 = eta.E.1,
                       RR = exp(eta.E.1 - eta.E.0))

  ## CI
  RR.sample <- lapply(1:nrow(object$CI.sample$index_sample), function(i) {
    if(verbose) {
      if(i %% 100 == 0) {
        cat("Sampling iteration:", i, " / ", nrow(object$CI.sample$index_sample), "\n")
      }
    }

    B_inner_sum_x0 <- Reduce("+", mapply("*", x0, object$CI.sample$index_sample[i,], SIMPLIFY = FALSE))
    B_inner_sum_x1 <- Reduce("+", mapply("*", x1, object$CI.sample$index_sample[i,], SIMPLIFY = FALSE))


    gridEl.x0 <- data.frame(E = B_inner_sum_x0,
                            l = 0:maxL)
    gridEl.x1 <- data.frame(E = B_inner_sum_x1,
                            l = 0:maxL)

    eta.E.0 <- sum(apply(gridEl.x0, 1, function(row.) {
      SurfaceEval(row.[1], cen, row.[2], object$CI.sample$alpha_f_sample[i,],
                            object$data$knots_f, object$data$Zf, Blag)
    }))

    eta.E.1 <- sum(apply(gridEl.x1, 1, function(row.) {
      SurfaceEval(row.[1], cen, row.[2], object$CI.sample$alpha_f_sample[i,],
                            object$data$knots_f, object$data$Zf, Blag)
    }))

    return(data.frame(eta.E.0 = eta.E.0,
                      eta.E.1 = eta.E.1,
                      RR = exp(eta.E.1 - eta.E.0)))
  })


  RR.sample <- data.table::rbindlist(RR.sample)
  return(list(RR.est = RR.est,
              RR.sample = RR.sample)
  )
}
