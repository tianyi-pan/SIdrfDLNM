#' Title summary() Method for Objects of Class 'sidDLNM_fit'
#'
#' @param object object of class \code{aceDLNM_fit}.
#' @param E.eval a vector of ACE where the ACE response-exposure-function is
#'   evaluated. The default is a sequence from 0.0025 quantile to 0.9975
#'   quantile of estimated E, if the argument is missing.
#' @param others.eval a data frame containing variables where the other smooth
#'   functions are evaluated.
#' @param plot whether or not to show the plots, default \code{FALSE}.
#' @param true.function the list containing the true functions.
#'   f(E) - f(E0).
#' @param ...
#'
#' @import dplyr
#' @return a list containing (1) the estimated functions evaluated at
#'   \code{E.eval} and \code{others.eval} or the default values, (2) AIC, and
#'   (3) plots.
#' @export
summary.sidDLNM_fit <- function(object, E.eval, l.eval, others.eval = NULL,
                                cen,
                                plot = FALSE, true.function, ...){
  ## point estimate
  pc <- object$data$pc
  alpha_f <- object$point$alpha_f
  kw <- object$data$kw
  kE <- object$data$kE

  if(missingArg(cen)) {
    cen <- 3
    # cat("set centering value as 3")
  }
  B_inner_sum <- Reduce("+", mapply("*", object$data$B_inner, object$point$index, SIMPLIFY = FALSE))


  Q1 <- quantile(B_inner_sum[,1], 0.25)
  Q3 <- quantile(B_inner_sum[,1], 0.75)

  maxL <- object$data$maxL
  E.mode.quantile <- quantile(B_inner_sum, c(0.0025, 0.9975))
  if(missingArg(E.eval)) E.eval <- seq(E.mode.quantile[1], E.mode.quantile[2], length.out = 100)
  if(missingArg(l.eval)) l.eval <- seq(0, maxL)

  Blag <- sapply(seq(0, maxL), function(l0) sidDLNM:::Bsplinevec2(l0, object$data$knots_w, 4))

  if(!(Q1 %in% E.eval)) E.eval <- c(Q1, E.eval)
  if(!(Q3 %in% E.eval)) E.eval <- c(Q3, E.eval)
  E.eval <- sort(E.eval)


  # sidDLNM:::SurfaceEval(E0, cen, l0, alpha_f,
  #             object$data$knots_f, object$data$Zf, Blag)
  # (basisfunEl(E0, l) - basisfunElcen(l0)) %*% alpha_f

  # basisfunEl <- function(E0, l0) {
  #   sidDLNM:::Bsplinevec2Con(E0, object$data$knots_f, 4, object$data$Zf) %x% sidDLNM:::Bsplinevec2(l0, object$data$knots_w, 4)
  # }
  # basisfunElcen <- function(l0) {
  #   sidDLNM:::Bsplinevec2Con(cen, object$data$knots_f, 4, object$data$Zf) %x% sidDLNM:::Bsplinevec2(l0, object$data$knots_w, 4)
  # }

  gridEl <- expand.grid(E = E.eval,
                        l = l.eval)

  surface.mode <- apply(gridEl, 1, function(row.) {
    sidDLNM:::SurfaceEval(row.[1], cen, row.[2], alpha_f,
                          object$data$knots_f, object$data$Zf, Blag)
  })

  surface.est <- data.frame(gridEl)
  surface.est$mode <- surface.mode

  # filter(surface.est, l == 1) %>% ggplot(aes(x = E, y = mode)) + geom_line()
#
#   library(plotly)
#   p <- plot_ly(
#     surface.est, x= ~E, y= ~l, z= ~mode,
#     type='mesh3d', intensity = ~mode,
#     colors= colorRamp(rainbow(5))
#   )
#   p


  ## CI
  R.CI <- nrow(object$CI.sample[[1]])

  ## implemented in R
  #   surface.sample <- lapply(1:R.CI, function(i){
  #     alpha_f_sample <- object$CI.sample$alpha_f_sample[i,]
  #
  #     surface.i <- apply(gridEl, 1, function(row.) {
  #       # (basisfunEl(row.[1], row.[2]) - basisfunElcen(row.[2])) %*% alpha_f_sample
  #       sidDLNM:::SurfaceEval(row.[1], cen, row.[2], alpha_f_sample,
  #                             object$data$knots_f, object$data$Zf, Blag)
  #     })
  #
  #     return(cbind(gridEl, surface.i, rep(i, length(surface.i))))
  #   })
  # surface.sample <- do.call(rbind.data.frame, surface.sample)
  # colnames(surface.sample) <- c("E","l","est", "iter")
  # surface.CI <- group_by(surface.sample, E, l) %>%
  #   summarize(ll = quantile(est, 0.025), ul = quantile(est, 0.975))

  ## faster implementation in cpp
  surface.sample <- sidDLNM:::SurfaceCI(as.matrix(gridEl), object$CI.sample$alpha_f_sample, cen,
                        object$data$knots_f, object$data$Zf, Blag)


  surface.est <- cbind(surface.est, t(apply(surface.sample, 1, function(row.) quantile(row., c(0.025, 0.975)))))
  colnames(surface.est)[c(4,5)] <- c("ll", "ul")

  surface.sample <- cbind(gridEl, surface.sample)

  out <- list(samples = list(surface = surface.sample),
              est = list(surface = surface.est))



  if(!missingArg(true.function)) {
    if(!is.null(true.function$surface)){
      surface.true <- apply(gridEl, 1, function(row.) {
        true.function$surface(row.[1], row.[2]) - true.function$surface(cen, row.[2])
      })
      out$est$surface$true <- surface.true
    }
  }
  out$data$Q1 <- Q1
  out$data$Q3 <- Q3
  out$data$E.eval <- E.eval
  out$data$l.eval <- l.eval


  ### Other smooth terms
  if(!is.null(object$smooth$others)) {
    # others.eval
    SS <- object$smooth$others

    smooth.name <- unlist(lapply(SS, "[[", "term")) # var name of smooth terms
    if(is.null(others.eval)) {
      others.eval <- data.frame(lapply(object$modeldata[,smooth.name, drop=FALSE],
                                       function(col){
                                         # support random effects
                                         if(class(col) == "factor") return(factor(rep(col[1], 500), levels = levels(col)))
                                         else return(seq(min(col), max(col), length.out = 500))
                                       }))
    }
    Xpred <- lapply(SS, function(a){
      preddat <- data.frame(others.eval[,a$term])
      colnames(preddat) <- a$term
      mgcv::PredictMat(a, data = preddat)
    })
    Xpred <- Reduce(cbind, Xpred)
    Xrandpred <- as.matrix(Xpred %*% object$data$UR %*% object$data$Dpi) ## reparametrized
    Xfixpred <- as.matrix(Xpred %*% object$data$UF)
    Xfixpred <- Xfixpred[,apply(Xfixpred, 2, function(col.) sum(col. != 0) > 0), drop = FALSE]


    betaF_R <- object$point$betaF[-seq_len(ncol(object$data$Xfix) - sum(object$data$m))] # unpenalized terms
    splitR_index <- c() # index for betaR
    for (i in seq_along(object$data$r)) {
      splitR_index <- c(splitR_index, rep(i, object$data$r[i]))
    }
    splitF_index <- c() # index for betaF
    for (i in seq_along(object$data$m)) {
      splitF_index <- c(splitF_index, rep(i, object$data$m[i]))
    }


    Xrandpred_split <- lapply(split(seq_len(ncol(Xrandpred)), splitR_index),
                              function(j) Xrandpred[,j])
    betaR_split <- split(object$point$betaR, splitR_index)

    Xfixpred_split <- lapply(split(seq_len(ncol(Xfixpred)), splitF_index),
                             function(j) Xfixpred[,j])
    betaF_split <- split(betaF_R, splitF_index)


    ## obtain fitted value
    fittedR <- mapply(function(a,b) as.matrix(a) %*% as.vector(b),
                      Xrandpred_split, betaR_split, SIMPLIFY = FALSE)
    fittedF <- mapply(function(a,b) as.matrix(a) %*% as.vector(b),
                      Xfixpred_split, betaF_split, SIMPLIFY = FALSE)


    if(length(fittedF) == 0) {
      ## if there is no unpenalized component in smooth
      fitted <- fittedR
    } else {
      insertname <- setdiff(names(fittedR), names(fittedF)) # fill the fittedF to match the length
      if(length(insertname) != 0) {
        for (jj in insertname) fittedF[[jj]] <- matrix(0, nrow = nrow(fittedR[[1]]), ncol = 1)
      }
      ## order by name
      fittedF <- fittedF[order(names(fittedF))]

      fitted <- mapply(function(a,b) a+b, fittedR, fittedF, SIMPLIFY = FALSE)
    }

    fitted <- do.call("rbind", fitted)

    smooth.mode <- data.frame(x = as.matrix(as.vector(sapply(others.eval, as.numeric)), ncol = 1),
                              var = rep(colnames(others.eval), each = nrow(others.eval)),
                              mode = fitted)

    fitted.CI <- lapply(1:R.CI, function(i){
      betaF_sample <- object$CI.sample$betaF_sample[i,]
      betaR_sample <- object$CI.sample$betaR_sample[i,]
      betaR_sample_split <- split(betaR_sample, splitR_index)
      betaF_sample_split <- split(betaF_sample[-seq_len(ncol(object$data$Xfix) - sum(object$data$m))],
                                  splitF_index)

      fittedR_sample <- mapply(function(a,b) as.matrix(a) %*% as.vector(b),
                               Xrandpred_split, betaR_sample_split, SIMPLIFY = FALSE)
      fittedF_sample <- mapply(function(a,b) as.matrix(a) %*% as.vector(b),
                               Xfixpred_split, betaF_sample_split, SIMPLIFY = FALSE)
      if(length(fittedF_sample) == 0) {
        ## if there is no unpenalized component in smooth
        fitted_sample <- fittedR_sample
      } else {
        insertname <- setdiff(names(fittedR_sample), names(fittedF_sample)) # fill the fittedF to match the length
        if(length(insertname) != 0) {
          for (jj in insertname) fittedF_sample[[jj]] <- matrix(0, nrow = nrow(fittedR_sample[[1]]), ncol = 1)
        }
        ## order by name
        fittedF_sample <- fittedF_sample[order(names(fittedF_sample))]

        fitted_sample <- mapply(function(a,b) a+b, fittedR_sample, fittedF_sample, SIMPLIFY = FALSE)
      }
      fitted_sample <- do.call("rbind", fitted_sample)

      return(data.frame(x = as.matrix(as.vector(sapply(others.eval, as.numeric)), ncol = 1),
                        var = rep(smooth.name, each = nrow(others.eval)),
                        est = fitted_sample,
                        iter = rep(i, length(fitted_sample))))
    })
    fitted.CI <- do.call(rbind.data.frame, fitted.CI)

    ## get quantile
    smooth.df <- group_by(fitted.CI, x, var) %>%
      summarize(ll = quantile(est, 0.025), ul = quantile(est, 0.975))

    smooth.df <- left_join(smooth.mode, smooth.df, by = c("x","var"))

    if(!missingArg(true.function)) {
      if(!is.null(true.function$smooth)){
        true.function$smooth <- Vectorize(true.function$smooth)
        smooth.df$true <- true.function$smooth(smooth.df$x, smooth.df$var)

        # if(center) {
        #   modemean <- group_by(smooth.df, var) %>% summarize(modemean = mean(mode)) %>% select(modemean) %>% as.numeric()
        #   truemean <- group_by(smooth.df, var) %>% summarize(truemean = mean(true)) %>% select(modemean) %>% as.numeric()
        #
        #   smooth.df$true <- smooth.df$true - truemean
        #   smooth.df$mode <- smooth.df$mode - modemean
        #   smooth.df$ul <- smooth.df$ul - modemean
        #   smooth.df$ll <- smooth.df$ll - modemean
        # }
      }
    }

    # remove re terms
    re_id <- which(unlist(lapply(lapply(SS, attr, "class"), "[[", 1)) == "random.effect")

    if(length(re_id) >= 1) smooth.df <- dplyr::filter(smooth.df, !var %in% unique(smooth.df$var)[re_id])

    out$est$smooth = smooth.df

  }

  if(plot) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))

    ## estimated surface
    par(mar = c(2, 2, 2, 2))
    with(out$est$surface,
         persp(E.eval, l.eval, matrix(mode, nrow = length(E.eval), ncol = length(l.eval)),
               xlab = "Exposure", ylab = "Lag", zlab = "", # main = "Estimated Surface",
               phi=25,cex.axis=0.7,cex.lab=0.8,nticks=9,
               col = "lightblue",
               cex.main = 0.8,
               theta=230)
         )
    par(mar = c(5, 4, 4, 2)) # reset to default

    ## lag curve
    df.plot <- filter(out$est$surface, E == Q1)
    with(df.plot, plot(l, mode, type = "l", ylim = c(min(c(ll)), max(c(ul))), ylab = "est", xlab = "lag",
                       main = paste0("Exposure = ", round(as.numeric(Q1),1), " (25% quantile)")))
    with(df.plot, lines(l, ll, lty = "dashed"))
    with(df.plot, lines(l, ul, lty = "dashed"))
    if(!missingArg(true.function)) with(df.plot, lines(l, true, col = "blue"))

    df.plot <- filter(out$est$surface, E == Q3)
    with(df.plot, plot(l, mode, type = "l", ylim = c(min(c(ll)), max(c(ul))),
                       ylab = "est", xlab = "lag",
                       main = paste0("Exposure = ", round(as.numeric(Q3),1), " (75% quantile)")))
    with(df.plot, lines(l, ll, lty = "dashed"))
    with(df.plot, lines(l, ul, lty = "dashed"))
    if(!missingArg(true.function)) with(df.plot, lines(l, true, col = "blue"))


    ## exposure curve
    df.plot <- filter(out$est$surface, l == 0)
    with(df.plot, plot(E, mode, type = "l", ylim = c(min(c(ll)), max(c(ul))), ylab = "est", xlab = "lag",
                       main = paste0("Lag = 0")))
    with(df.plot, lines(E, ll, lty = "dashed"))
    with(df.plot, lines(E, ul, lty = "dashed"))
    if(!missingArg(true.function)) with(df.plot, lines(E, true, col = "blue"))

    df.plot <- filter(out$est$surface, l == 2)
    with(df.plot, plot(E, mode, type = "l", ylim = c(min(c(ll)), max(c(ul))), ylab = "est", xlab = "lag",
                       main = paste0("Lag = 0")))
    with(df.plot, lines(E, ll, lty = "dashed"))
    with(df.plot, lines(E, ul, lty = "dashed"))
    if(!missingArg(true.function)) with(df.plot, lines(E, true, col = "blue"))



    df.plot <- filter(out$est$surface, l == 5)
    with(df.plot, plot(E, mode, type = "l", ylim = c(min(c(ll)), max(c(ul))), ylab = "est", xlab = "lag",
                       main = paste0("Lag = 5")))
    with(df.plot, lines(E, ll, lty = "dashed"))
    with(df.plot, lines(E, ul, lty = "dashed"))
    if(!missingArg(true.function)) with(df.plot, lines(E, true, col = "blue"))


    if(!is.null(object$smooth$others)) {
      for(var. in unique(smooth.df$var)){
        smooth.df. <- filter(smooth.df, var == var.)
        with(smooth.df., plot(x, mode, type = "l", ylim = c(min(ll), max(ul)), ylab = "est", xlab = var.))
        with(smooth.df., lines(x, ll, lty = "dashed"))
        with(smooth.df., lines(x, ul, lty = "dashed"))
        if(!missingArg(true.function)){
          if(!is.null(true.function$smooth)) with(smooth.df., lines(x, true, col = "blue"))
        }
      }
    }

    devAskNewPage(oask)
  }



  invisible(out)
}
