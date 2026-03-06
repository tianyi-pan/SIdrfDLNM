#' Title Fit a single index distributed response function distributed lag non-linear model
#' (SID-DLNM)
#'
#' @param formula the formula for the effect of exposure at time t regarding the
#'   exposure history. Use sX() for the distributed lag term. For example y ~
#'   sX(t, x).
#' @param smooth the formula for the other smooth functions. The smooths are
#'   specified by mgcv::s. For example ~ s(z). Note also that the Gaussian
#'   random effect models are accommodated, by setting s(bs = "re").
#' @param unpen.smooth the formula for the other unpenalized smooth functions.
#' @param fe.varying the formula for time-varying covariates in linear
#'   components. For example ~ z
#' @param fe.cont the formula for time-invariant covariates in linear
#'   components. For example ~ z
#' @param dat the data frame containing the variables.
#' @param maxL the maximum lag, default 14.
#' @param kw the dimension for B-spline for weight function, default 20.
#' @param kE the dimension for B-spline for ACE response-exposure-function,
#'   default 20.
#' @param E.min the lower bound of ACE. The lower bound is set as the minimum of
#'   the exposure if E.min is missing, since this the range we are typically
#'   interested in.
#' @param E.max the upper bound of ACE. The upper bound is set as the maximum of
#'   the exposure if E.max is missing.
#' @param ifAIC whether or not to calculate the AIC, default TRUE.
#' @param pc the point constraint for ACE response-exposure-function, i.e. f(x0)
#'   = 0. Use sum-to-zero constraint if \code{pc = NULL}. Default \code{pc =
#'   NULL}.
#' @param par.start the starting values for BFGS in the outer optimization. If
#'   the argument is missing, the starting values are the fitted values from
#'   \code{mgcv::bam}.
#' @param par.fix the values of smoothing/dispersion parameters which is fixed
#'   in the outer optimization. For example, set \code{par.fix = c(10,NA,NA)}
#'   for fixing log(theta)=10 if the prior knowledge believes there is no
#'   over-dispersion.
#' @param upper.bound the upper bound of parameters in BFGS.
#' @param fit whether or not to fit the model. If \code{fit = FALSE}, the model 
#'   will be evalutated at the given parameters in \code{par.start}. 
#' @param CI the confidence level, default 0.95.
#' @param CI.R the number of sampling in the sampling method for CI, default
#'   1000.
#' @param CI.seed the random seed in the sampling method for CI, default 123.
#' @param eta whether or not to report the CI for the linear predictor eta =
#'   log(mu), default \code{FALSE}.
#' @param verbose whether or not to print the progress and diagnostic
#'   information, such as the gradients and function values in inner and middle
#'   optimizations.
#'
#' @return Object of class \code{sidDLNM_fit}. S3 function \code{summary},
#'   \code{residuals}, and \code{predict}.
#' @importFrom mgcv s
#' @importFrom data.table as.data.table
#' @importFrom data.table setorder
#' @export
SIdrfDLNM <- function(formula,
                    smooth = NULL,
                    unpen.smooth = NULL,
                    fe.varying = NULL,
                    fe.cont = NULL,
                    offset = NULL,
                    dat = NULL,
                    maxL = 14, kE = 20, kw = 10,
                    E.min,
                    E.max,
                    ifAIC = FALSE,
                    marginalAIC = FALSE,
                    ifNCV = FALSE,
                    kNCV = 0,
                    NCV.nthreads = 1,
                    conL = FALSE,
                    conLorder = 2,
                    pc = NULL,
                    par.start, par.fix, upper.bound,
                    fit = TRUE,
                    index.init.fix = TRUE,
                    con.index.init,
                    cindex,
                    CI = 0.95,
                    CI.R = 1000, CI.seed = 123,
                    eta = FALSE,
                    check.BFGS = FALSE,
                    maxit.LAML = 150,
                    tmp.file = NULL,
                    verbose = TRUE) {

  if(length(formula) != 3) stop("Incorrect formula. Please set a formula for example y~sX(t,x).")
  sXobject <- eval(formula[[3]])

  ## change character columns to factor. support random effects defined by mgcv::s(bs = "re")
  chr_col <- which(sapply(dat, class) == "character")
  if(length(chr_col) >= 1) {
    for(col. in chr_col) dat[,col.] <- factor(dat[,col.])
  }



  sXdat <- dat

  ## change colnames in dataframe as x, t and y
  x.names.list <- sXobject$x

  colnames(sXdat)[which(colnames(sXdat) == sXobject$t)] <- "t"
  colnames(sXdat)[which(colnames(sXdat) == as.character(formula[[2]]))] <- "y"


  ## standardize if not
  verbose.std <- FALSE
  for (xx in x.names.list) {
    if((mean(sXdat[[xx]]) > 0.1) | (abs(sd(sXdat[[xx]]) - 1) > 0.1) ) {
      sXdat[[xx]] <- (sXdat[[xx]] - mean(sXdat[[xx]])) / sd(sXdat[[xx]])
      verbose.std <- TRUE
    }
  }
  if(verbose & verbose.std) cat("Standardize the exposures. \n")

  ## byvar
  byvar <- sXobject$by


  maxLreal <- maxL+1

  # shift <- min(min(sXdat$x1,sXdat$x2,sXdat$x3), 0)
  shift <- 0
  # sXdat$x1 <- sXdat$x1 - shift
  # sXdat$x2 <- sXdat$x2 - shift
  # sXdat$x3 <- sXdat$x3 - shift


  formula.list <- list(formula = formula,
                       fe.cont = fe.cont,
                       fe.varying = fe.varying,
                       smooth = smooth)

  ### CONSTRUCTIONS #######
  ## time non-varying for group
  if(length(unique(sXdat$t)) < nrow(sXdat)) {

    if(is.null(fe.cont) & is.null(byvar)) stop("The exposure process data are duplicated at some time point. Please provide fe.cout.")

    if(!is.null(fe.cont)){
      # group_name.fe.cont <- fe.cont[[2]]
      # if(length(group_name.fe.cont) >= 2) {
      #   group_name.fe.cont <- as.character(group_name.fe.cont[2:length(group_name.fe.cont)])
      # } else {
      #   group_name.fe.cont <- as.character(group_name.fe.cont)
      # }
      group_name.fe.cont <- all.vars(fe.cont)
    } else {
      group_name.fe.cont <- NULL
    }
    if(!is.null(byvar)) {
      group_name.byvar <- byvar
    } else {
      group_name.byvar <- NULL
    }

    group_name <- c(group_name.fe.cont, group_name.byvar)

    ## split data
    sXdatlist <- split(as.data.table(sXdat), by = group_name) # use split in data.table to preserve order

  } else {
    group_name <- NULL
    sXdatlist <- list(sXdat)
  }

  ## set t starting at 1 and sort
  min_t <- min(sXdat$t)
  sXdatlist <- lapply(sXdatlist, function(sXdati) {
    sXdati$t <- sXdati$t - min_t + 1
    data.table::setorder(sXdati, t)
    return(data.frame(sXdati))
  })





  DLprepare <- lapply(sXdatlist, function(sXdati){
    DLi <- lapply(x.names.list, function(xx) {

      x <- sXdati[[xx]]
      y <- sXdati$y
      t <- sXdati$t
      Nti <- length(y) - maxL

      B_inner1 <- tsModel::Lag(x, 0:maxL)

      removed.t <- t[1:maxL]
      t <- t[-(1:maxL)] # delete the first maxL days
      y <- y[-(1:maxL)] # delete the first maxL days
      x <- x[-(1:maxL)] # delete the first maxL days
      B_inner1 <- B_inner1[-(1:maxL), ]


      return(list(B_inner = B_inner1,
                  removed.t = removed.t))
    })

    return(list(B_inner = lapply(DLi, "[[", "B_inner"),
                removed.t =  DLi[[1]][["removed.t"]]
    ))
  })

  # B_inner <- list(B_inner1, B_inner2, B_inner3)
  B_inner.list <- Reduce(function(x, y) Map(function(...){
    do.call(rbind, list(...))
  }, x, y), lapply(DLprepare, "[[", "B_inner"))


  M <- length(B_inner.list) # number of exposures
  if(missingArg(E.max))  E.max <- max(apply(dat[unlist(x.names.list)], 1, function(s) sqrt(sum(s^2))))
  if(missingArg(E.min))  E.min <- -1.0*E.max
  removed.t <- lapply(DLprepare, "[[", "removed.t")



  ## remove the starting time points
  sXdatlist <- mapply(function(sXdati, removed.ti) return(sXdati[-which(sXdati$t %in% removed.ti),]),
                      sXdatlist, removed.t, SIMPLIFY = FALSE)


  sXdat <- do.call("rbind", sXdatlist)


  # remove NA rows
  na.id <- is.na(sXdat$y)
  sXdat <- sXdat[!na.id, ]

  # B_inner <- list(B_inner1[!na.id, ], B_inner2[!na.id, ], B_inner3[!na.id, ])
  B_inner <- lapply(B_inner.list, function(xx) xx[!na.id, ])




  ### prepare NCV. construct a list for neighborhood of i
  sXdat$ii <- 1:nrow(sXdat)
  if(!is.null(group_name)) {
    NCVsXdat <- split(as.data.table(sXdat), by = group_name)
  } else {
    NCVsXdat <- list(as.data.table(sXdat))
  }

  nei.list <- lapply(NCVsXdat, function(NCVsXdati) {
    ni <- nrow(NCVsXdati)
    lapply(1:ni, function(ii)
      NCVsXdati$ii[seq(max(ii - kNCV, 1), min(ii + kNCV, ni))]
    )
  })
  nei.list <- do.call("c", nei.list)



  ## 1. time-varying fixed effects
  if(!is.null(fe.varying)) {
    Xlin <- stats::model.matrix(fe.varying,data=sXdat)[,-1,drop=F]
  } else {
    Xlin <- matrix(1, nrow = nrow(sXdat))
  }
  if(!is.null(unpen.smooth)) {
    unpen.smooth <- lapply(lapply(attr(terms(unpen.smooth),"term.labels"), function(text) parse(text = text)), eval)
    unpen.SS <- lapply(lapply(unpen.smooth,mgcv::smoothCon,data=sXdat,absorb.cons = TRUE),'[[',1) ## construct unpenalized smooth terms
    unpen.Xlist <- lapply(unpen.SS,'[[','X')
    unpen.X <- Reduce(cbind,unpen.Xlist) ## design matrix
    Xlin <- cbind(Xlin, unpen.X)
  }
  ### code following https://github.com/awstringer1/mam/blob/master/R/mam.R
  ## start following ##
  ## 2.2 smooth term
  numsmooth <- 0 # initialize
  if(!is.null(smooth)){
    smooth <- lapply(lapply(attr(terms(smooth),"term.labels"), function(text) parse(text = text)), eval)

    SS <- lapply(lapply(smooth,mgcv::smoothCon,data=sXdat,absorb.cons = TRUE),'[[',1) ## construct smooth terms
    numsmooth <- length(smooth) # Number of smooth terms

    EE <- lapply(lapply(lapply(SS,'[[','S'),'[[',1),eigen) ## eigen decomposition for penalty matrix

    p <- sapply(lapply(EE,'[[','vectors'),ncol) ## dimension of penalty matrix
    r <- sapply(lapply(EE,'[[','values'),function(x) sum(x>.Machine$double.eps)) ## rank of penalty matrix
    m <- p-r ## dim of null space (minus intercept)
    URlist <- mapply(function(x,y) x[ ,1:y],lapply(EE,'[[','vectors'),r,SIMPLIFY = FALSE)
    UFlist <- mapply(
      function(x,y,z) {
        if (y<z) return(x[ ,(1+y):z])
        newsparsemat(z,z)
      },lapply(EE,'[[','vectors'),r,p,SIMPLIFY = FALSE)
    URlist <- lapply(URlist,cbind) # Ensure they stay matrices
    UFlist <- lapply(UFlist,cbind) # Ensure they stay matrices

    UR <- Matrix::bdiag(URlist)
    UF <- Matrix::bdiag(UFlist)
    # if m=1 UF gets coerced to numeric
    if (!is.matrix(UF)) UF <- cbind(UF)

    Dpi <- Matrix::Diagonal(sum(r),1 / sqrt(Reduce(c,lapply(lapply(EE,'[[','values'),function(x) x[x>.Machine$double.eps]))))

    Xlist <- lapply(SS,'[[','X')
    X <- Reduce(cbind,Xlist) ## design matrix

    Xrand <- as.matrix(X %*% UR %*% Dpi) ## reparametrized
    Xfix <- as.matrix(X %*% UF)
    dups <- !base::duplicated(t(Xfix)) & apply(Xfix,2,function(x) !all(x==0)) # Remove the duplicated intercepts
    if (length(dups) > 1) Xfix <- Xfix[ ,which(dups)]

    model.choice = "with.smooth"
  } else{
    model.choice = "without.smooth"
  }


  # 2.3 add the intercept
  if(exists("Xfix")) {
    Xfix <- cbind(Xlin,Xfix) ## linear effect + unpenalized columns
  } else {
    Xfix <- Xlin
  }
  if (any(apply(Xfix,2,function(x) all(x==1)))) Xfix <- Xfix[,-(apply(Xfix,2,function(x) all(x==1)))]

  if (!is.null(fe.cont)){
    Xgroup <- stats::model.matrix(fe.cont,data=sXdat)[,-1,drop=FALSE]
    if (any(apply(Xgroup,2,function(x) all(x==1)))) Xgroup <- Xgroup[,-(apply(Xgroup,2,function(x) all(x==1)))] # remove the intercept
    Xfix <- cbind(1, Xgroup, Xfix) # add the intercept to the first column
  } else {
    Xfix <- cbind(1,Xfix) # add the intercept to the first column
  }


  ### END following https://github.com/awstringer1/mam/blob/master/R/mam.R #######

  ## 2.4 offset
  if(!is.null(offset)) {
    Xoffset.original <- as.vector(stats::model.matrix(offset,data=sXdat)[,-1,drop=F])
    Xoffset.original.min <- min(Xoffset.original)
    Xoffset <- Xoffset.original - Xoffset.original.min
  } else {
    Xoffset.original <- as.vector(rep(0, nrow(sXdat)))
    Xoffset.original.min <- 0
    Xoffset <- as.vector(rep(0, nrow(sXdat)))
  }



  ### Model Fitting


  N <- nrow(sXdat)
  y <- sXdat$y
  t <- sXdat$t


  if(missingArg(upper.bound)) upper.bound <- c(10, 20, 20)
  lower.bound <- c(-2,-8,-8)


  SSf <- mgcv::smoothCon(s(E, bs = "bs", k = kE), absorb.cons = FALSE,
                         data = data.frame(E = seq(E.min, E.max, length.out = N)))[[1]]

  if(is.null(pc)) {
    ## sum-to-zero constraint reparameterization for SSf
    SSfCon <- mgcv::smoothCon(s(E, bs = "bs", k = kE), absorb.cons = TRUE,
                         data = data.frame(E = seq(E.min, E.max, length.out = N)))[[1]]

    V <- t(SSf$X) %*% as.vector(rep(1,nrow(SSf$X)))

  } else {
    ## point constraint reparameterization for SSf
    # SSfCon <- mgcv::smoothCon(s(E, bs = "bs", k = kE, pc = pc), absorb.cons = TRUE,
    #                      data = data.frame(E = seq(E.min, E.max, length.out = N)))[[1]]
    ## NOTE: (s(E, bs = "bs", k = kE, pc = pc) will give a different reparametrization.
    # max(abs(SSf$X %*% Zf - SSfCon$X)) is not zero!
    # We build the design matrix and penalty matrix by hand

    V <- matrix(Bsplinevec2(pc, SSf$knots, 4), ncol = 1)

  }


  QRf <- qr(V)
  Qf <- qr.Q(QRf, complete = TRUE)
  Zf <- Qf[,2:ncol(Qf)]
  ## Check whether the design matrices are identical
  # max(abs(SSf$X %*% Zf - SSfCon$X))
  # max(abs(t(Zf) %*% SSf$S[[1]] %*% Zf - SSfCon$S[[1]]))


  if(is.null(pc)) {
    SfCon <- SSfCon$S[[1]]
  } else {
    ## NOTE: (s(E, bs = "bs", k = kE, pc = pc) will give a different reparametrization.
    # max(abs(SSf$X %*% Zf - SSfCon$X)) is not zero!
    # We build the design matrix and penalty matrix by hand
    SfCon <- t(Zf) %*% SSf$S[[1]] %*% Zf
  }


  SSw <- mgcv::smoothCon(s(l, bs = "bs", k = kw), absorb.cons = FALSE,
                         data = data.frame(l = seq(0, maxL, length.out = N)))[[1]]

  # SSw$knots
  Sf <- t(Zf) %*% SSf$S[[1]] %*% Zf %x% diag(1, kw) # beta_xi * vector(beta_l), ... , beta_x(kE-1) * vector(beta_l)
  Sw <- diag(1, kE-1) %x% SSw$S[[1]]

  alpha_f.init.default <- rep(0,(kE-1)*kw)

  betaF.init.default <- rep(0, ncol(Xfix))

  # delta.init.default <- rep(1, M - 1)


  ## full rank Sw and Sf

  EES <- eigen(Sw/norm(Sw, type = "F") + Sf/norm(Sf, type = "F"))
  pS <- ncol(EES$vector)
  rS <- sum(EES$values > .Machine$double.eps)
  mS <- pS - rS
  URS <- EES$vectors[,seq_len(rS)]
  # UFS <- EES$vectors[,rS + seq_len(mS)]
  SwR <- t(URS) %*% Sw %*% URS
  SfR <- t(URS) %*% Sf %*% URS




  ## prepare index weight
  Iindex <- diag(1, nrow = M, ncol = M) # the so-called basis function for index weight

  # cindex <- c(1,0,0) # for alpha_1 > 0
  # cindex <- as.vector(rep(1,M)) # for alpha_1 + alpha_2 + alpha_3 > 0
  if(missingArg(cindex)) cindex <- rep(1,M)
  Vindex <- Iindex %*% cindex
  Qfindex <- qr.Q(qr(Vindex), complete = TRUE)
  Zfindex <- Qfindex[,2:ncol(Qfindex)]
  Bindex <- cbind(cindex, Iindex %*% Zfindex)

  delta.init.default <- solve(Bindex, rep(1, M))
  delta.init.default <- delta.init.default / delta.init.default[1]
  delta.init.default <- delta.init.default[-1]
  ## the yielded index.init is c(1,1,1)

  if(!missingArg(con.index.init)) {
    if(length(con.index.init) != (M-1)) con.index.init <- delta.init.default
    delta.init.default <- con.index.init
    if(verbose) cat("con.index.init:", delta.init.default, "\n")
  }

  ## A list to store the optimization results at given smoothing log_theta.given, log_smoothing_f.given, log_smoothing_w.given.
  if(model.choice == "with.smooth") {
    betaR.init.default <- rep(0,ncol(Xrand))
    LAMLenv <- list(par = NULL,
                    fn = NULL,
                    gr = NULL,
                    mod = list(alpha_f.mod = alpha_f.init.default,
                              betaF.mod = betaF.init.default,
                              betaR.mod = betaR.init.default,
                              delta.mod = delta.init.default),
                    iter = 0,
                    PLg = 0)
  } else {
    # LAMLenv <- list(par = NULL,
    #             fn = NULL,
    #             gr = NULL,
    #             mod = list(phi.mod = phi.init.default, alpha_f.mod = alpha_f.init.default,
    #                         betaF.mod = betaF.init.default))
    cat("NOT supported. TODO")
    return(0)
  }

  initLAMLenv <- function() {
    LAMLenv <<- list(par = NULL,
                    fn = NULL,
                    gr = NULL,
                    mod = list(alpha_f.mod = alpha_f.init.default,
                               betaF.mod = betaF.init.default,
                               betaR.mod = betaR.init.default,
                               delta.mod = delta.init.default),
                    iter = 0,
                    PLg = 0)
  }

  ## upper bound
  if(model.choice == "with.smooth") {
    upper.bound <- c(upper.bound, rep(20, numsmooth))
    lower.bound <- c(lower.bound, rep(-8, numsmooth))
  }


  build <- function(par) {
    log_theta.given <- par[1]
    log_smoothing_f.given <- par[2]
    log_smoothing_w.given <- par[3]
    if(model.choice == "with.smooth") logsmoothing.given <- par[-(1:3)]

    ## Inner opt
    alpha_f.init <- LAMLenv$mod$alpha_f.mod
    if(any(is.nan(alpha_f.init))) alpha_f.init <- alpha_f.init.default
    betaF.init <- LAMLenv$mod$betaF.mod
    if(any(is.nan(betaF.init))) betaF.init <- betaF.init.default
    delta.init <- LAMLenv$mod$delta.mod
    if(any(is.nan(delta.init))) delta.init <- delta.init.default

    if(model.choice == "with.smooth"){
      betaR.init <- LAMLenv$mod$betaR.mod
      if(any(is.nan(betaR.init))) betaR.init <- betaR.init.default
      mod.build <- sidDLNMbuild(y, B_inner, SSf$knots, SSw$knots, Sw, SwR, Sf, SfR,
                                 Xrand, Xfix, Zf, Xoffset, r,
                                 Bindex,
                                 alpha_f.init,
                                 delta.init,
                                 log_theta.given, log_smoothing_f.given, log_smoothing_w.given,
                                 betaR.init, betaF.init, logsmoothing.given)
    }
    if(verbose) cat("Done.")
    mod.build
  }

  ## profile likelihood
  getLAML <-  function(par, ifgradient = TRUE) {
    if(verbose) cat("parameters:", par, "\n")
    log_theta.given <- par[1]
    log_smoothing_f.given <- par[2]
    log_smoothing_w.given <- par[3]
    if(model.choice == "with.smooth") logsmoothing.given <- par[-(1:3)]
    ## Inner opt
    if(is.null(LAMLenv$par) | !all(par == LAMLenv$par)) {
      ## Inner opt
      alpha_f.init <- LAMLenv$mod$alpha_f.mod
      if(any(is.nan(alpha_f.init))) alpha_f.init <- alpha_f.init.default
      betaF.init <- LAMLenv$mod$betaF.mod
      if(any(is.nan(betaF.init))) betaF.init <- betaF.init.default
      delta.init <- LAMLenv$mod$delta.mod
      if(any(is.nan(delta.init))) delta.init <- delta.init.default

      if(LAMLenv$iter > 1) {
        if(LAMLenv$PLg > 0.1) {
          if(verbose) cat("reset init \n")
          set.seed(LAMLenv$iter)
          alpha_f.init <- alpha_f.init.default + rnorm(n = length(alpha_f.init.default), mean = 0, sd = 1)
          betaF.init <- betaF.init.default + rnorm(n = length(betaF.init.default), mean = 0, sd = 1)
          delta.init <- delta.init.default + rnorm(n = length(delta.init.default), mean = 0, sd = 1)
        }
      }

      if(model.choice == "with.smooth"){
        betaR.init <- LAMLenv$mod$betaR.mod
        if(any(is.nan(betaR.init))) betaR.init <- betaR.init.default

        LAML.results <- sidDLNMopt(mod.address,
                                    ad.address,
                                    alpha_f.init,
                                    delta.init,
                                    log_theta.given, log_smoothing_f.given, log_smoothing_w.given,
                                    betaR.init, betaF.init, logsmoothing.given,
                                    ifgradient, 
                                    verbose)
        LAMLenv$par <<- par
        LAMLenv$mod$alpha_f.mod <<- LAML.results$alpha_f.mod
        LAMLenv$mod$delta.mod <<- LAML.results$con_index_par.mod
        LAMLenv$mod$betaF.mod <<- LAML.results$betaF.mod
        LAMLenv$mod$betaR.mod <<- LAML.results$betaR.mod
        LAMLenv$fn <<- LAML.results$LAML.fn
        LAMLenv$gr <<- LAML.results$LAML.gradient
        LAMLenv$iter <<- LAMLenv$iter + 1
        LAMLenv$PLg <<- LAML.results$PLg
      } else {
        cat("NOT supported. TODO")
      }

    }
  }


  LAML.fn <- function(par){
    getLAML(par)
    if(verbose) cat("LAML: ", LAMLenv$fn, "\n\n")
    if(!is.null(LAMLenv$fn)){
      if(any(is.nan(LAMLenv$fn))) {
        if(verbose) cat("reset starting values for inner... \n\n")
        initLAMLenv()
        getLAML(par)
        if(verbose) cat("LAML: ", LAMLenv$fn, "\n\n")
        if(any(is.nan(LAMLenv$fn))) {
          LAMLenv$fn <- 1e8
        }
      }
    }
    LAMLenv$fn
  }


  ## gradient of LAML
  LAML.gr <- function(par){
    getLAML(par)
    if(verbose) cat("Gradient: ", LAMLenv$gr, "\n\n")
    if(!is.null(LAMLenv$gr)) {
      if(any(is.nan(LAMLenv$gr))) {
        if(verbose) cat("reset starting values for inner... \n\n")
        initLAMLenv()
        getLAML(par)
        if(verbose) cat("Gradient: ", LAMLenv$gr, "\n\n")
      }
    }
    LAMLenv$gr
  }



  if(missingArg(par.fix)){
    par.fix <- rep(NA, 3)
    if(model.choice == "with.smooth") par.fix <- rep(NA, 3+numsmooth)
  }
  par.fix.id <- !sapply(par.fix, is.na)


  ### model fitting
  if(verbose) {
    cat("Start model fitting... \n")
  }
  start <- Sys.time()


  if(missingArg(par.start)) {
    par.start <- c(6, 6, 6)
    if(model.choice == "with.smooth") par.start <- c(par.start, rep(6, numsmooth))



    tryCatch(
    expr = {
        y.mean <- mean(sXdat$y, na.rm = TRUE)

        suppressWarnings(
          dDLNM.fit <- dDLNM(formula.list, x.names.list, kE, kw, sXdat, Bindex, M, maxL, verbose, delta.init.default)
        )
        mod.mgcv <- dDLNM.fit$mod.mgcv

        par.start[1] <- max(family(mod.mgcv)$getTheta(), log(2*y.mean)) # choose a larger theta if theta from bam() is smaller than 2*mean(y). Var = mu + mu^2/theta. i.e. init. var = 1.5*mu
        par.start[2] <- log(mod.mgcv$sp)[1]
        par.start[3] <- log(mod.mgcv$sp)[2]
        if(length(par.start) > 3) par.start[4:length(par.start)] <- log(mod.mgcv$sp)[-c(1,2)]
        if(par.start[1] > 7) par.start[1] <- 7 # do not set a huge log(theta)
        if(par.start[2] > 18) par.start[2] <- 18 # do not set a huge log(smoothing_f)
        if(par.start[3] > 15) par.start[3] <- 15 # do not set a huge log(smoothing_w)

        for (i in 1:length(par.start)) {
          if((par.start[i] > (upper.bound[i] - 1)) | (par.start[i] < lower.bound[i])) par.start[i] <- 7
        }
        if(verbose) cat("Use results from dlnm package and mgcv::bam as the initial guess for BFGS. \n ")
      },
      warning = function(w) {
        if(verbose) {
          cat("some warnings from bam() for setting the initial guess. Its influence on the model fitting is ignorable, but you can set par.start by hand. \n")
          warning(w)
        }
      },
      error = function(e) {
        cat("use default starting values for BFGS ... \n")
      }
    )
  }

  par.fn <- par.fix
  par.fn[!par.fix.id] <- par.start[!par.fix.id]
  address.list <- build(par.fn)
  mod.address <- address.list$address.eigen
  ad.address <- address.list$address.cppad

  hessian <- FALSE
  if(fit) {
    opt.LAML <- optim(par.start[!par.fix.id],
                    fn = function(par.){
                                        par.fn <- par.fix
                                        par.fn[!par.fix.id] <- par.
                                        return(LAML.fn(par.fn))
                                        }, # objective function
                    gr = function(par.){
                                        par.fn <- par.fix
                                        par.fn[!par.fix.id] <- par.
                                        return(LAML.gr(par.fn)[!par.fix.id])
                                        }, # gradient function
                    method = "L-BFGS-B",
                    lower = lower.bound[!par.fix.id],
                    upper = upper.bound[!par.fix.id],
                    control = list(trace = verbose,
                                  maxit = maxit.LAML),
                    hessian = hessian
                    )

    if(all(opt.LAML$par == par.start[!par.fix.id])){
      ## get stuck in the starting values
      par.start[!par.fix.id] <- par.start[!par.fix.id]/2
      opt.LAML <- optim(par.start[!par.fix.id],
                        fn = function(par.){
                                            par.fn <- par.fix
                                            par.fn[!par.fix.id] <- par.
                                            return(LAML.fn(par.fn))
                                            }, # objective function
                        gr = function(par.){
                                            par.fn <- par.fix
                                            par.fn[!par.fix.id] <- par.
                                            return(LAML.gr(par.fn)[!par.fix.id])
                                            }, # gradient function
                        method = "L-BFGS-B",
                        lower = lower.bound[!par.fix.id],
                        upper = upper.bound[!par.fix.id],
                        control = list(trace = verbose,
                                      maxit = maxit.LAML),
                        hessian = hessian)
    }

  } else {
    opt.LAML <- NULL
    getLAML(par.start, ifgradient = FALSE)
    if(verbose) cat("Model fitted at par.start: ", par.start, "\n")
  }
  opttime <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  if(verbose){
    cat("Finished model fitting. It took ", round(opttime,5), " seconds.\n", sep = "")
  }


  ### obtain estimate ########
  ## point estimate
  log_theta.opt <- LAMLenv$par[1]
  log_smoothing_f.opt <- LAMLenv$par[2]
  log_smoothing_w.opt <- LAMLenv$par[3]
  alpha_f.opt <- LAMLenv$mod$alpha_f.mod
  betaF.opt <- LAMLenv$mod$betaF.mod
  delta.opt <- LAMLenv$mod$delta.mod

  index.opt <- c(Bindex %*% c(1, delta.opt) / sqrt(as.numeric(t(c(1, delta.opt)) %*% t(Bindex) %*% Bindex %*% c(1, delta.opt))))

  if(sum(index.opt * cindex) < 0.01) warning("The index parameter is close to the boundary cindex * alpha = 0
                                                it might cause an identifiability problem.
                                              Possible Solution: Change cindex and refit the model.
                                              Be cautious when setting cindex = c(1, 0, 0, ...), as this assumes the first exposure has a NON-NULL effect!")

  out = list(opt = opt.LAML,
             env = LAMLenv,
             CI.sample = NULL,
             point = list(log_theta = log_theta.opt,
                          log_smoothing_f = log_smoothing_f.opt,
                          log_smoothing_w = log_smoothing_w.opt,
                          alpha_f = alpha_f.opt,
                          betaF = betaF.opt,
                          delta  = delta.opt,
                          index = index.opt),
             smooth = list(wl = SSw, fE = NULL))

  if(is.null(pc)) {
    out$smooth$fE <- SSfCon
  } else {
    ## NOTE: (s(E, bs = "bs", k = kE, pc = pc) will give a different reparametrization.
    # max(abs(SSf$X %*% Zf - SSfCon$X)) is not zero!
    # We build the design matrix and penalty matrix by hand
    out$smooth$fE <- SSf
  }

  if(model.choice == "with.smooth") {
    log_smoothing.opt <- LAMLenv$par[-(1:3)]
    betaR.opt <- LAMLenv$mod$betaR.mod
    out$point$betaR = betaR.opt
    out$point$log_smoothing = log_smoothing.opt
  }

  ## CI
  if(verbose) {
    cat("Start sampling for CI ... \n")
    start <- Sys.time()
  }
  if(model.choice == "with.smooth") {
    sampled <- sidDLNMCI(mod.address,
                         CI.R, CI.seed,
                         eta, verbose)
  } else {
    # sampled <- sidDLNMCI_nosmooth(y, B_inner, SSf$knots, SwI, SfI, Dw,
    #                               Xfix, Zfnew, Xoffset,
    #                                  LAMLenv$mod$alpha_f.mod,
    #                                   phi.opt,
    #                                   log_theta.opt, log_smoothing_f.opt, log_smoothing_w.opt,
    #                                   betaF.opt,
    #                                   CI.R, CI.seed,
    #                                   eta, delta.method, verbose)
    cat("NOT supported. TODO")
  }

  eta.est <- sampled$eta_point
  eta_E.est <- sampled$eta_E_point
  eta_other.est <- sampled$eta_other_point

  out$eta_E = data.frame(est = eta_E.est)
  out$eta_other = data.frame(est = eta_other.est)
  out$eta = data.frame(est = eta.est)


  if(eta) {
    eta.CI <- apply(sampled$eta_sample_mat, 2, function(row.) quantile(row., c((1-CI)/2, 1-(1-CI)/2)))
    out$eta = data.frame(est = eta.est,
                         ll = eta.CI[1,],
                         ul = eta.CI[2,])
    eta_E.CI <- apply(sampled$eta_E_sample_mat, 2, function(row.) quantile(row., c((1-CI)/2, 1-(1-CI)/2)))
    eta_other.CI <- apply(sampled$eta_other_sample_mat, 2, function(row.) quantile(row., c((1-CI)/2, 1-(1-CI)/2)))
    out$eta_E = data.frame(est = eta_E.est,
                           ll = eta_E.CI[1,],
                           ul = eta_E.CI[2,])
    out$eta_other = data.frame(est = eta_other.est,
                               ll = eta_other.CI[1,],
                               ul = eta_other.CI[2,])
  }
  if(verbose){
    runningtime <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    cat("Finished sampling for CI. It took ", round(runningtime,5), " seconds.\n", sep = "")
  }


  out$CI.sample <- sampled
  out$CI <- CI

  ## index
  index.CI <- apply(out$CI.sample$index_sample, 2, function(col.) quantile(col., c(0.025, 0.975)))
  out$index <- data.frame(est = index.opt,
                          ll = index.CI[1, ,drop = TRUE],
                          ul = index.CI[2, ,drop = TRUE])
  ## return formula
  out$formula <- formula.list


  ## return data
  out$data <- list(maxL = maxL,
                   B_inner = B_inner,
                   E.max = E.max,
                   E.min = E.min,
                   kE = kE,
                   kw = kw,
                   y = y,
                   x = sXdat[unlist(x.names.list)],
                   t = t,
                   pc = pc # point constraint
                   )
  out$data$Sw = Sw
  out$data$Sf = Sf
  out$data$SwR = SwR
  out$data$SfR = SfR
  out$data$Xfix = Xfix
  out$data$knots_f = SSf$knots
  out$data$knots_w = SSw$knots
  out$data$Zf = Zf
  out$data$shift = shift
  out$data$offset = list(Xoffset.original = Xoffset.original,
                         Xoffset.original.min = Xoffset.original.min,
                         Xoffset = Xoffset)


  if(model.choice == "with.smooth") {
    out$data$Xrand = Xrand
    out$data$UR = UR
    out$data$UF = UF
    out$data$Dpi = Dpi
    out$data$r = r
    out$data$m = m
    out$smooth$others = SS
  }
  # if(!is.null(unpen.smooth)) {
  #   out$smooth$unpen.smooth = unpen.SS
  # }

  out$inputdata <- dat
  out$modeldata <- sXdat
  out$opttime <- opttime
  out$penalty <- "second"


  ## check convergence

  out$Hessian_inner <- sampled$Hessian_inner
  out$eigval_Hessian_inner <- eigen(sampled$Hessian_inner)$values


  if(marginalAIC) {

    if(verbose) cat("start calculating Hessian of LAML \n")
    par.start[!par.fix.id] <- opt.LAML$par
    H.LAML <- optimHess(par.start[!par.fix.id],
                        fn = function(par.){
                          par.fn <- par.fix
                          par.fn[!par.fix.id] <- par.
                          return(LAML.fn(par.fn))
                        }, # objective function
                        gr = function(par.){
                          par.fn <- par.fix
                          par.fn[!par.fix.id] <- par.
                          return(LAML.gr(par.fn)[!par.fix.id])
                        }, # gradient function
                        control = list(trace = verbose)
    )
    if(verbose) cat("finished! \n")

  } else {
    H.LAML <- diag(1, length(par.start))
    if(verbose) cat("Hessian of LAML is not calculated! Set marginalAIC = TRUE if it is needed.")
  }

  ## check convergence of BFGS
  # H.LAML <- opt.LAML$hessian
  out$H.LAML <- H.LAML


  if(check.BFGS) {
      e.H <- eigen(H.LAML)

      out$eigen.hessian <- e.H

      evals <- e.H$values

      if(abs(prod(evals)) < 1e-3) evals <- evals + sqrt(sum((out$env$gr)^2))

      suggest.step <- e.H$vectors %*% diag(1/abs(evals)) %*% t(e.H$vectors) %*% out$env$gr

      out$suggest.step <- suggest.step
      if( sqrt(sum((out$env$gr)^2)) > 0.2) cat("BFGS might not converge. You could try other par.start and rerun the model.")

    if(min(out$eigval_Hessian_inner) < 0.0001) {
      warning("The optimization algorithm might not converge. Try rerunning the model with par.start = ", c(out$opt$par - out$suggest.step), "\n")
    }

  } else {
    if(min(out$eigval_Hessian_inner) < 0.0001) {
      warning("The optimization algorithm might not converge. Try rerunning the model with par.start = ", c(out$opt$par), "\n")
    }
  }

  if(!is.null(tmp.file)) {
    out.tmp <- structure(out, class = "sidDLNM_fit")
    save(out.tmp, file = tmp.file)
    cat("Temporary file saved to ", tmp.file, "\n")
  }
  ## AIC
  if(ifAIC) {
    if(verbose) {
      cat("Calculate AIC ... \n")
      start <- Sys.time()
    }
    H.LAML.par <- H.LAML[-1,-1]

    # e.H.par <- eigen(H.LAML.par)
    # evals.par <- e.H.par$values
    # evals.par <- abs(evals.par)
    # cprior <- 1e-3 # following mgcv
    # Vrho <- Matrix::solve(e.H.par$vectors %*% diag(evals.par + cprior) %*% t(e.H.par$vectors))
    cprior <- 1e-3 # following mgcv
    Vrho <- Matrix::solve(H.LAML.par + diag(rep(cprior, nrow(H.LAML.par))))
    AIC <- ConditionalAICsidDLNM(mod.address, ad.address, Vrho, verbose, marginalAIC)
    out$AIC <- AIC

    if(verbose){
      runningtime <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      cat("Finished calculating AIC. It took ", round(runningtime,5), " seconds.\n", sep = "")
    }
  }

  if(!is.null(tmp.file)) {
    out.tmp <- structure(out, class = "sidDLNM_fit")
    save(out.tmp, file = tmp.file)
    cat("Temporary file saved to ", tmp.file, "\n")
  }

  if(ifNCV) {
    if(verbose) {
      cat("Compute NCV ... \n")
      start <- Sys.time()
    }
    NCVresults <- NCVsidDLNM(mod.address, nei.list, verbose, NCV.nthreads)
    out$NCVresults <- NCVresults

    if(verbose){
      runningtime <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      cat("Finished computing NCV. It took ", round(runningtime,5), " seconds.\n", sep = "")
    }
  }

  if(!is.null(tmp.file)) {
    out.tmp <- structure(out, class = "sidDLNM_fit")
    save(out.tmp, file = tmp.file)
    cat("Temporary file saved to ", tmp.file, "\n")
  }

  structure(out, class = "sidDLNM_fit") # S3 class
}
