#' Title Obtain predictions from the fitted model by sidDLNM
#'
#' @param object object of class \code{sidDLNM_fit}.
#' @param dat the data frame containing the variables.
#' @param CI.seed random seed for sampling, default \code{123}.
#' @param CI confidence level, default \code{0.95}.
#' @param ...
#'
#' @return
#' @importFrom mgcv s
#' @importFrom data.table as.data.table
#' @importFrom data.table setorder
#' @export
predict.sidDLNM_fit <- function(object, dat, CI.seed = 123, CI = 0.95, ...) {

  if(missingArg(dat)) {
    if (!is.null(object$eta)) {
      return(object$eta)
    } else {
      dat <- object$inputdata
    }
  }

  smooth = object$formula$smooth
  formula = object$formula$formula
  fe.cont = object$formula$fe.cont
  fe.varying = object$formula$fe.varying
  kw = object$data$kw
  kE = object$data$kE
  CI = object$CI
  maxL = object$data$maxL
  maxLreal = maxL+1






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


  ## byvar
  byvar <- sXobject$by



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
  E.max <- object$data$E.max
  E.min <- object$data$E.min
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






  ## smooths from mgcv
  SSw <- object$smooth$wl

  SSfCon <- object$smooth$fE







  Xfix <- object$data$Xfix
  knots_f <- object$data$knots_f

  knots_w <- object$data$knots_w
  Zf <- object$data$Zf



  ## 1. time-varying fixed effects
  if(!is.null(fe.varying)) {
    Xlin <- stats::model.matrix(fe.varying,data=sXdat)[,-1,drop=F]
  } else {
    Xlin <- matrix(1, nrow = nrow(sXdat))
  }
    #   if(!is.null(unpen.smooth)) {
    #     unpen.smooth <- lapply(lapply(attr(terms(unpen.smooth),"term.labels"), function(text) parse(text = text)), eval)
    #     unpen.SS <- lapply(lapply(unpen.smooth,mgcv::smoothCon,data=sXdat,absorb.cons = TRUE),'[[',1) ## construct unpenalized smooth terms
    #     unpen.Xlist <- lapply(unpen.SS,'[[','X')
    #     unpen.X <- Reduce(cbind,unpen.Xlist) ## design matrix
    #     Xlin <- cbind(Xlin, unpen.X)
    #   }

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
  ## TODO: update offset
  # if(!is.null(offset)) {
  #   Xoffset.original <- as.vector(stats::model.matrix(offset,data=sXdat)[,-1,drop=F])
  #   Xoffset.original.min <- min(Xoffset.original)
  #   Xoffset <- Xoffset.original - Xoffset.original.min
  # } else {
  #   Xoffset.original <- as.vector(rep(0, nrow(sXdat)))
  #   Xoffset.original.min <- 0
  #   Xoffset <- as.vector(rep(0, nrow(sXdat)))
  # }
  Xoffset <- object$data$offset$Xoffset







  y <- sXdat$y



  Bindex <- object$data$Bindex
  if(is.null(Bindex)) {
    cindex <- rep(1, M)
    Iindex <- diag(1, nrow = M, ncol = M) # the so-called basis function for index weight
    Vindex <- Iindex %*% cindex
    Qfindex <- qr.Q(qr(Vindex), complete = TRUE)
    Zfindex <- Qfindex[,2:ncol(Qfindex)]
    Bindex <- cbind(cindex, Iindex %*% Zfindex)
  }




  mod.build <- sidDLNMbuild(y, B_inner, knots_f, knots_w, object$data$Sw, object$data$SwR, object$data$Sf, object$data$SfR,
                                 Xrand, Xfix, Zf, Xoffset, r,
                                 Bindex,
                                 object$point$alpha_f,
                                 object$point$delta,
                                 object$point$log_theta, object$point$log_smoothing_f, object$point$log_smoothing_w,
                                 object$point$betaR, object$point$betaF, object$point$log_smoothing)

  mod.address <- mod.build$address.eigen

  R.CI <- nrow(object$CI.sample[[1]])


  sampled <- sidDLNMCI(mod.address,
                        R.CI, CI.seed,
                        TRUE, FALSE,
                        object$Hessian_inner)
  eta.CI <- apply(sampled$eta_sample_mat, 2, function(row.) quantile(row., c((1-CI)/2, 1-(1-CI)/2)))


  eta.est <- sampled$eta_point
  eta_E.est <- sampled$eta_E_point
  eta_other.est <- sampled$eta_other_point

  out <- list(eta = data.frame(est = eta.est,
                               ll = eta.CI[1,],
                               ul = eta.CI[2,]),
              sample = sampled$eta_sample_mat)

  return(out)
}
