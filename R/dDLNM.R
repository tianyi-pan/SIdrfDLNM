#' Title fit a DLNM with fixed index parameters using dlnm package
#'
#' @param formula.list
#' @param kE
#' @param sXdat
#' @param M
#' @param verbose
#'
#' @returns
#' @importFrom tsModel Lag
#' @import mgcv
#' @import dlnm
#' @export
dDLNM <- function(formula.list, x.names.list, kE, kw, sXdat, Bindex, M, maxL, verbose, con.index.init) {

  if(verbose) cat("Fitting a DRF-DLNM using dlnm package and mgcv::bam() to obtain init ... \n")
  alpha <- c(1,con.index.init)
  # kk <- sqrt(sum(alpha^2))
  # alpha <- alpha/kk  ## so now ||alpha||=1
  kk <- as.numeric(t(alpha) %*% t(Bindex) %*% Bindex %*% alpha)
  alpha <- Bindex %*% alpha/kk  ## so now ||alpha||=1
  a <- as.vector(as.matrix(sXdat[unlist(x.names.list)]) %*% alpha)

  sXdat$a <- a
  Q <- Lag(a, 0:maxL)
  ## dlnm package only support "ps" and "cr".
  ## ps uses difference penalty
  cb <- crossbasis(Q[,1],lag=c(0,maxL),argvar=list(fun='ps',df=kE-1),
                   arglag=list(fun='ps', df=kw))

  cbPengam <- cbPen(cb)


  formula.list.mgcv <- Filter(Negate(is.null), formula.list[names(formula.list) %in% c("smooth","fe.varying")])


  formula.mgcv.pen <- "y~cb"
  if(length(formula.list.mgcv) > 0) {
    formula.other.mgcv <- paste(formula.list.mgcv, collapse = "+")
    formula.other.mgcv <- gsub("~", "", formula.other.mgcv)
    formula.other.mgcv <- gsub("\"", "'", formula.other.mgcv)
    formula.mgcv.pen <- paste0(formula.mgcv.pen, "+", formula.other.mgcv)
  }
  formula.mgcv.pen <- as.formula(formula.mgcv.pen)



  b.opt <- bam(formula.mgcv.pen,
               family=nb(),
               data = sXdat,
               paraPen=list(cb=cbPengam),
               method = "REML")
              
              

  if(verbose) cat("DRF-DLNM: Done! \n")

  return(list(mod.mgcv = b.opt)
  )
}





