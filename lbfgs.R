lbfgs <- function(test_fun, param) {

  if ( (is.null(param$autolog)) || (param$autolog==TRUE) ) {
    store_hist <<- TRUE # start history recording (specific to lbfgs)
    nbcalls <<- 0
    glob_xhist<<- matrix(,1,ncol=param$dim)
    glob_fhist<<- matrix(,1,1)
  }
  
  dim <- length(param$xinit)  
  if ( dim != param$dim) stop("dim and xinit do not have same dimension")
  LB <- param$LB
  UB <- param$UB
  # expand LB and UB to the number of dimensions
  if ((length(LB)==1) & (dim>1)) LB <- rep(LB[1],dim)
  if ((length(UB)==1) & (dim>1)) UB <- rep(UB[1],dim)
  # control parameters of the optimizer
  # trace : integer between 0 and 6 (see ?optim or source code)
  ctrl_param = list(trace=param$trace,maxit=param$maxit, lmm = param$lmm, factr=param$factr , pgtol=param$pgtol)
  
  optres <- optim(param$xinit, test_fun, method="L-BFGS-B",lower=LB, upper=UB, 
                 control=ctrl_param)
  
  res <- list(xhist=glob_xhist , fhist=glob_fhist , x_best=optres$par, f_best=optres$value)
  if ( (is.null(param$autolog)) || (param$autolog==TRUE) ) {
    store_hist <<- FALSE  # end history recording
  }

  return(res)
}
