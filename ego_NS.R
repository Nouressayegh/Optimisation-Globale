source("./normal_search.R")
library(DiceDesign)
library(DiceKriging)

ego_NS <- function(test_fun, param) {
  debugMode <- FALSE
  printCurve <- TRUE
  budget <- param$budget
  zdim <- param$dim  
  LB <- param$LB
  UB <- param$UB
  if (length(LB)==1 & zdim>1) LB <- rep(LB[1],zdim)
  if (length(UB)==1 & zdim>1) UB <- rep(UB[1],zdim)
  
  x_hist <- matrix(, budget, zdim)
  y_hist <- matrix(, budget, 1)
  cat("******* Start EGO \n")
  
  ninit <- 10*zdim
  Xinit <- data.frame(matrix(LB,nrow=ninit,ncol=zdim,byrow=T) + matrix(UB-LB,nrow=ninit,ncol=zdim,byrow=T)*lhsDesign(n = ninit,dimension = zdim)$design)
  Yinit <- apply(X = Xinit,MARGIN = 1,FUN = test_fun)
  imin <- which.min(Yinit)
  fmin <- Yinit[imin]
  xmin <- as.numeric(Xinit[imin,])
  x_hist[1:ninit,]<-as.matrix(Xinit)
  y_hist[1:ninit,1]<-Yinit
  if (!debugMode) {
    capture.output(GPmodel <-km(design = Xinit,response = Yinit,covtype="matern3_2",lower = rep(0.8,zdim),multistart = 20,nugget=1.e-8),file = nullfile())}
  else {
    GPmodel <-km(design = Xinit,response = Yinit,covtype="matern3_2",lower = rep(0.8,zdim),multistart = 20)    
  }
  
  cat("*** initial DoE and GP done\n")
  
  meanofGP <- function(x){
    x <- matrix(x,ncol=zdim)
    z <- predict(object = GPmodel,newdata=data.frame(x),type="UK")
    return(z$mean)
  }
  
  
  mEI <- function(x) {
    x <- matrix(x,ncol=zdim)
    y <- predict(object = GPmodel,newdata=data.frame(x),type="UK")
    EI<-(fmin-y$mean)*pnorm((fmin-y$mean)/y$sd)+y$sd*dnorm((fmin-y$mean)/y$sd)
    return(-EI)
  }
  internalFun <-mEI
  for (i in 1:(budget-ninit)) {
    
    # optimize with NS the acquisition criterion to define the next iterate
    
    paramNS <- list(LB=LB,UB = UB,budget = 150, dim=zdim, xinit=runif(n = zdim,min = LB,max = UB),sigma=1) 
    optresNS <- normal_search(internalFun, paramNS)
    ###############################################################
    while (is.null(optresNS$x_best)) {
      paramNS$xinit<-runif(n = dim,min = LB,max = UB)
      optresNS <- normal_search(internalFun, paramNS)
    }
    ###############################################################
    
    fnext <- test_fun(optresNS$x_best) 
    if (fnext<fmin){
      fmin <- fnext
      xmin <- optresNS$x_best
    }   
    x_hist[i+ninit,]<-optresNS$x_best
    y_hist[i+ninit]<-fnext
    cat("*** iteration ",i,"\n")
    cat("    xnext=",optresNS$x_best,"\n")
    cat("    fnext=",fnext," f_predicted :",meanofGP(optresNS$x_best),"\n")
    cat("    best acquisition criterion=",optresNS$f_best,"\n")
    
    newX <- rbind(GPmodel@X,optresNS$x_best)
    newY <- rbind(GPmodel@y,fnext)
    if (!debugMode) {
      capture.output(GPmodel <-km(design = newX,response = newY,covtype="matern3_2",lower = rep(0.5,zdim),multistart = 20,nugget=1.e-8),file = nullfile())}
    else {
      GPmodel <-km(design = newX,response = newY,covtype="matern3_2",lower = rep(0.5,zdim),multistart = 20,nugget=1.e-8)    
    }
    
  }
  cat("******* EGO iterations completed\n")
  cat("fmin=",fmin,"\n")
  cat("xmin=",xmin,"\n")
  if (printCurve) {plot(x = seq(1,GPmodel@n),y = GPmodel@y,type="l",xlab="no. calls to f",ylab="f")}
  
  res <- list(xhist=x_hist, fhist=y_hist, x_best=xmin, f_best=fmin)
  return(res)
} 