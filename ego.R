####     An EGO with decomposed steps for teaching
##  
#  2 acquisition functions are given: the mean of the GP or the EI (actually the EI is to be coded). 
#    internalFun points to the acquisition function.
#
#   Requirements (put them in you main file):
library(DiceDesign)
library(DiceKriging)
# source("./cmaes.R")


##### EGO function
egoSimple <- function(test_fun, param) {
  # recover parameters
  debugMode <- FALSE
  printCurve <- TRUE
  budget <- param$budget
  # xinit <- param$xinit  # one could provide the DoE
  zdim <- param$dim  
  LB <- param$LB
  UB <- param$UB
  # expand LB and UB to the number of dimensions
  if (length(LB)==1 & zdim>1) LB <- rep(LB[1],zdim)
  if (length(UB)==1 & zdim>1) UB <- rep(UB[1],zdim)
  
  x_hist <- matrix(, budget, zdim)
  y_hist <- matrix(, budget, 1)
  
  ##############
  cat("******* Start EGO \n")
  
  # create a LHS DoE
  ninit <- 10*zdim
  # Xinit <- data.frame(LB + (UB-LB)*lhsDesign(n = ninit,dimension = zdim)$design) # this version works only well all bounds are the same
  Xinit <- data.frame(matrix(LB,nrow=ninit,ncol=zdim,byrow=T) + matrix(UB-LB,nrow=ninit,ncol=zdim,byrow=T)*lhsDesign(n = ninit,dimension = zdim)$design)
  # calculate associated objective function
  Yinit <- apply(X = Xinit,MARGIN = 1,FUN = test_fun)
  imin <- which.min(Yinit)
  fmin <- Yinit[imin]
  xmin <- as.numeric(Xinit[imin,])
  x_hist[1:ninit,]<-as.matrix(Xinit)
  y_hist[1:ninit,1]<-Yinit
  # make a kriging model. 
  # I use a fairly high lower bound on the thetas here to prevent degenerated models
  if (!debugMode) {
    capture.output(GPmodel <-km(design = Xinit,response = Yinit,covtype="matern3_2",lower = rep(0.8,zdim),multistart = 20,nugget=1.e-8),file = nullfile())}
  else {
    GPmodel <-km(design = Xinit,response = Yinit,covtype="matern3_2",lower = rep(0.8,zdim),multistart = 20)    
  }
  
  cat("*** initial DoE and GP done\n")
  
  ##### GP mean as a function
  #   GPmodel and zdim passed as global variables, oh well ... 
  meanofGP <- function(x){
    x <- matrix(x,ncol=zdim)
    z <- predict(object = GPmodel,newdata=data.frame(x),type="UK")
    return(z$mean)
  }
  
  ##### -1*EI function
  #   GPmodel and fmin are passed as global variables, ouch!
  mEI <- function(x) {
    x <- matrix(x,ncol=zdim)
    y <- predict(object = GPmodel,newdata=data.frame(x),type="UK")
    EI<-(fmin-y$mean)*pnorm((fmin-y$mean)/y$sd)+y$sd*dnorm((fmin-y$mean)/y$sd)
    mEI<--EI
    mEI
  }
  
  # set the global internalFun variable to the meanofGP for minimizing the GP mean
  internalFun <- mEI
  
  
  for (i in 1:(budget-ninit)) {
    
    # optimize with CMA-ES the acquisition criterion to define the next iterate
    paramCMA <- list(LB=LB,UB = UB,budget = 3000, dim=zdim, xinit=runif(n = zdim,min = LB,max = UB),sigma=2.) # param for cmaes  
    optresCMA <- cmaes(internalFun, paramCMA)
    fnext <- test_fun(optresCMA$x_best) # the one call to the true objective function
    if (fnext<fmin){
      fmin <- fnext
      xmin <- optresCMA$x_best
    }   
    x_hist[i+ninit,]<-optresCMA$x_best
    y_hist[i+ninit]<-fnext
    # print out new iterate
    cat("*** iteration ",i,"\n")
    cat("    xnext=",optresCMA$x_best,"\n")
    cat("    fnext=",fnext," f_predicted :",meanofGP(optresCMA$x_best),"\n")
    cat("    best acquisition criterion=",optresCMA$f_best,"\n")
    # this would plot the convergence of the internal optimization
    # plot(optresCMA$ymeanhist,type="l", xlab="no. of iterations",ylab="f of xmean")
    
    # update DoE
    newX <- rbind(GPmodel@X,optresCMA$x_best)
    newY <- rbind(GPmodel@y,fnext)
    # build a new GP (which involves another internal optimization, typically a likelihood maximization)
    # the printout of the km function is thrown to the garbage for better lisibility, remove the capture.output for debug
    if (!debugMode) {
      capture.output(GPmodel <-km(design = newX,response = newY,covtype="matern3_2",lower = rep(0.5,zdim),multistart = 20,nugget=1.e-8),file = nullfile())}
    else {
      GPmodel <-km(design = newX,response = newY,covtype="matern3_2",lower = rep(0.5,zdim),multistart = 20,nugget=1.e-8)    
    }
    
  } # end EGO loops
  cat("******* EGO iterations completed\n")
  # a bit of reporting
  cat("fmin=",fmin,"\n")
  cat("xmin=",xmin,"\n")
  if (printCurve) {plot(x = seq(1,GPmodel@n),y = GPmodel@y,type="l",xlab="no. calls to f",ylab="f")}
  
  res <- list(xhist=x_hist, fhist=y_hist, x_best=xmin, f_best=fmin)
  return(res)
} # end function