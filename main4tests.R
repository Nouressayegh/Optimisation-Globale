### Main code for repeated tests of optimizers
### but it can also be used for only 1 call to an optimizer
### Author: Rodolphe Le Riche
#
#   INPUTS:
#   Look into "I) Choice of the function", then 
#   "II) Choice of the optimizer and associated parameters "
#   below for hints about the main test campaign parameters to choose.
#
#   OUTPUTS:
#   a nameOfOptim_nameOfFunction_(additionalInfo).RData file 
#     for postprocessing (plotting) with postproc_tests_optimizers.R
#   a nameOfOptim_nameOfFunction_(additionalInfo).pdf file with some
#     some rough plots
#
#   Note : glob_... are global variables

T1<-Sys.time()
# code R

rm(list=ls())

# library("DiceOptim")
source("./test_functions.R")
source("./random_search.R")
source("./normal_search.R")
source("./of_wrapper.R")
source("./lbfgsglob.R")
source("./cmaes.R")
source("./NS_lbfgs.R")
source("./RS_lbfgs.R")
source("./egoSimple.R")

### I) Choice of the function  #############################
#
nameoffun <- "ackley" ## objective function, "quadratic", 
                     # "ackley", "michalewicz","sphere", "schwefel","rastrigin",tunnel", ...
zdim <-2
zLB <- -5 # lower bound on parameters
zUB <- 5 # upper bound on parameters
# don't change the global variables below unless you know what you are doing
glob_xstar <<- rep(2.5,zdim) # position of optimum on SOME of the test functions
glob_noisy <- FALSE  # make test functions noisy
glob_tau <- 1 # std dev of test functions' noise
# glob_estim_noise <- FALSE #this should go in KNF parameters

### II) Choice of optimizer and associated parameters ######
#
nameofoptim <- "egoSimple" # "normal_search" , "normal_search" , "lbfgs", "egoSimple" , "cmaes", "NS_lbfgs"
budget <- 30*zdim #1500 #70*zdim
no_test <<- 5
RandomInit <- FALSE # default
# making RandomInit TRUE creates different randomly chosen initial
# points during the tests. The variability is larger and in some cases
# one may wish to control the initial point, then make RandomInit FALSE
if (nameofoptim=="random_search") {
  param <- list(LB=zLB,UB = zUB,budget = budget,dim=zdim) # param for random_search 
} else if  (nameofoptim=="normal_search") {
  param <- list(LB=zLB,UB = zUB,budget = budget,dim=zdim, xinit=rep(-4,zdim),sigma=1) # param for normal_search
  RandomInit <- TRUE 
} else if  (nameofoptim=="lbfgsglob") {
  param <- list(LB=zLB,UB = zUB,maxit = floor(budget/zdim),dim=zdim, xinit=rep((zLB+zUB)/2,zdim), trace=0, lmm =1, factr=100,pgtol=1.e-2) # param for lbfgs
  RandomInit <- TRUE 
} else if  (nameofoptim=="cmaes") {
  param <- list(LB=zLB,UB = zUB,budget = budget, dim=zdim, xinit=rep(-4, zdim),sigma=1) # param for cmaes  
  RandomInit <- FALSE  
} else if  (nameofoptim=="egoSimple") {
  param <- list(LB=zLB,UB = zUB,budget = budget,dim=zdim) # param for EGO
} else if  (nameofoptim=="NS_lbfgs") {
  param <- list(LB=zLB,UB = zUB,budget = budget,dim=zdim, xinit=rep(-4,zdim),sigma=2.,maxit = 12,trace=0, lmm =1, factr=100,pgtol=1.e-2) # param for NS_lbfgs
  RandomInit <- FALSE 
} else if  (nameofoptim=="RS_lbfgs") {
  param <- list(LB=zLB,UB = zUB,budget = budget,dim=zdim, xinit=rep(-4,zdim),sigma=2.,maxit = 12,trace=0, lmm =1, factr=100,pgtol=1.e-2) # param for RS_lbfgs
  RandomInit <- FALSE  
} else {
  stop("unknown nameofoptim")
}

###############################################################
eval(parse(text=paste("fun <- ",nameoffun)))
eval(parse(text=paste("opt <- ",nameofoptim)))

store_hist <<- FALSE # TRUE for lbfgs, EGO (from DiceOptim, not the house), cmaes, maybe others, 
                     # WHEN USED AS MAIN optimizers (as opposed to internal optimizers), look inside files

cat("******* Start testing ",nameofoptim," on ",nameoffun,"\n")

hist <- list()
lastbudget <- vector()

for(i in 1:no_test){
  cat("   test no.",i,"\n")

  if  (RandomInit) {
    param$xinit=runif(n = zdim,min = zLB,max = zUB)
  }
  
  ##### the call to the optimizer
  optres <- opt(ofwrapper, param)
  ##### end of call to the optimizer

  lastbudget <- c(lastbudget,length(optres$fhist))
  hist[[i]] <- optres

}

cat("******* Done testing ",nameofoptim," on ",nameoffun,"\n")

###### Extract matrices and info for later plotting  ########

fbesthist <- list()
for (i in 1:no_test){
  fbesthist[[i]] <- hist[[i]]$fhist
}
fbesthist <- lapply(fbesthist,cummin)
y_lim_max_best <- max(unlist(lapply(fbesthist,max)))

y_lim_max <- -.Machine$double.xmax
y_lim_min <- .Machine$double.xmax
for(i in 1:no_test){
  themin <- min(hist[[i]]$fhist)
  themax <- max(hist[[i]]$fhist)
  if (themin < y_lim_min) { y_lim_min <- themin }
  if (themax > y_lim_max) { y_lim_max <- themax }
}

###### Plotting ########
aaa <- paste(nameofoptim, nameoffun,zdim, sep="_")
if(nameofoptim == "normal_search") {
  aaa <- paste(aaa, "sigma",param$sigma, sep="_")
  } 
aaapdf <- paste(aaa,"pdf", sep=".")

pdf(aaapdf,width=5, height=4 )

 plot(hist[[1]]$fhist,type="l",col=1, xlim = c(0, max(lastbudget)), 
      ylim=c(y_lim_min, y_lim_max), ylab="f",xlab="no. calls to f")
 if (no_test>1) {
   for (i in 2:no_test){
   lines(hist[[i]]$fhist,type="l",col=i)
   }
 }

 plot(fbesthist[[1]],type="l",col=1, xlim = c(0,max(lastbudget)), 
#      ylim=c(y_lim_min, y_lim_max_best), ylab="fbest",xlab="no. calls to f",log="y")
      ylim=c(y_lim_min, y_lim_max_best), ylab="fbest",xlab="no. calls to f")
 if (no_test>1) {
   for (i in 2:no_test){
   lines(fbesthist[[i]],type="l",col=i)
   }
 }
 
 if (!is.null(hist[[1]]$sigmahist)) {
   plot(hist[[1]]$sigmahist,type="l",col=1, 
        ylab="CMA-ES step size",xlab="no. of iterations",ylim=c(0,1.5*max(hist[[1]]$sigmahist)))
   for (i in 1:no_test){
     lines(hist[[i]]$sigmahist,col=i)
   }
 }

if (!is.null(hist[[1]]$eigenhist)) {
  for(i in 1:no_test) {
    matplot(hist[[i]]$eigenhist, type="l", lty=1, xlab="no. of iterations", ylab="eigenvalues")
    }
}
 
 if (length(hist[[1]]$noise_var_hist) > 0) {
   plot(hist[[1]]$noise_var_hist,type="l",col=1, 
        ylab="estimated noise var.",xlab="no. of iterations",ylim=c(0,1.5*max(hist[[1]]$noise_var_hist)))
   for (i in 1:no_test){
     lines(hist[[i]]$noise_var_hist,col=i)
   }
 }

 
xconvhist <- list()

if (nameofoptim != "cmaes") {
  if (!is.null(hist[[1]]$xhist)) {
    for (i in 1:no_test){
      A<-hist[[i]]$xhist-matrix(glob_xstar,nrow(hist[[i]]$xhist),zdim,byrow=T)
      xconvhist[[i]] <- sqrt(apply(A*A,1,sum))
      if (i == 1) {
        plot(xconvhist[[i]],type="l",col=1,ylab="||xhist- xstar||",xlab="no. of iterations",ylim=c(0,1.2*max(xconvhist[[i]])))
      } else {
        lines(xconvhist[[i]],type="l",col=i)
      }
    }
  }
} else {
  # for cmaes 
  if (!is.null(hist[[1]]$xmeanhist)) {
    for (i in 1:no_test){
      A<-hist[[i]]$xmeanhist-matrix(glob_xstar,nrow(hist[[i]]$xmeanhist),zdim,byrow=T)
      xconvhist[[i]] <- sqrt(apply(A*A,1,sum))
      if (i == 1) {
        plot(xconvhist[[i]],type="l",col=1,ylab="||xmean- xstar||",xlab="no. of iterations",ylim=c(0,1.2*max(xconvhist[[i]])))
      } else {
        lines(xconvhist[[i]],type="l",col=i)
      }
    }
  }
}
  
if(zdim==2) {
  n.grid <- 50
  x.grid <- seq(-5, 5, ,n.grid)
  X.grid <- expand.grid(x.grid, x.grid)
  y.grid <- apply(X.grid, 1, fun)
  for (i in 1:no_test){
    contour(x.grid, x.grid, matrix(y.grid, n.grid), nlevels=20,lwd=0.5, xlab="x1", ylab="x2")
    points(hist[[i]]$xhist[, 1], hist[[i]]$xhist[, 2], pch=4, col="blue")
  }
}

dev.off()

# save results in a file
ff <- paste(aaa, "RData", sep=".")
save(hist,fbesthist,xconvhist,file=ff)

