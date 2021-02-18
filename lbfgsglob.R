library(DiceDesign)

lbfgsglob <- function (test_fun, param){
  

  budget <- param$budget
  xinit <- param$xinit
  dim <- length(xinit)  

  current_par <- xinit
  current_value <- test_fun(xinit)
  best_par <- xinit
  best_value <- test_fun(xinit)

  PlanExpLHS <- lhsDesign(budget,dim)$design
  PlanOp<-maximinSA_LHS(PlanExpLHS)$design

  param_lbfgs <- param
  for (i in 1:budget){
      param_lbfgs$xinit <- PlanOp[i,]
      res_lbfgs <- lbfgs(ofwrapper, param_lbfgs)
      if (res_lbfgs$f_best< best_value){
        best_value<-res_lbfgs$f_best
        best_par <- res_lbfgs$x_best
        xhist<-res_lbfgs$xhist
        fhist<-res_lbfgs$fhist
      }
    }

res <- list() 
res$xhist <- rbind(xhist)
res$fhist <- rbind(fhist)
res$x_best <- best_value
res$f_best <- best_par
return(res)
}    