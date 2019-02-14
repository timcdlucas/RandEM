#' nonparametric REM (NPREM) and Parametric REM (PREM)
#' @description Estimates the density of animals from camera trap data. REM is a nonparamteric method but could be derived by assuming the counts follow a Poisson or Negative distribtuion
#' @param V The speed data, including zero observations
#' @param Y The count data
#' @param n Sample size of count data
#' @param m sample size of speed data
#' @param g average group size for species moving in groups. Species moving individually have average group size, g=1
#' @param kappa auxiliary parameter from negative binomial distribution
#' @param B Samples in the bootstrap
#' @param theta Angle of detection of camera trap
#' @param r Distance of detection of camera trap
#' @param t Time period for which the camera was functioning (set to 1 day)
#' @param A The total area an animal can cover for contact to occur in a sector-shaped detection zone
#' @param species The species of interest
#' @param CountData The count data object
#' @param SpeedData The speed data object


#NONPARAMETRIC REM++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if(g == 1){
             EstimatedDensityNPREM =  (mean(Y)/A/ mean(V))
         }else{
             EstimatedDensityNPREM =  (mean(Y)/A/ mean(V))*g
         }
  EstimatedDensityNPREM

  #Estimated variance of density using the Delta method or Taylor series approximation
  StandardErrorNPREM = sqrt(EstimatedDensityNPREM^2*(var(Y)/n/mean(Y)^2 + var(V)/m/mean(V)^2))
  
  #Estimated standard error on the log scale
  LogStandardErrorNPREM = StandardErrorNPREM/EstimatedDensityNPREM
  
   

#BOOTSTRAP NONPARAMETRIC REM++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     NPREMBootstrap = function(Y, V, B) {
                                n = length(Y)
                                m = length(V)
                                EstimatedDensity = numeric(B)
 
   #sample the counts and speeds B times with replacement and estimate the density-----------------------------                                                          
                           for (nsim in 1: B) {
                                 YY = sample(Y, n, replace=TRUE)
                                 VV = sample(V, m, replace=TRUE)
                              
                                      if(g==1){
                                          EstimatedDensity[nsim] = mean(YY) / A / mean(VV)
                                      }else{
                                          EstimatedDensity[nsim] = (mean(YY) / A / mean(VV))*g
                                        }
                                    }
                                          EstimatedDensity
                      }

BootstrapEstimatedDensityNPREM  = NPREMBootstrap(Y, V, B)
BootstrapStandardErrorNPREM     = sqrt(var(BootstrapEstimatedDensityNPREM))
LogBootstrapStandardErrorNPREM  = sd(log(BootstrapEstimatedDensityNPREM))

           #Percentile confidence intervals---------------------------
              qqPercentile   = quantile(BootstrapEstimatedDensityNPREM ,c(0.025,0.975))   
              PercentileQ025 =  qqPercentile[1]
              PercentileQ975 =  qqPercentile[2]
              
          #Bias-corrected confidence intervals------------------------ 
              ObservedDensity = EstimatedDensityNPREM 
              b         = qnorm((sum(BootstrapEstimatedDensityNPREM  > ObservedDensity)+sum(BootstrapEstimatedDensityNPREM ==ObservedDensity)/2)/length(BootstrapEstimatedDensityNPREM))
              alpha     = 0.05                              
              z         = qnorm(c(alpha/2,1-alpha/2))      
              p         = pnorm(z-2*b)                      
              qqBiasCorrected   = quantile(BootstrapEstimatedDensityNPREM ,p=p)                
              BiasCorrectedQ025 = qqBiasCorrected[1]
              BiasCorrectedQ925 = qqBiasCorrected[2]

          #Fieller's confidence intervals------------------------------
              Fieller       = ttestratio(Y, A*V)
              Fieller.low   = Fieller$conf.int[1]
              Fieller.upper = Fieller$conf.int[2]
              
              
    Estimates = cbind(EstimatedDensityNPREM, mean(BootstrapEstimatedDensityNPREM), median(BootstrapEstimatedDensityNPREM))
    colnames(Estimates) = c("EstimatedDensityNPREM", "MeanBootstrapEstimatedDensityNPREM", "MedianBootstrapEstimatedDensityNPREM")
              
    StandardErrorEstimates  = cbind(StandardErrorNPREM, BootstrapStandardErrorNPREM, LogStandardErrorNPREM, LogBootstrapStandardErrorNPREM)
              



#PARAMETRIC REM++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
  PREMPoisson = function(param,Y,V){
                    EstimatedDensity = exp(param)
                    lambda           = A*mean(V)*EstimatedDensity
                    -sum(dpois(Y, lambda,log=TRUE))}  
  
  PREMNegativeBinomial = function(param,Y,V){
                            EstimatedDensity = exp(param[1])
                            Kappa            = exp(param[2])
                            lambda           = A*mean(V)*EstimatedDensity
                            -sum(dnbinom(Y, mu=lambda, size = 1/Kappa, log=TRUE))}  
  
  # Fit NB model or Poisson if var(y) < mean(y)
  if (var(Y) < mean(Y)) {
    optimization        = optim(0, PREMPoisson,, Y,V, hessian=TRUE, method = "Brent",lower=-2, upper=1000,
                           control=list(reltol=1e-15,maxit=25000))
                           kappa.hat = 0
  }else{
    optimization       = optim(c(0,0), PREMNegativeBinomial,, Y,V, hessian=TRUE, method ="Nelder-Mead",
                          control=list(reltol=1e-15,maxit=25000))
                          kappa.hat  = exp(optimization$par[2])
  }
  
  
 EstimatedDensityPREM  = exp(optimization$par[1])
 StandardErrorPREM     = sqrt(EstimatedDensityPREM^2 * solve(optimization$hessian)[1,1])
  
