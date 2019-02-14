#' integrated Random Encounter Model (iREM)
#' @description Estimates the density of animals from camera trap data. iREM and NPREM will coincide when expected values equal the sample means of counts and speeds. Speeds are assumed to follow a ZAGA or gamma model
#' @param V1 The speed data that exclude zeros
#' @param Y The count data
#' @param n Sample size of count data
#' @param m sample size of speed data
#' @param A The total area an animal can cover for contact to occur in a sector-shaped detection zone
#'
#' @export
#' @example
#' out <- iREM(wwapspeed$V1, wwapdata$Y, length(wwapspeed$V1), length(wwapdata$Y), 10)

iREM <- function(V1, Y, m, n, A){

    # Fit NB model or Poisson if var(Y) < mean(Y)
      if (var(Y) < mean(Y)) {
    
        maximizationGamma  = optim(c(0,0,0), PoissonGammaiREM, , Y,V1, A, hessian=TRUE,
                                control=list(reltol=1e-15,maxit=25000))
        KappaGamma = 0               
      } else {

        maximizationGamma   = optim(c(0,0,0,0), NBGAMMAiREM,, Y,V1, A, hessian=TRUE,
                                control=list(reltol=1e-15,maxit=25000))
        KappaGamma  = exp(maximizationGamma$par[4])                
      }


     ExpectedSpeedGamma            = exp(maximizationGamma$par[1])
     ShapeParameterGamma           = exp(maximizationGamma$par[2])
     EstimatedDensityiREMGamma     = exp(maximizationGamma$par[3])
     ScaleParameterGamma           = ExpectedSpeedGamma / ShapeParameterGamma

     #Estimated standard errors from Hessian-------------------------------
     HessianExpectedSpeedGamma         =  sqrt(ExpectedSpeedGamma^2 *diag(MASS::ginv(maximizationGamma$hessian))[1])
     HessianShapeParameterGamma        =  sqrt(ShapeParameterGamma^2 *diag(MASS::ginv(maximizationGamma$hessian))[2])
     HessianEstimatedDensityiREMGamma  =  sqrt(EstimatedDensityiREMGamma^2 *diag(MASS::ginv(maximizationGamma$hessian))[3])
     HessianKappaGamma                 =  sqrt(KappaGamma^2 *diag(MASS::ginv(maximizationGamma$hessian))[4])

     #Approximate standard error using theoretical standard error evaluated at estimated density from iREM--------------------      
     varExpectedSpeedGamma           = ExpectedSpeedGamma/ShapeParameterGamma/m
     EstimatedLambdaGamma            = A * EstimatedDensityiREMGamma * ExpectedSpeedGamma
     TheoreticalErrorGamma           = sqrt(EstimatedDensityiREMGamma^2 * ((1 + KappaGamma * EstimatedLambdaGamma)/(n*EstimatedLambdaGamma) + varExpectedSpeedGamma/(ExpectedSpeedGamma^2)))

     #Output ----------------------------------------------------------------        
     AllEstimatesGamma               = c(ExpectedSpeedGamma, ShapeParameterGamma, EstimatedDensityiREMGamma, KappaGamma)
     AllErrorsGamma                  = c(HessianExpectedSpeedGamma, HessianShapeParameterGamma, HessianEstimatedDensityiREMGamma,HessianKappaGamma)

     HessianAndTheoreticalGamma      = c(HessianEstimatedDensityiREMGamma, TheoreticalErrorGamma) 

     return(
       list(parameters = AllEstimatesGamma, standarderrors = AllErrorsGamma, gamma = HessianAndTheoreticalGamma)
     )
        
}





#-------------------------------------------------------------------
#Poisson and gamma iREM 
PoissonGAMMAiREM = function(param, Y, V1, A){
#-------------------------------------------------------------------
         ExpectedSpeed      = exp(param[1])
         ShapeParameter     = exp(param[2])
         Density            = exp(param[3])
         ScaleParameter     = ExpectedSpeed*ShapeParameter
         lambda             = A*ExpectedSpeed*Density
         -sum(dpois(Y, lambda, log=TRUE))-
         sum(dgamma(V1, ScaleParameter, ShapeParameter, log=TRUE))}


#-------------------------------------------------------------------
#Negative binomial and gamma iREM
NBGAMMAiREM = function(param, Y, V1, A){
#-------------------------------------------------------------------
        ExpectedSpeed      = exp(param[1])
        ShapeParameter     = exp(param[2])
        Density            = exp(param[3])
        Kappa              = exp(param[4])
        ScaleParameter     = ExpectedSpeed*ShapeParameter
        lambda             = A*ExpectedSpeed*Density
        -sum(dnbinom(Y, mu=lambda, size=1/Kappa, log=TRUE))-
        sum(dgamma(V1, ScaleParameter, ShapeParameter, log=TRUE))}


