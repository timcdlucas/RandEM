#' integrated Random Encounter Model (iREM)
#' @description Estimates the density of animals from camera trap data. iREM and NPREM will coincide when expected values equal the sample means of counts and speeds. Speeds are assumed to follow a ZAGA or gamma model
#' @param V1 The speed data that exclude zeros
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


#-------------------------------------------------------------------
#Poisson and gamma iREM 
PoissonGAMMAiREM = function(param, Y, V1){
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
NBGAMMAiREM = function(param, Y, V1){
#-------------------------------------------------------------------
        ExpectedSpeed      = exp(param[1])
        ShapeParameter     = exp(param[2])
        Density            = exp(param[3])
        Kappa              = exp(param[4])
        ScaleParameter     = ExpectedSpeed*ShapeParameter
        lambda             = A*ExpectedSpeed*Density
        -sum(dnbinom(Y, mu=lambda, size=1/Kappa, log=TRUE))-
        sum(dgamma(V1, ScaleParameter, ShapeParameter, log=TRUE))}



#--------------------------------------------------------------------
#optimization and output
#--------------------------------------------------------------------


# Fit NB model or Poisson if var(Y) < mean(Y)
   if (var(Y) < mean(Y)) {
  
       maximizationGamma  = optim(c(0,0,0), PoissonGammaiREM,, Y,V1, hessian=TRUE,
                            control=list(reltol=1e-15,maxit=25000))
                            KappaGamma = 0               
    }else{

      maximizationGamma   = optim(c(0,0,0,0), NBGAMMAiREM,, Y,V1, hessian=TRUE,
                            control=list(reltol=1e-15,maxit=25000))
                            KappaGamma  = exp(maximizationGamma$par[4])                
    }


     ExpectedSpeedGamma            = exp(maximizationGamma$par[1])
     ShapeParameterGamma           = exp(maximizationGamma$par[2])
     EstimatedDensityiREMGamma     = exp(maximizationGamma$par[3])
     ScaleParameterGamma           = ExpectedSpeedGamma/ShapeParameterGamma

     #Estimated standard errors from Hessian-------------------------------
     HessianExpectedSpeedGamma         =  sqrt(ExpectedSpeedGamma^2 *diag(ginv(maximizationGamma$hessian))[1])
     HessianShapeParameterGamma        =  sqrt(ShapeParameterGamma^2 *diag(ginv(maximizationGamma$hessian))[2])
     HessianEstimatedDensityiREMGamma  =  sqrt(EstimatedDensityiREMGamma^2 *diag(ginv(maximizationGamma$hessian))[3])
     HessianKappaGamma                 =  sqrt(KappaGamma^2 *diag(ginv(maximizationGamma$hessian))[4])

     #Approximate standard error using theoretical standard error evaluated at estimated density from iREM--------------------      
     varExpectedSpeedGamma           = ExpectedSpeedGamma/ShapeParameterGamma/m1
     EstimatedLambdaGamma            = A * EstimatedDensityiREMGamma * ExpectedSpeedGamma
     TheoreticalErrorGamma           = sqrt(EstimatedDensityiREMGamma^2 * ((1 + KappaGamma * EstimatedLambdaGamma)/(n*EstimatedLambdaGamma) + varExpectedSpeedGamma/(ExpectedSpeedGamma^2)))

     #Output ----------------------------------------------------------------        
     AllEstimatesGamma               = c(ExpectedSpeedGamma, ShapeParameterGamma, EstimatedDensityiREMGamma, KappaGamma)
     AllErrorsGamma                  = c(HessianExpectedSpeedGamma, HessianShapeParameterGamma, HessianEstimatedDensityiREMGamma,HessianKappaGamma)

     HessianAndTheoreticalGamma      = c(HessianEstimatedDensityiREMGamma, TheoreticalErrorGamma) 

