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
         maximization      = optim(c(0,0,0,0), PoissonZAGAiREM,, Y, V1, A, hessian=TRUE,
                             control=list(reltol=1e-15,maxit=25000))
                             Kappa = 0
                        
         maximizationGamma  = optim(c(0,0,0), PoissonGammaiREM,, Y, V1, A, hessian=TRUE,
                             control=list(reltol=1e-15,maxit=25000))
                             KappaGamma = 0               
     }else{
         maximization     = optim(c(0,0,0,0,0), NBZAGAiREM, , Y, V1, A, hessian=TRUE,
                             control=list(reltol=1e-15,maxit=25000))
                             Kappa  = exp(maximization$par[5])
                        
         maximizationGamma   = optim(c(0,0,0,0), NBGAMMAiREM, , Y, V1, A, hessian=TRUE,
                                               control=list(reltol=1e-15,maxit=25000))
                            KappaGamma  = exp(maximizationGamma$par[4])                
     }


#ZAGA Estimates ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ExpectedSpeed            = exp(maximization$par[1])
        ShapeParameter           = exp(maximization$par[2])
        EstimatedDensityiREM     = exp(maximization$par[3])
        WeightOfZeros            = 1/(1 + exp(-maximization$par[4]))
        Mu                       = ExpectedSpeed/(1 - WeightOfZeros)

        #Estimated standard errors from Hessian-------------------------------
        HessianExpectedSpeed         =  sqrt(ExpectedSpeed^2 *diag(MASS::ginv(maximization$hessian))[1])
        HessianShapeParameter        =  sqrt(ShapeParameter^2 *diag(MASS::ginv(maximization$hessian))[2])
        HessianEstimatedDensityiREM  =  sqrt(EstimatedDensityiREM^2 *diag(MASS::ginv(maximization$hessian))[3])
        HessianWeightOfZeros         =  sqrt(WeightOfZeros^2 *diag(MASS::ginv(maximization$hessian))[4])
        HessianKappa                 =  sqrt(Kappa^2 *diag(MASS::ginv(maximization$hessian))[5])
        HessianMu                    =  sqrt(c((1/(1-WeightOfZeros))^2 *diag(MASS::ginv(maximization$hessian))[1]  + 
                                        2*(((1/(1-WeightOfZeros))*(-ExpectedSpeed/(1-WeightOfZeros)^2))*MASS::ginv(maximization$hessian)[1,4])+
                                        ((-ExpectedSpeed/(1-WeightOfZeros)^2)^2*diag(MASS::ginv(maximization$hessian))[4])))
         
        #Approximate standard error using theoretical standard error evaluated at estimated density from iREM--------------------      
        varExpectedSpeed           = ((1-WeightOfZeros)* Mu^2*(WeightOfZeros + ShapeParameter^2 ))/m
        EstimatedLambda            = A * EstimatedDensityiREM * ExpectedSpeed
        TheoreticalError           = sqrt(EstimatedDensityiREM^2 * ((1 + Kappa * EstimatedLambda)/(n*EstimatedLambda) + varExpectedSpeed/(ExpectedSpeed^2)))
        
        #Output ----------------------------------------------------------------        
        AllEstimates               = c(ExpectedSpeed, ShapeParameter, EstimatedDensityiREM, WeightOfZeros, Kappa, Mu)
        AllErrors                  = c(HessianExpectedSpeed, HessianShapeParameter, HessianEstimatedDensityiREM, HessianWeightOfZeros,
                                       HessianKappa, HessianMu)
        
        HessianAndTheoretical      = c(HessianEstimatedDensityiREM, TheoreticalError) 
        
     return(
       list(parameters = AllEstimates, standarderrors = AllErrors, gamma = HessianAndTheoretical)
     )
        
}

        
        
        
        
        
        
#GAMMA Estimates++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        ExpectedSpeedGamma            = exp(maximizationGamma$par[1])
        ShapeParameterGamma           = exp(maximizationGamma$par[2])
        EstimatedDensityiREMGamma     = exp(maximizationGamma$par[3])
        ScaleParameterGamma           = ExpectedSpeedGamma/ShapeParameterGamma
        
        #Estimated standard errors from Hessian-------------------------------
        HessianExpectedSpeedGamma         =  sqrt(ExpectedSpeedGamma^2 *diag(MASS::ginv(maximizationGamma$hessian))[1])
        HessianShapeParameterGamma        =  sqrt(ShapeParameterGamma^2 *diag(MASS::ginv(maximizationGamma$hessian))[2])
        HessianEstimatedDensityiREMGamma  =  sqrt(EstimatedDensityiREMGamma^2 *diag(MASS::ginv(maximizationGamma$hessian))[3])
        HessianKappaGamma                 =  sqrt(KappaGamma^2 *diag(MASS::ginv(maximizationGamma$hessian))[4])
       
        #Approximate standard error using theoretical standard error evaluated at estimated density from iREM--------------------      
        varExpectedSpeedGamma           = ExpectedSpeedGamma/ShapeParameterGamma/m1
        EstimatedLambdaGamma            = A * EstimatedDensityiREMGamma * ExpectedSpeedGamma
        TheoreticalErrorGamma           = sqrt(EstimatedDensityiREMGamma^2 * ((1 + KappaGamma * EstimatedLambdaGamma)/(n*EstimatedLambdaGamma) + varExpectedSpeedGamma/(ExpectedSpeedGamma^2)))
        
        #Output ----------------------------------------------------------------        
        AllEstimatesGamma               = c(ExpectedSpeedGamma, ShapeParameterGamma, EstimatedDensityiREMGamma, KappaGamma)
        AllErrorsGamma                  = c(HessianExpectedSpeedGamma, HessianShapeParameterGamma, HessianEstimatedDensityiREMGamma,HessianKappaGamma)
        
        HessianAndTheoreticalGamma      = c(HessianEstimatedDensityiREMGamma, TheoreticalErrorGamma) 





#-------------------------------------------------------------------
#Poisson and ZAGA iREM 
PoissonZAGAiREM = function(param, Y, V, A){
#-------------------------------------------------------------------
                  ExpectedSpeed      = exp(param[1])
                  ShapeParameter     = exp(param[2])
                  Density            = exp(param[3])
                  WeightOfZeros      = 1/(1 + exp(-param[4]))
                  Mu                 = ExpectedSpeed/(1-WeightOfZeros)
                  lambda             = A*ExpectedSpeed*Density
                 -sum(dpois(Y, lambda, log=TRUE))-
                  sum(gamlss.dist::dZAGA(V, Mu, ShapeParameter, WeightOfZeros, log=TRUE))}


#-------------------------------------------------------------------
#Negative binomial and ZAGA iREM
NBZAGAiREM = function(param, Y, V, A){
#-------------------------------------------------------------------
                  ExpectedSpeed      = exp(param[1])
                  ShapeParameter     = exp(param[2])
                  Density            = exp(param[3])
                  WeightOfZeros      = 1/(1 + exp(-param[4]))
                  Kappa              = exp(param[5])
                  Mu                 = ExpectedSpeed/(1-WeightOfZeros)
                  lambda             = A*ExpectedSpeed*Density
                 -sum(dnbinom(Y, mu=lambda, size=1/Kappa, log=TRUE))-
                  sum(gamlss.dist::dZAGA(V, Mu, ShapeParameter, WeightOfZeros, log=TRUE))}




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



        
        
