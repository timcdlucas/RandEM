

# calcDensity is the main function to calculate density.
# It takes parameters z, alpha, theta, r, animalSpeed, t
# count - The number of camera/acoustic counts or captures.
# alpha - Call width in radians.
# theta - Sensor width in radians.
# r - Sensor range in metres.
# animalSpeed - Average animal speed in metres per second.
# t - Length of survey in sensor seconds i.e. number of sensors x survey duration.
#
# calcAbundance calculates abundance rather than density and requires an extra parameter
# area - In metres squared. The size of the region being examined.


#' Internal function to calculate profile width as described in the text
#'
#'@inheritParams gremAbundance
#'@name calcProfileWidth

calcProfileWidth <- function(alpha, theta, r){
        if(alpha > 2*pi | alpha < 0) 
		stop('alpha is out of bounds. alpha should be in interval 0<a<2*pi')
        if(theta > 2*pi | theta < 0) 
		stop('theta is out of bounds. theta should be in interval 0<a<2*pi')

	if(alpha > pi){
	        if(alpha < 4*pi - 2*theta){
		        p <- r*(theta - cos(alpha/2) + 1)/pi
                } else if(alpha <= 3*pi - theta){
                        p <- r*(theta - cos(alpha/2) + cos(alpha/2 + theta))/pi
                } else {
                        p <- r*(theta + 2*sin(theta/2))/pi
                }
        } else {
        	if(alpha < 4*pi - 2*theta){
                        p <- r*(theta*sin(alpha/2) - cos(alpha/2) + 1)/pi
 		} else {
                        p <- r*(theta*sin(alpha/2) - cos(alpha/2) + cos(alpha/2 + theta))/pi
                }
        }
        return(p)
}


#' Calculate population density using the gREM.
#'
#' Calculate population density from count data using the gREM. 
#'
#' Note that the REM of Rowcliffe et al. 2008 can be used by setting alpha to 
#'    2*pi and the Gas Model of Yapp 1956 can be used setting alpha and theta
#'    to 2*pi.
#'
#'
#'@references Lucas, T. C. D., Moorcroft, E. A., Freeman, R., Rowcliffe, J. M.,
#'    Jones, K. E. (2015), A generalised random encounter model for estimating 
#'    animal density with remote sensor data. Methods in Ecology and Evolution. 
#'    doi: 10.1111/2041-210X.12346
#'
#' Rowcliffe, J., Field, J., Turvey, S. & Carbone, C. (2008) Estimating animal 
#'    density using camera traps without the need for individual recognition. 
#'    Journal of Applied Ecology, 45, 1228-1236.
#'
#' Yapp, W. (1956) The theory of line transects. Bird Study, 3, 93-104.
#'@inheritParams gremAbundance
#'@seealso  \code{\link{gremAbundance}}
#'@name gremDensity
#'@export

gremDensity <- function(count, alpha, theta, r, v, tm){
        # Check the parameters are suitable.
        if(count <= 0 | !is.numeric(count)) stop('Count must be a positive number.')
        if(v <= 0 | !is.numeric(v)) stop('animalSpeed must be a positive number.')
        if(tm <= 0 | !is.numeric(tm)) stop('Time, t, must be a positive number.')

        # Calculate profile width, then density.
        p <- calcProfileWidth(alpha, theta, r)
        D <- count / {v * tm * p}
        return(D)
}


#' Calculate abudance using the gREM.
#'
#' Calculate abundance from count data using the gREM. 
#'
#' It is necessary to define a study area to calculate abundance rather than
#'    density. The easiest way to do this is to define the study area before
#'    data collection and then place camera traps/acoustic detectors/etc.
#'    randomly inside the area. It is more difficult to define the area
#'    studied post hoc.
#'
#' Note that the REM of Rowcliffe et al. 2008 can be used by setting alpha to 
#'    zero and the Gas Model of Yapp 1956 can be used setting alpha and theta
#'    to zero.
#'
#'@references Lucas, T. C. D., Moorcroft, E. A., Freeman, R., Rowcliffe, J. M.,
#'    Jones, K. E. (2015), A generalised random encounter model for estimating 
#'    animal density with remote sensor data. Methods in Ecology and Evolution. 
#'    doi: 10.1111/2041-210X.12346
#'
#' Rowcliffe, J., Field, J., Turvey, S. & Carbone, C. (2008) Estimating animal 
#'    density using camera traps without the need for individual recognition. 
#'    Journal of Applied Ecology, 45, 1228-1236.
#'
#' Yapp, W. (1956) The theory of line transects. Bird Study, 3, 93-104.
#'@seealso  \code{\link{gremDensity}}
#'@param count Number of detections.
#'@param alpha Call width in radians.
#'@param theta Detector width in radians.
#'@param r Sensor detection radius in metres.
#'@param v Average animal speed in metres per second.
#'@param tm Total survey time. This is the amount of time the sensors are
#'    active multiplied by the number of sensors used.
#'@param area The size of the study area in metres squared.
#'@name gremAbundance
#'@export

gremAbundance <- function(count, alpha, theta, r, v, tm, area){
        if(area <= 0 | !is.numeric(area)) stop('Area must be a positive number')
        D <- gremDensity(count, alpha, theta, r, v, tm)
        A <- D*area
        return(A)
}
