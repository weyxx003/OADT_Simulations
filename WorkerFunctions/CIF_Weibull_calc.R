
########################################################################################################
####
#### Obtaining 'true' value for a component of the CIF at a given point



CIF_Weibull_calc <- function(shape1, scale1, shape2, scale2, shape3, scale3, shape4, scale4, TheCIF, t_point){

	current_scale <- ifelse(TheCIF == 1, scale1,
	                        ifelse(TheCIF == 2, scale2,
	                               ifelse(TheCIF == 3, scale3, 
	                                      ifelse(TheCIF == 4, scale4, 666))))
	

	current_shape <- ifelse(TheCIF == 1, shape1,
	                        ifelse(TheCIF == 2, shape2,
	                               ifelse(TheCIF == 3, shape3, 
	                                      ifelse(TheCIF == 4, shape4, 666))))
	
	fnl <- integrate(function(y){
	  sapply(y, function(y){
	    current_shape / current_scale * ((y / current_scale)) ^ (current_shape - 1) * 
	      exp(-integrate(function(x){
	        (shape1 / scale1 * ((x / scale1) ^ (shape1 - 1)) +
	           shape2 / scale2 * ((x / scale2) ^ (shape2 - 1)) +
	           shape3 / scale3 * ((x / scale3) ^ (shape3 - 1)) +
	           shape4 / scale4 * ((x / scale4) ^ (shape4 - 1))) },
	        lower = 0, upper = y)$value)
	    })
	  }, 
	  lower = 0, upper = t_point)$value
	
	return(fnl)
}
