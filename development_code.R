##############################################################################################################
##############################################################################################################
######
###### Topic: Development of simulation code for evaluating the performance of the offer acceptance decision
######		tool (OADT).  This script aims to develop the framework before moving individual components
######		to seperate functions and actually running the simulations.
######
##############################################################################################################
##############################################################################################################


library(survival)


source("A:/SRTR/CIF_Simulations/CIF_Weibull_calc.R")
source("A:/SRTR/CIF_Simulations/candidate_iter.R")
source("A:/SRTR/CIF_Simulations/sim_iter.R")



# the sample size for the simulation
n <- 10000


# I need to simulate survival times for the individual components of the CIF.  We will work with Weibull distributions
#	since the likelihood of a removal reason usually increases (DD Tx and WL removal/mortality) or decreases (LD tx)

# decreased donor transplant baseline parameters
dd_shape    <- 1.05
dd_scale_bl <- 5

# living donor transplant baseline parameters
ld_shape    <- 0.95
ld_scale_bl <- 20

# removal baseline parameters
rem_shape    <- 1
rem_scale_bl <- 20

# wl mortality baseline parameters
wlm_shape    <- 1.05
wlm_scale_bl <- 20


########################################################################################################################
#
# Need to specify the effect of candidate age and candidate allocation priority on the likelihood of a deceased donor
#	transplant.  Similarly, need to specify the effect of diabetes on the likelihood of removal and wl mortality
#

diabetes_effect   <- 0.35
allocation_effect <- -0.50
age_effect        <- 0.10


########################################################################################################################
#
# Need to specify the baseline hazard and covariate effects for the posttransplant survival models and the
#	covariate effects (assuming an exponential distributed survival functions)
#

dd_bl   <- -log(0.850) / 3
ld_bl   <- -log(0.925) / 3
rem_bl  <- -log(0.700) / 3


c_age_effect      <- log(2) / 4
c_diabetes_effect <- log(1.25)
d_quality_effect  <- 4 / 100


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

n1 <- 250
n2 <- 500





test1_start <- proc.time()
test1 <- sim_iter(1,
                  n1,
                  dd_bl,
                  ld_bl,
                  rem_bl,
                  c_age_effect,
                  c_diabetes_effect,
                  d_quality_effect,
                  dd_shape,
                  dd_scale_bl,
                  ld_shape,
                  ld_scale_bl,
                  rem_shape,
                  rem_scale_bl,
                  wlm_shape,
                  wlm_scale_bl,
                  0.50,
                  100,
                  allocation_effect,
                  age_effect,
                  diabetes_effect)
test1_end <- proc.time()


test2_start <- proc.time()
test2 <- sim_iter(1,
                  n2,
                  dd_bl,
                  ld_bl,
                  rem_bl,
                  c_age_effect,
                  c_diabetes_effect,
                  d_quality_effect,
                  dd_shape,
                  dd_scale_bl,
                  ld_shape,
                  ld_scale_bl,
                  rem_shape,
                  rem_scale_bl,
                  wlm_shape,
                  wlm_scale_bl,
                  0.50,
                  100,
                  allocation_effect,
                  age_effect,
                  diabetes_effect)
test2_end <- proc.time()







test3_start <- proc.time()
test3 <- sim_iter_vec(1,
                      n1,
                      dd_bl,
                      ld_bl,
                      rem_bl,
                      c_age_effect,
                      c_diabetes_effect,
                      d_quality_effect,
                      dd_shape,
                      dd_scale_bl,
                      ld_shape,
                      ld_scale_bl,
                      rem_shape,
                      rem_scale_bl,
                      wlm_shape,
                      wlm_scale_bl,
                      0.50,
                      100,
                      allocation_effect,
                      age_effect,
                      diabetes_effect)
test3_end <- proc.time()







test4_start <- proc.time()
test4 <- sim_iter_vec(1,
                      n2,
                      dd_bl,
                      ld_bl,
                      rem_bl,
                      c_age_effect,
                      c_diabetes_effect,
                      d_quality_effect,
                      dd_shape,
                      dd_scale_bl,
                      ld_shape,
                      ld_scale_bl,
                      rem_shape,
                      rem_scale_bl,
                      wlm_shape,
                      wlm_scale_bl,
                      0.50,
                      100,
                      allocation_effect,
                      age_effect,
                      diabetes_effect)
test4_end <- proc.time()






(test1_end[3] - test1_start[3]) / 60
(test2_end[3] - test2_start[3]) / 60
(test3_end[3] - test3_start[3]) / 60













