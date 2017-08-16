########################################################################################################
####
#### Run an iteration of a simulation with the vectorized version of 'candidate_iter_vec'
####


# SEED <- 1; n <- n1; allocation_share <- 0.50; age_share <- 100


sim_iter_vec <- function(SEED,
                         n,
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
                         allocation_share,    # this is the amount of sharing within allocation priority
                         age_share,           # this is the amount of sharing within age
                         allocation_effect,
                         age_effect,
                         diabetes_effect){
  
  
  # setting the seed...
  set.seed(SEED)
  
  
  # Need to generate data for the survival models
  
  donor_quality      <- runif(n, min = -50, max = 50)
  candidate_age      <- runif(n, min = -2, max = 2)
  candidate_diabetes <- rbinom(n, 1, prob = 0.35)
  
  post_dd_scale  <- dd_bl * exp(c_age_effect * candidate_age + 
                                  c_diabetes_effect * candidate_diabetes + 
                                  d_quality_effect * donor_quality)
  post_ld_scale  <- ld_bl
  post_rem_scale <- rem_bl * exp(c_age_effect * candidate_age + 
                                   c_diabetes_effect * candidate_diabetes)
  
  
  post_dd_time_pre  <- rexp(n, post_dd_scale)
  post_ld_time_pre  <- rexp(n, post_ld_scale)
  post_rem_time_pre <- rexp(n, post_rem_scale)
  
  
  post_dd_time  <- ifelse(post_dd_time_pre  < 5, post_dd_time_pre, 5L)
  post_ld_time  <- ifelse(post_ld_time_pre  < 5, post_ld_time_pre, 5L)
  post_rem_time <- ifelse(post_rem_time_pre < 5, post_rem_time_pre, 5L)
  
  post_dd_event  <- ifelse(post_dd_time_pre  < 5, 1, 0)
  post_ld_event  <- ifelse(post_ld_time_pre  < 5, 1, 0)
  post_rem_event <- ifelse(post_rem_time_pre < 5, 1, 0)
  
  
  # creating data frames of the vectors so that it is easier to use 'predict' functions
  post_dd_data <- data.frame(time = post_dd_time,
                             event = post_dd_event,
                             age = candidate_age,
                             diabetes = candidate_diabetes,
                             d_quality = donor_quality)
  
  post_ld_data <- data.frame(time = post_ld_time,
                             event = post_ld_event)
  
  post_rem_data  <- data.frame(time = post_rem_time,
                               event = post_rem_event,
                               age = candidate_age,
                               diabetes = candidate_diabetes)
  
  
  
  post_dd_model  <- coxph(Surv(time, event) ~ age + factor(diabetes) + d_quality, data = post_dd_data)
  post_ld_model  <- coxph(Surv(time, event) ~ 1, data = post_ld_data)
  post_rem_model <- coxph(Surv(time, event) ~ age + factor(diabetes), data = post_rem_data)
  
  
  
  ### Need to modify the baseline parameters based on randomly generated covariates
  x_diabetes   <- rbinom(n, 1, prob = 0.35)
  x_age        <- runif(n, min = -2, max = 2)
  x_allocation <- runif(n, min = -2, max = 2)
  d_quality    <- runif(n, min = -50, max = 50)
  
  
  
  dd_scale  <- dd_scale_bl * exp(allocation_effect * x_allocation + age_effect * x_age)
  ld_scale  <- ld_scale_bl
  rem_scale <- rem_scale_bl * exp(diabetes_effect * x_diabetes + age_effect * x_age)
  wlm_scale <- wlm_scale_bl * exp(diabetes_effect* x_diabetes + age_effect * x_age)
  
  
  
  
  
  dd_time  <- rweibull(n, shape = dd_shape, scale = dd_scale)
  ld_time  <- rweibull(n, shape = ld_shape, scale = ld_scale)
  rem_time <- rweibull(n, shape = rem_shape, scale = rem_scale)
  wlm_time <- rweibull(n, shape = wlm_shape, scale = wlm_scale)
  
  
  
  times <- apply(cbind(dd_time, ld_time, rem_time, wlm_time), 1, min)
  times <- times * {times < 3} + 3 * {times >= 3}
  event <- 1 * {{dd_time < apply(cbind(ld_time, rem_time, wlm_time), 1, min)} & {times < 3}} +
    2 * {{ld_time < apply(cbind(dd_time, rem_time, wlm_time), 1, min)} & {times < 3}} +
    3 * {{rem_time < apply(cbind(dd_time, ld_time, wlm_time), 1, min)} & {times < 3}} +
    4 * {{wlm_time < apply(cbind(dd_time, ld_time, rem_time), 1, min)} & {times < 3}}
  
  
  temp_data <- data.frame(times = times, 
                          event_factor = factor(event),
                          diabetes = x_diabetes,
                          age = x_age,
                          allocation = x_allocation,
                          donor_quality = d_quality)
  
  
  
  result <- sapply(1:nrow(temp_data), function(x){ 
    candidate_iter_vec(temp_data, x, 3,
                       allocation_share,
                       age_share,
                       post_dd_model,
                       post_ld_model,
                       post_rem_model,
                       c_age_effect,
                       c_diabetes_effect,
                       d_quality_effect,
                       dd_shape,
                       dd_scale,
                       ld_shape,
                       ld_scale,
                       rem_shape,
                       rem_scale,
                       wlm_shape,
                       wlm_scale) 
  })
  
  return(result)
}

