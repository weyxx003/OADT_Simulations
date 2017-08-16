########################################################################################################
####
#### Estimating the post-offer patient mortality
####


# dataset <- temp_data; pt_index <- x; evaluated_t <- 3; allocation_sharing <- allocation_share; age_sharing <- age_share
# dd_model <- post_dd_model; ld_model <- post_ld_model; rem_model <- post_rem_model; cand_age_effect <- c_age_effect
# cand_diabetes_effect <- c_diabetes_effect; don_quality_effect <- d_quality_effect; true_dd_shape <- dd_shape
# true_dd_scale <- dd_scale; true_ld_shape <- ld_shape; true_ld_scale <- ld_scale; true_rem_shape <- rem_shape
# true_rem_scale <- rem_scale; true_wlm_shape <- wlm_shape; true_wlm_scale <- wlm_scale


candidate_iter <- function(dataset, 
                           pt_index, 
                           evaluated_t, 
                           allocation_sharing,
                           age_sharing,
                           dd_model,
                           ld_model,
                           rem_model,
                           cand_age_effect,
                           cand_diabetes_effect,
                           don_quality_effect,
                           true_dd_shape,
                           true_dd_scale,
                           true_ld_shape,
                           true_ld_scale,
                           true_rem_shape,
                           true_rem_scale,
                           true_wlm_shape,
                           true_wlm_scale){
  
  
  # 'tmpD1' provides the characteristics of the given candidate
	tmpD1 <- dataset[pt_index,,drop = FALSE]


	# 'tmp_data1' is the set of offers within the locality of the candidate
	tmp_data1 <- subset(dataset, 
	                    {allocation >= (tmpD1$allocation - allocation_sharing)} & 
	                      {allocation <= (tmpD1$allocation + allocation_sharing)} & 
	                      {diabetes == tmpD1$diabetes})

	
	# Adaptively choosing the age cutoff to select the nearest 100 offers in terms of age (results in 'tmp_data2')
	abs_age_diff <- abs(tmp_data1$age - tmpD1$age)
	
	if(length(abs_age_diff) >= age_sharing){ abs_age_cutoff <- abs_age_diff[order(abs_age_diff)][age_sharing] }
	if(length(abs_age_diff) <  age_sharing){ abs_age_cutoff <- max(abs_age_diff) }
	
	
	
	tmp_data2 <- subset(tmp_data1, 
	                    {age >= (tmpD1$age - abs_age_cutoff)} & 
	                      {age <= (tmpD1$age + abs_age_cutoff)})


	# Estimating the CIF for the candidate
	tmp_cif <- survfit(Surv(times, event_factor) ~ 1, data = tmp_data2)


	# points in the CIF the require evaluation
	cif_ts_unique <- unique(tmp_cif$time)
	cif_ts        <- cif_ts_unique[order(cif_ts_unique)]
	cif_ts        <- cif_ts[cif_ts <= evaluated_t]



	########################################################################################################
	####
	#### Setup for the prediction DD/LD/REM survival
	####

	# creating the data sets to obtain survival estimates
	dd_pred_data <- data.frame(age = rep(tmpD1$age, sum(tmp_data2$event_factor == 1)),
	                           diabetes = rep(factor(tmpD1$diabetes), sum(tmp_data2$event_factor == 1)),
	                           d_quality = tmp_data2$donor_quality[tmp_data2$event_factor == 1])


	dd_pred_data_median <- data.frame(age = tmpD1$age,
	                                  diabetes = factor(tmpD1$diabetes),
	                                  d_quality = median(tmp_data2$donor_quality[tmp_data2$event_factor == 1]))


	rem_pred_data <- data.frame(age = tmpD1$age,
	                            diabetes = factor(tmpD1$diabetes))


	rem_pred_data_median <- data.frame(age = tmpD1$age,
	                                   diabetes = factor(tmpD1$diabetes))

	########################################################################################################
	####
	#### Estimating the survival functions for the averaged effects
	####

	dd_preds  <- survfit(dd_model, newdata = dd_pred_data, se.fit = FALSE)
	ld_preds  <- survfit(ld_model)
	rem_preds <- survfit(rem_model, newdata = rem_pred_data, se.fit = FALSE)


	dd_surv_preds <- sapply(cif_ts, function(x){
	  if({sum(dd_preds$time <= (evaluated_t - x)) > 0} & {sum(as.numeric(as.character(tmp_data2$event_factor)) == 1) > 1}){ 
	    fnl <- mean(dd_preds$surv[sum(dd_preds$time <= (evaluated_t - x)), ]) 
	  }
	  if({sum(dd_preds$time <= (evaluated_t - x)) > 0} & {sum(as.numeric(as.character(tmp_data2$event_factor)) == 1) == 1}){ 
	    fnl <- mean(dd_preds$surv[sum(dd_preds$time <= (evaluated_t - x))])
	  }
	  if(sum(dd_preds$time <= (evaluated_t - x)) == 0){ fnl <- 1 }
	  return(fnl)
	  })

	ld_surv_preds <- sapply(cif_ts, function(x){
	  if(sum(ld_preds$time <= (evaluated_t - x)) >  0){ 
	    fnl <- ld_preds$surv[sum(ld_preds$time <= (evaluated_t - x))] 
	  }
	  if(sum(ld_preds$time <= (evaluated_t - x)) == 0){ fnl <- 1 }
	  return(fnl)
	  })

	rem_surv_preds <- sapply(cif_ts, function(x){
	  if(sum(rem_preds$time <= (evaluated_t - x)) >  0){ 
	    fnl <- rem_preds$surv[sum(rem_preds$time <= (evaluated_t - x))]
	  }
	  if(sum(rem_preds$time <= (evaluated_t - x)) == 0){ fnl <- 1 }
	  return(fnl)
	  })
	
	########################################################################################################
	####
	#### Estimating the survival functions for the true effects

	dd_surv_true <- sapply(cif_ts, function(x){
	  txs <- sapply(1:nrow(dd_pred_data), function(y){ 
	    dd_bl * 
	      exp(cand_age_effect * tmpD1$age + 
	            cand_diabetes_effect * as.numeric(as.character(tmpD1$diabetes)) + 
	            don_quality_effect * dd_pred_data$d_quality[y]) })
	  fnl <- mean(pexp((evaluated_t - x), txs, lower.tail = FALSE))
	  return(fnl)
	  })


	ld_surv_true <- sapply(cif_ts, function(x){
	  fnl <- pexp((evaluated_t - x), ld_bl, lower.tail = FALSE)
	  return(fnl)
	  })


	rem_surv_true <- sapply(cif_ts, function(x){
	  txs <- rem_bl * exp(cand_age_effect * tmpD1$age + 
	                        cand_diabetes_effect * as.numeric(as.character(tmpD1$diabetes)))
	  fnl <- pexp((evaluated_t - x), txs, lower.tail = FALSE)
	  return(fnl)
	  })

	########################################################################################################
	####
	#### Estimating the survival functions for the median effects

	dd_preds_median  <- survfit(dd_model, newdata = dd_pred_data_median, se.fit = FALSE)
	ld_preds_median  <- survfit(ld_model)
	rem_preds_median <- survfit(rem_model, newdata = rem_pred_data_median, se.fit = FALSE)


	dd_surv_preds_median <- sapply(cif_ts, function(x){
	  if(sum(dd_preds_median$time <= (evaluated_t - x)) >  0){ 
	    fnl <- dd_preds_median$surv[sum(dd_preds_median$time <= (evaluated_t - x))]
	  }
	  if(sum(dd_preds_median$time <= (evaluated_t - x)) == 0){ fnl <- 1 }
	  return(fnl)
	  })

	ld_surv_preds_median <- sapply(cif_ts, function(x){
	  if(sum(ld_preds_median$time <= (evaluated_t - x)) >  0){ 
	    fnl <- ld_preds_median$surv[sum(ld_preds_median$time <= (evaluated_t - x))] 
	  }
	  if(sum(ld_preds_median$time <= (evaluated_t - x)) == 0){ fnl <- 1 }
	  return(fnl)
		})

	rem_surv_preds_median <- sapply(cif_ts, function(x){
	  if(sum(rem_preds_median$time <= (evaluated_t - x)) >  0){ 
	    fnl <- rem_preds_median$surv[sum(rem_preds_median$time <= (evaluated_t - x))] 
	  }
	  if(sum(rem_preds_median$time <= (evaluated_t - x)) == 0){ fnl <- 1 }
	  return(fnl)
	  })


	########################################################################################################
	####
	#### Obtaining the estimated and 'true' quantities


	preds_all    <- 0
	preds_median <- 0
	preds_true   <- 0
	cif_ests     <- NULL
	for(i in 1:length(cif_ts)){
	  
	  if(i == 1){ 
	    p_dd  <- tmp_cif$pstate[1,1]
			p_ld  <- tmp_cif$pstate[1,2]
			p_rem <- tmp_cif$pstate[1,3]
			p_wlm <- tmp_cif$pstate[1,4]

	
			p_o_dd <- CIF_Weibull_calc(true_dd_shape, true_dd_scale[pt_index], 
			                           true_ld_shape, true_ld_scale, 
			                           true_rem_shape, true_rem_scale[pt_index], 
			                           true_wlm_shape, true_wlm_scale[pt_index], 
			                           1, cif_ts[i])
			
			p_o_ld <- CIF_Weibull_calc(true_dd_shape, true_dd_scale[pt_index], 
			                           true_ld_shape, true_ld_scale, 
			                           true_rem_shape, true_rem_scale[pt_index], 
			                           true_wlm_shape, true_wlm_scale[pt_index],
			                           2, cif_ts[i])

			p_o_rem <- CIF_Weibull_calc(true_dd_shape, true_dd_scale[pt_index], 
			                            true_ld_shape, true_ld_scale, 
			                            true_rem_shape, true_rem_scale[pt_index], 
			                            true_wlm_shape, true_wlm_scale[pt_index], 
			                            3, cif_ts[i])

			p_o_wlm <- CIF_Weibull_calc(true_dd_shape, true_dd_scale[pt_index], 
			                            true_ld_shape, true_ld_scale, 
			                            true_rem_shape, true_rem_scale[pt_index], 
			                            true_wlm_shape, true_wlm_scale[pt_index],
			                            4, cif_ts[i])
		}
		if(i > 1){
		  
			p_dd  <- tmp_cif$pstate[i,1] - tmp_cif$pstate[i - 1,1]
			p_ld  <- tmp_cif$pstate[i,2] - tmp_cif$pstate[i - 1,2]
			p_rem <- tmp_cif$pstate[i,3] - tmp_cif$pstate[i - 1,3]
			p_wlm <- tmp_cif$pstate[i,4] - tmp_cif$pstate[i - 1,4]


			p_o_dd <- CIF_Weibull_calc(true_dd_shape, true_dd_scale[pt_index], 
			                           true_ld_shape, true_ld_scale, 
			                           true_rem_shape, true_rem_scale[pt_index], 
			                           true_wlm_shape, true_wlm_scale[pt_index], 
			                           1, cif_ts[i]) - 
			  CIF_Weibull_calc(true_dd_shape, true_dd_scale[pt_index], 
			                   true_ld_shape, true_ld_scale, 
			                   true_rem_shape, true_rem_scale[pt_index], 
			                   true_wlm_shape, true_wlm_scale[pt_index],
			                   1, cif_ts[i - 1])

			p_o_ld <- CIF_Weibull_calc(true_dd_shape, true_dd_scale[pt_index], 
			                           true_ld_shape, true_ld_scale, 
			                           true_rem_shape, true_rem_scale[pt_index], 
			                           true_wlm_shape, true_wlm_scale[pt_index],
			                           2, cif_ts[i]) - 
			  CIF_Weibull_calc(true_dd_shape, true_dd_scale[pt_index], 
			                   true_ld_shape, true_ld_scale, 
			                   true_rem_shape, true_rem_scale[pt_index], 
			                   true_wlm_shape, true_wlm_scale[pt_index], 
			                   2, cif_ts[i - 1])

			p_o_rem <- CIF_Weibull_calc(true_dd_shape, true_dd_scale[pt_index], 
			                            true_ld_shape, true_ld_scale, 
			                            true_rem_shape, true_rem_scale[pt_index], 
			                            true_wlm_shape, true_wlm_scale[pt_index],
			                            3, cif_ts[i]) - 
			  CIF_Weibull_calc(true_dd_shape, true_dd_scale[pt_index], 
			                   true_ld_shape, true_ld_scale, 
			                   true_rem_shape, true_rem_scale[pt_index], 
			                   true_wlm_shape, true_wlm_scale[pt_index],
			                   3, cif_ts[i - 1])

			p_o_wlm <- CIF_Weibull_calc(true_dd_shape, true_dd_scale[pt_index], 
			                            true_ld_shape, true_ld_scale, 
			                            true_rem_shape, true_rem_scale[pt_index], 
			                            true_wlm_shape, true_wlm_scale[pt_index],
			                            4, cif_ts[i]) - 
			  CIF_Weibull_calc(true_dd_shape, true_dd_scale[pt_index], 
			                   true_ld_shape, true_ld_scale, 
			                   true_rem_shape, true_rem_scale[pt_index], 
			                   true_wlm_shape, true_wlm_scale[pt_index],
			                   4, cif_ts[i - 1])
		}

		preds_all_tmp    <- p_dd * (1 - dd_surv_preds[i]) + p_ld * (1 - ld_surv_preds[i]) + p_rem * (1 - rem_surv_preds[i]) + p_wlm
		preds_median_tmp <- p_dd * (1 - dd_surv_preds_median[i]) + p_ld * (1 - ld_surv_preds_median[i]) + 
						p_rem * (1 - rem_surv_preds_median[i]) + p_wlm
		preds_true_tmp   <- p_o_dd * (1 - dd_surv_true[i]) + p_o_ld * (1 - ld_surv_true[i]) + 
						p_o_rem * (1 - rem_surv_true[i]) + p_o_wlm

		preds_all    <- preds_all + preds_all_tmp
		preds_median <- preds_median + preds_median_tmp
		preds_true   <- preds_true + preds_true_tmp

		cif_ests <- rbind(cif_ests, c(p_dd, p_o_dd, p_ld, p_o_ld, p_rem, p_o_rem, p_wlm, p_o_wlm))
	}

	cif_sqr_err <- c(mean( (cif_ests[,1] - cif_ests[,2]) ^ 2),
	                 mean( (cif_ests[,3] - cif_ests[,4]) ^ 2),
	                 mean( (cif_ests[,5] - cif_ests[,6]) ^ 2),
	                 mean( (cif_ests[,7] - cif_ests[,8]) ^ 2))
	
	cif_abs_err <- c(mean( abs(cif_ests[,1] - cif_ests[,2]) ),
	                 mean( abs(cif_ests[,3] - cif_ests[,4]) ),
	                 mean( abs(cif_ests[,5] - cif_ests[,6]) ),
	                 mean( abs(cif_ests[,7] - cif_ests[,8]) ))
	
	cif_bias <- c(mean( cif_ests[,1] - cif_ests[,2] ),
	              mean( cif_ests[,3] - cif_ests[,4] ),
	              mean( cif_ests[,5] - cif_ests[,6] ),
	              mean( cif_ests[,7] - cif_ests[,8] ))
	
	surv_bias <- c(mean( dd_surv_preds - dd_surv_true ),
	               mean( ld_surv_preds - ld_surv_true ),
	               mean( rem_surv_preds - rem_surv_true ))
	
	surv_sqr_err <- c(mean( (dd_surv_preds - dd_surv_true) ^ 2 ),
	                  mean( (ld_surv_preds - ld_surv_true) ^ 2 ),
	                  mean( (rem_surv_preds - rem_surv_true) ^ 2 ))
	
	surv_abs_err <- c(mean( abs(dd_surv_preds - dd_surv_true) ),
	                  mean( abs(ld_surv_preds - ld_surv_true) ),
	                  mean( abs(rem_surv_preds - rem_surv_true) ))
	
	median_surv_bias <- c(mean( dd_surv_preds_median - dd_surv_true ),
	                      mean( ld_surv_preds_median - ld_surv_true ),
	                      mean( rem_surv_preds_median - rem_surv_true ))
	
	median_surv_sqr_err <- c(mean( (dd_surv_preds_median - dd_surv_true) ^ 2 ),
	                         mean( (ld_surv_preds_median - ld_surv_true) ^ 2 ),
	                         mean( (rem_surv_preds_median - rem_surv_true) ^ 2 ))
	
	median_surv_abs_err <- c(mean( abs(dd_surv_preds_median - dd_surv_true) ),
	                         mean( abs(ld_surv_preds_median - ld_surv_true) ),
	                         mean( abs(rem_surv_preds_median - rem_surv_true) ))
	
	fnl <- c(preds_all, preds_median, preds_true, cif_bias, cif_sqr_err, cif_abs_err,
	         surv_bias, surv_sqr_err, surv_abs_err, median_surv_bias, median_surv_sqr_err, median_surv_abs_err)
	
	
	return(fnl)
}

