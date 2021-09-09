# -----------------------------------------------------------------------------
# 2021-01-03 SL library set up
# -----------------------------------------------------------------------------
grid_params <- list(max_depth = c(2, 4, 6, 8),
                    eta = c(0.001, 0.01, 0.1, 0.2, 0.3),
                    nrounds = c(20, 50))
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
params_default <- list(nthread = getOption("sl.cores.learners", 1))
xgb_learners <- apply(grid, MARGIN = 1, function(params_tune) {
  do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))})

stack <- make_learner_stack('Lrnr_glm_fast', 
                            'Lrnr_mean',
                            'Lrnr_gam',
                            'Lrnr_ranger',
                            'Lrnr_glmnet')

fancy_stack <- make_learner(Stack, 
                            xgb_learners[[1]],
                            xgb_learners[[10]],
                            xgb_learners[[20]],
                            xgb_learners[[30]],
                            xgb_learners[[40]],
                            stack)

# -----------------------------------------------------------------------------
# CV-TMLE, TMLE, IPW, AIPW, CV-IPW, CV-AIPW
# -----------------------------------------------------------------------------
# inputs: 
#       data is a dataframe, the simulated data
#       node_list is a list consist of W, A, Y, each contains the variable(s) names
#       gbound the lower bound for propensity score
run_6est <- function(data, node_list, gbound = 0.025){
  
  tmle_spec = tmle_ATE(treatment_level = 1, control_level = 0)
  tmle_task <- tmle_spec$make_tmle_task(data, node_list)
  
  lrnr_sl_a <- make_learner(Lrnr_sl,
                            learners = fancy_stack,
                            outcome_type = 'binomial',
                            metalearner= make_learner(Lrnr_solnp,
                                                      learner_function = metalearner_logistic_binomial,
                                                      loss_function = loss_loglik_binomial))
  lrnr_sl_y <- make_learner(Lrnr_sl,
                            learners = fancy_stack,
                            outcome_type = 'continuous',
                            metalearner = make_learner(Lrnr_nnls))
  
  learner_list <- list(Y = lrnr_sl_y, A = lrnr_sl_a)
  
  bounded_factor_list <- list(
    define_lf(LF_emp, "W"),
    define_lf(LF_fit, "A", learner = learner_list[["A"]],bound = gbound),#
    define_lf(LF_fit, "Y", learner = learner_list[["Y"]], type = 'mean')
  )
  
  bounded_likelihood_def <- Likelihood$new(bounded_factor_list)
  bounded_likelihood <- bounded_likelihood_def$train(tmle_task)
  
  # 1. compute CV-TMLE
  targeted_likelihood <- Targeted_Likelihood$new(bounded_likelihood)
  tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
  
  cvtmle_fit <- fit_tmle3(tmle_task, targeted_likelihood,
                          tmle_params,targeted_likelihood$updater)
  
  # 2. compute TMLE
  targeted_likelihood_no_cv <- Targeted_Likelihood$new(bounded_likelihood, 
                                                       updater = list(cvtmle = FALSE))
  tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood_no_cv)
  
  tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood_no_cv, 
                        tmle_params,targeted_likelihood_no_cv$updater)
  
  # extract cv Q(1,W), Q(0,W) and g(1|W) to compute cv-ipw and cv-aipw
  ate <- tmle_params[[1]]
  cf_task1 <- ate$cf_likelihood_treatment$cf_tasks[[1]]
  cf_task0 <- ate$cf_likelihood_control$cf_tasks[[1]]
  
  QbarAW_cv <- bounded_likelihood$get_likelihood(tmle_task, "Y", fold_number="validation")
  Qbar1W_cv <- bounded_likelihood$get_likelihood(cf_task1, "Y", fold_number="validation")
  Qbar0W_cv <- bounded_likelihood$get_likelihood(cf_task0, "Y", fold_number="validation")
  g1W_cv <- bounded_likelihood$get_likelihood(cf_task1, "A", fold_number="validation")
  
  
  # 3. CV-IPW(Hajek)
  a <- as.numeric(data$a)
  y <- as.numeric(data$haz)
  psi_cvipw <- mean(a*y/g1W_cv)/mean(a/g1W_cv) - mean((1-a)*y/(1-g1W_cv))/mean((1-a)/(1-g1W_cv))
  ic_cvipw <- y * ((a/g1W_cv) - (1 - a)/(1 - g1W_cv)) - psi_cvipw
  
  # 4. CV-AIPW
  psi_cvaipw <- mean(a*(y-Qbar1W_cv)/g1W_cv + Qbar1W_cv) - mean((1-a)*(y-Qbar0W_cv)/(1-g1W_cv) + Qbar0W_cv)
  ic_cvaipw <- (y - QbarAW_cv) * ((a/g1W_cv) - (1 - a)/(1 - g1W_cv)) + (Qbar1W_cv - Qbar0W_cv) - psi_cvaipw
  
  # extract Q(1,W), Q(0,W) and g(1|W) to compute ipw and aipw
  QbarAW <- bounded_likelihood$get_likelihood(tmle_task, "Y")
  Qbar1W <- bounded_likelihood$get_likelihood(cf_task1, "Y")
  Qbar0W <- bounded_likelihood$get_likelihood(cf_task0, "Y")
  g1W <- bounded_likelihood$get_likelihood(cf_task1, "A")
  
  # 5. IPW(Hajek)
  psi_ipw <- mean(a*y/g1W)/mean(a/g1W) - mean((1-a)*y/(1-g1W))/mean((1-a)/(1-g1W))
  ic_ipw <- y * ((a/g1W) - (1 - a)/(1 - g1W)) - psi_ipw
  
  # 6. AIPW
  psi_aipw <- mean(a*(y-Qbar1W)/g1W + Qbar1W) - mean((1-a)*(y-Qbar0W)/(1-g1W) + Qbar0W)
  ic_aipw <- (y - QbarAW) * ((a/g1W) - (1 - a)/(1 - g1W)) + (Qbar1W - Qbar0W) - psi_aipw
  
  # save 
  # output_g1_hat <- paste0('~/Repo/...', "g1_hat_", studies[i], '.csv')
  # output_QA_hat <- paste0('~/Repo/...', "QA_hat_", studies[i], '.csv')
  # write.csv(QbarAW, output_QA_hat)
  # write.csv(g1W, output_g1_hat)
  
  return(list("tmle_fit" = tmle_fit,
              "cvtmle_fit" = cvtmle_fit,
              "psi_ipw" = psi_ipw,
              "psi_aipw" = psi_aipw,
              "psi_cvipw" = psi_cvipw,
              "psi_cvaipw" = psi_cvaipw,
              "ic_ipw" = ic_ipw,
              "ic_aipw" = ic_aipw,
              "ic_cvipw" = ic_cvipw,
              "ic_cvaipw" = ic_cvaipw,
              "Qbar1W" = Qbar1W,
              "Qbar0W" = Qbar0W))
}

# -----------------------------------------------------------------------------
# Run All Estimators
# -----------------------------------------------------------------------------
run_all_simu <- function(B, df, ss_estimate, nodes, w, a, y, gbound){
  results_cols <- c('i', 'ss', 'tmle', 'tmle_se', 'tmle_lower', 
                    'tmle_upper', 'cvtmle', 'cvtmle_se',
                    'cvtmle_lower', 'cvtmle_upper', 'iptw',
                    'iptw_se', 'iptw_lower', 'iptw_upper',
                    'gcomp',
                    'aiptw', 'aiptw_se', 'aiptw_lower', 'aiptw_upper',
                    'ctmle', 'ctmle_se', 'ctmle_lower',
                    'ctmle_upper', 'drop_cov', 
                    'cviptw', 'cviptw_se', 'cviptw_lower', 'cviptw_upper',
                    'cvaiptw', 'cvaiptw_se', 'cvaiptw_lower', 'cvaiptw_upper')
  
  results_df <- data.frame(matrix(NA, nrow = B, ncol = length(results_cols)))
  colnames(results_df) <- results_cols
  
  run_bootstrap <- foreach(b = 1:B, .combine = 'rbind') %do% {
    
    print("may the power be with you")
    
    # define variables
    w <- colnames(df) %w/o% c('haz', 'a')
    a <- 'a'
    y <- 'haz'
    nodes <- list(W = w,
                  A = a,
                  Y = y)
    
    results_df_row <- data.frame(matrix(NA, nrow = 1, ncol = length(results_cols)))
    colnames(results_df_row) <- results_cols
    #set.seed(123 + b)
    
    # Generate bootstrap data
    simu_output <- get_simu_data(df, w, a, y, ss_estimate$g_fit,
                                 ss_estimate$Q_fit, ss_estimate$rv)
    simu_data <- simu_output$data
    results_df_row$i <- b
    results_df_row$ss <- ss_estimate$psi
    
    # redefine variables in case variable drop in simulation
    w <- colnames(simu_data) %w/o% c('haz', 'a')
    nodes <- list(W = w,
                  A = 'a',
                  Y = 'haz')
    
    # run 
    allres <- run_6est(data = simu_data, 
                       node_list = nodes, 
                       gbound = gbound)

    # TMLE
    tmle_output <- allres$tmle_fit
    results_df_row$tmle <- tmle_output$summary$tmle_est
    results_df_row$tmle_se <- tmle_output$summary$se
    results_df_row$tmle_lower <- tmle_output$summary$lower
    results_df_row$tmle_upper <- tmle_output$summary$upper
    
    # G-computation
    results_df_row$gcomp <- tmle_output$summary$init_est
    
    rm(tmle_output)
    
    # CVTMLE 
    cvtmle_output <- allres$cvtmle_fit
    results_df_row$cvtmle <- cvtmle_output$summary$tmle_est
    results_df_row$cvtmle_se <- cvtmle_output$summary$se
    results_df_row$cvtmle_lower <- cvtmle_output$summary$lower
    results_df_row$cvtmle_upper <- cvtmle_output$summary$upper
    rm(cvtmle_output)
    
    # C-TMLE 
    q = cbind(allres$Q1, allres$Q0)
    y= simu_data$haz
    a= simu_data$a
    w <- simu_data %>%  select (-c('a', 'haz')) 
    
    ctmle_output <- ctmleDiscrete(Y = y, A = a, W = w, Q = q,
                                  preOrder = FALSE, detailed = TRUE, gbound = gbound)
    
    results_df_row$ctmle <- ctmle_output$est
    results_df_row$ctmle_se <-sqrt(ctmle_output$var.psi)  
    results_df_row$ctmle_lower <- ctmle_output$CI[1]
    results_df_row$ctmle_upper <- ctmle_output$CI[2]
    rm(ctmle_output)
    
    
    # 0413 Fixed AIPW and Stabilized IPW
    psi_ipw <- allres$psi_ipw
    psi_aipw <- allres$psi_aipw
    
    # IC for AIPTW to get SE and 95% confidence intervals
    ic_aiptw <- allres$ic_aipw
    se_aiptw <- sd(ic_aiptw)
    
    results_df_row$aiptw <- psi_aipw
    results_df_row$aiptw_se <- se_aiptw
    results_df_row$aiptw_upper <- psi_aipw + (1.96 * se_aiptw / sqrt(length(ic_aiptw)))
    results_df_row$aiptw_lower <- psi_aipw - (1.96 * se_aiptw / sqrt(length(ic_aiptw)))
    
    # IC for IPTW to get SE and 95% confidence intervals
    ic_iptw <- allres$ic_ipw
    se_iptw <- sd(ic_iptw)
    
    results_df_row$iptw <- psi_ipw
    results_df_row$iptw_se <- se_iptw
    results_df_row$iptw_upper <- psi_ipw + (1.96 * se_iptw / sqrt(length(ic_iptw)))
    results_df_row$iptw_lower <- psi_ipw - (1.96 * se_iptw / sqrt(length(ic_iptw)))
    
    
    # 0414 CV AIPW and Stabilized IPW !!!
    psi_cvipw <- allres$psi_cvipw
    psi_cvaipw <- allres$psi_cvaipw
    
    # IC for AIPTW to get SE and 95% confidence intervals
    ic_cvaiptw <- allres$ic_cvaipw
    se_cvaiptw <- sd(ic_cvaiptw)
    
    results_df_row$cvaiptw <- psi_cvaipw
    results_df_row$cvaiptw_se <- se_cvaiptw
    results_df_row$cvaiptw_upper <- psi_cvaipw + (1.96 * se_cvaiptw / sqrt(length(ic_cvaiptw)))
    results_df_row$cvaiptw_lower <- psi_cvaipw - (1.96 * se_cvaiptw / sqrt(length(ic_cvaiptw)))
    
    # IC for IPTW to get SE and 95% confidence intervals
    ic_cviptw <- allres$ic_cvipw
    se_cviptw <- sd(ic_cviptw)
    
    results_df_row$cviptw <- psi_cvipw
    results_df_row$cviptw_se <- se_cviptw
    results_df_row$cviptw_upper <- psi_cvipw + (1.96 * se_cviptw / sqrt(length(ic_cviptw)))
    results_df_row$cviptw_lower <- psi_cvipw - (1.96 * se_cviptw / sqrt(length(ic_cviptw)))
    
    # print(se_cvaiptw)
    # print(se_aiptw)
    return_list <- c('i' = results_df_row$i, 'ss' = results_df_row$ss, 
                     'tmle' = results_df_row$tmle, 'tmle_se' = results_df_row$tmle_se, 
                     'tmle_lower' = results_df_row$tmle_lower, 
                     'tmle_upper' = results_df_row$tmle_upper, 
                     'cvtmle' = results_df_row$cvtmle, 
                     'cvtmle_se' = results_df_row$cvtmle_se,
                     'cvtmle_lower' = results_df_row$cvtmle_lower, 
                     'cvtmle_upper' = results_df_row$cvtmle_upper, 
                     'iptw' = results_df_row$iptw,
                     'iptw_se' = results_df_row$iptw_se, 
                     'iptw_lower' = results_df_row$iptw_lower, 
                     'iptw_upper' = results_df_row$iptw_upper,
                     'gcomp' = results_df_row$gcomp,
                     'aiptw' = results_df_row$aiptw , 
                     'aiptw_se' = results_df_row$aiptw_se, 
                     'aiptw_lower' = results_df_row$aiptw_lower, 
                     'aiptw_upper' = results_df_row$aiptw_upper,
                     'ctmle' = results_df_row$ctmle, 
                     'ctmle_se' = results_df_row$ctmle_se,
                     'ctmle_lower' = results_df_row$ctmle_lower, 
                     'ctmle_upper' = results_df_row$ctmle_upper,
                     'drop_cov'= bootstrap_output$drop_cov,
                     'cviptw' = results_df_row$cviptw,
                     'cviptw_se' = results_df_row$cviptw_se, 
                     'cviptw_lower' = results_df_row$cviptw_lower, 
                     'cviptw_upper' = results_df_row$cviptw_upper,
                     'cvaiptw' = results_df_row$cvaiptw , 
                     'cvaiptw_se' = results_df_row$cvaiptw_se, 
                     'cvaiptw_lower' = results_df_row$cvaiptw_lower, 
                     'cvaiptw_upper' = results_df_row$cvaiptw_upper)
    
    return(return_list)
  }
  
  results_df <- as.data.frame(run_bootstrap)
  results_df[, 1] <- sub('.*\\.', '', results_df[, 1])
  
  traceback()
  gc()
  return(results_df)
}




