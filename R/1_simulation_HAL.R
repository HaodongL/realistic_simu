# -----------------------------------------------------------------------------
# Undersoomthed HAL function
# -----------------------------------------------------------------------------
# inputs: 
#       df is a dataframe
#       yname is a string
#       xname is a string
#       Nlam is a scalar, number of candidates lambda
#       type is a string, type of y
under_HAL <- function(df, yname, xname, Nlam, type = "gaussian"){
  
  # variables
  y <- as.numeric(as.matrix(df %>% select(all_of(yname))))
  x <- df %>% select(all_of(xname)) %>% 
              mutate_if(sapply(., is.factor), as.numeric)
  n <- nrow(df)
  
  # get initial fit
  print("---initial fitting---")
  tic()
  fit <- fit_hal(X = x,
                 Y = y, 
                 return_x_basis = TRUE,
                 family = type,
                 num_knots = num_knots_generator(
                   max_degree = ifelse(ncol(x) >= 20, 2, 3),
                   smoothness_orders = 1,
                   base_num_knots_0 = 500,
                   base_num_knots_1 = max(100, ceiling(sqrt(n)))
                 )
                )
  toc()

  # only non-zero direction
  init_coef <-fit$coefs[-1]
  nonzero_col <- which(init_coef != 0)
  
  # if all coefs are zero, skip undersmooth and use the initial fit
  if (length(nonzero_col) == 0){
    res <- list("lambda_init" = fit$lambda_star,
                "lambda_under"= fit$lambda_star)
  }else{
    # refit on new lambda sequence
    print("---candidates fitting---")
    tic()
    us_lambda <- fit$lambda_star*10^seq(from=0, to=-3, length=Nlam)
    us_fit <- glmnet(fit$x_basis, y, lambda=us_lambda, family = type, standardize = FALSE)
    toc()
    
    # evaluate refits
    if (type == "gaussian"){
      pred_mat <- predict(us_fit, fit$x_basis)
      preds_init <- predict(fit, new_data = x)
    }else if (type == "binomial"){
      pred_mat <- predict(us_fit, fit$x_basis, type = "response")
      preds_init <- predict(fit, new_data = x) 
    }
    
    resid_mat <- pred_mat-y
    basis_mat <- as.matrix(fit$x_basis)
    basis_mat <- as.matrix(basis_mat[, nonzero_col])
    
    # estimates of sd in each direction using initial fit
    resid_init <- preds_init-y
    sd_est <- rep(NA, ncol(basis_mat))
    
    for (i in 1:ncol(basis_mat)){
      u <- basis_mat[,i]
      score_init <- resid_init * u
      sd_est[i] <- sd(score_init)
    }
    
    # check the criterion 
    print("---checking criterion---")
    tic()
    max_score <- get_maxscore(basis_mat = basis_mat, 
                              resid_mat = resid_mat,
                              sd_est = sd_est, 
                              Nlam = Nlam, us_fit = us_fit)
    toc()
    
    # get the first lambda that satisfies the criteria
    lambda_under <- us_lambda[max_score <= 1/(sqrt(n)*log(n))][1]
    
    # collect and save results
    coef_mat <- as.matrix(us_fit$beta)
    
    df_under <- data.frame("lambda" = NA,
                           "l1_norm" = NA,
                           "n_coef" = NA)
    
    for (j in 1:Nlam){
      df_under[j,1] = us_lambda[j]
      df_under[j,2] = sum(abs(coef_mat[,j]))
      df_under[j,3] = sum(coef_mat[,j] != 0)
    }
    
    res <- list("lambda_init" = fit$lambda_star,
                "lambda_under" = lambda_under,
                "df_under" = df_under)
  }
  return(res)
}

# # make plot
# plot(x = df_under$lambda, y = df_under$l1_norm)
# abline(v = lambda_under, col="red", lwd=3, lty=2)


# -----------------------------------------------------------------------------
# Simulate data
# -----------------------------------------------------------------------------
# inputs: 
#       df is a dataframe
#       w is a set of covariates' names 
#       a is a string, the name of treatment variable
#       y is a string, the name of outcome variable
#       type is a string, outcome type
#       g_fit is the undersmoothed HAL fit of g model
#       Q_fit is the undersmoothed HAL fit of Q model
#       rv is a scalar, the residual variance

get_simu_data <- function(df, w, a, y, g_fit, Q_fit, rv){
  simu_data <- dplyr::sample_n(df, size = nrow(df), replace = TRUE)
  
  # A Bootstrapping: undersmoothed HAL
  
  # generate predictions for A using g fit on original data and binomial sampling
  covariates <- simu_data %>% select(all_of(w)) %>% 
                          mutate_if(sapply(., is.factor), as.numeric)
  a_preds <- predict(g_fit, new_data = covariates)
  
  # Save bootstrapped intervention A to dataframe
  simu_data$a <- rbinom(length(a_preds), 1, prob = a_preds)
  
  # Y Bootstrapping: Super Learner
  # generate Y predictions using Q fit on original data and normal error
  covariates <- simu_data %>% select(all_of(c(w, a))) %>% 
                          mutate_if(sapply(., is.factor), as.numeric)
  y_preds <- predict(Q_fit, new_data = covariates)
  y_preds_error <- rnorm(length(y_preds), mean = 0, sd = sqrt(rv))
  y_preds <- y_preds + y_preds_error
  
  # Save bootstrapped outcome Y to dataframe
  simu_data$haz <- y_preds
  
  
  # Fit SL for g and Q on bootstrap data (for influence curves)
  
  #0630, drop constant cols
  drop_cov = 0
  for (covs in names(simu_data) %w/o% c('haz', 'a')) {
    if (length(unique(simu_data[[covs]])) == 1) {
      simu_data = simu_data %>% select(-all_of(covs))
      drop_cov = drop_cov + 1
    }
  }
  w =  colnames(simu_data) %w/o% c('haz', 'a')
  
  output <- list('data' = simu_data, 'drop_cov'=drop_cov)
  return(output)
}

# -----------------------------------------------------------------------------
# Simple substitution estimator for true ATE
# -----------------------------------------------------------------------------
ss_estimator <- function(df, w, a, y){
  
  set.seed(123)
  
  # Fit undersmoothed HAL for g on original data
  res_g <- under_HAL(df=df, yname = a, xname = w, Nlam = 100, type = "binomial")
  
  print(paste0(studies[i], " g lambda_init: ", res_g$lambda_init))
  print(paste0(studies[i], " g lambda_under: ", res_g$lambda_under))
  
  covariates <- df %>% select(all_of(w)) %>% 
                       mutate_if(sapply(., is.factor), as.numeric)
  n <- nrow(df)
  
  # tic()
  g_fit <- fit_hal(X = covariates,
                   Y = as.numeric(as.matrix(df %>% select(all_of(a)))),
                   family = "binomial",
                   num_knots = num_knots_generator(
                                   max_degree = ifelse(ncol(covariates) >= 20, 2, 3),
                                   smoothness_orders = 1,
                                   base_num_knots_0 = 500,
                                   base_num_knots_1 = 100
                                                       #ceiling(sqrt(n))
                                   ),
                   fit_control = list(
                     cv_select = FALSE,
                     n_folds = 10,
                     foldid = NULL,
                     use_min = TRUE,
                     lambda.min.ratio = 1e-4,
                     prediction_bounds = "default"
                   ),
                   lambda = res_g$lambda_under
                   )
  # toc()
  
  # Fit undersmoothed HAL for Q on original data
  res_Q <- under_HAL(df=df, yname = y, xname = c(w, a), Nlam = 100, type = "gaussian")
  
  print(paste0(studies[i], " Q lambda_init: ", res_Q$lambda_init))
  print(paste0(studies[i], " Q lambda_under: ", res_Q$lambda_under))
  
  covariates <- df %>% select(all_of(c(w, a))) %>% 
                       mutate_if(sapply(., is.factor), as.numeric)
  
  Q_fit <- fit_hal(X = covariates,
                   Y = as.numeric(as.matrix(df %>% select(all_of(y)))),
                   family = "gaussian",
                   num_knots = num_knots_generator(
                                   max_degree = ifelse(ncol(covariates) >= 20, 2, 3),
                                   smoothness_orders = 1,
                                   base_num_knots_0 = 500,
                                   base_num_knots_1 = 100
                                                       #ceiling(sqrt(n))
                                   ),
                   fit_control = list(
                     cv_select = FALSE,
                     n_folds = 10,
                     foldid = NULL,
                     use_min = TRUE,
                     lambda.min.ratio = 1e-4,
                     prediction_bounds = "default"
                   ),
                   lambda = res_Q$lambda_under
  )
  
  
  # Calculate true ATE by generating a large sample (n = 5*10^4)
  print("---calculate the truth---")
  
  list_ate <- rep(NA, 10)
  list_rss <- rep(NA, 10)
  
  # Divide and conquer
  tic()
  for (m in 1:10){
    print(paste0("calc ate: ", m, "/", 10))
      
    res_ate <- calc_ate(df = df, 
                        g_fit = g_fit, 
                        Q_fit = Q_fit, 
                        n_sample = 5*10^3)
      
    list_ate[m] <- res_ate$psi_ss
    list_rss[m] <- res_ate$rss
  }
  toc()
  
  # truth
  psi_ss <- mean(list_ate)
    
  # residual variance
  rv <- sum(list_rss)/(5*10^4)
  
  # Return Q and g and ATE
  ss <- c('psi' = psi_ss, 'g_fit' = list(g_fit), 'Q_fit' = list(Q_fit),
          'rv' = rv)
  return(ss)
}


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
# undersoomthed HAL helper function
get_maxscore <- function(basis_mat, resid_mat, sd_est, Nlam, us_fit){
  
  score_all <- matrix(NA, nrow = Nlam, 
                      ncol = ncol(basis_mat))
  
  for (i in 1:ncol(basis_mat)){
    u <- basis_mat[,i]
    score_mat <- resid_mat * u / sd_est[i]
    score_mean <- apply(score_mat, 2, mean)
    score_all[,i] <- score_mean
  }
  
  # absolute value
  max_score <- apply(abs(score_all), 1, max)
  return(max_score)
}

# helper function to calculate the true ATE
calc_ate <- function(df, g_fit, Q_fit, n_sample = 5*10^3){
  # sample W from emp
  large_data <- dplyr::sample_n(df, size = n_sample, replace = TRUE)
  
  # generate predictions for A using g HAL fit
  covariates <- large_data %>% select(all_of(w)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  a_preds <- predict(g_fit, new_data = covariates)
  
  # generate A
  large_data$a <- rbinom(length(a_preds), 1, prob = a_preds)
  
  covariates <- large_data %>% select(all_of(c(w, a))) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  # generate Y
  df0 <- data.frame(covariates)
  df1 <- data.frame(covariates)
  df0$a = 0
  df1$a = 1
  
  # QbarAW
  y_preds <- predict(Q_fit, new_data = covariates)
  
  # Qbar1W
  covariates1 <- df1 %>% select(all_of(c(w, a))) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  y1_preds <- predict(Q_fit, new_data = covariates1)
  
  # Qbar0W
  covariates0 <- df0 %>% select(all_of(c(w, a))) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  y0_preds <- predict(Q_fit, new_data = covariates0)
  
  # ate
  psi_ss <- mean(y1_preds - y0_preds)
  rss <- sum((y_preds - large_data %>% select(all_of(y)))^2)
  
  return(list("psi_ss" = psi_ss,
              "rss" = rss))
}


num_knots_generator <- function(max_degree, smoothness_orders, base_num_knots_0 = 500,
                                base_num_knots_1 = 200) {
  if (all(smoothness_orders > 0)) {
    return(sapply(seq_len(max_degree), function(d) {
      round(base_num_knots_1 / 2^(d - 1))
    }))
  }
  else {
    return(sapply(seq_len(max_degree), function(d) {
      round(base_num_knots_0 / 2^(d - 1))
    }))
  }
}

# tuning the num_knots (very slow)
tune_knots <- function(df, 
                       yname, 
                       xname, 
                       num_knots_candi = c(150, 200, 250, 300, 350),
                       type = "gaussian"){
  # variables
  y <- as.numeric(as.matrix(df %>% select(all_of(yname))))
  x <- df %>% select(all_of(xname)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  df_tune <- data.frame("num_knots" = NA,
                        "mse" = NA,
                        "n_coef" = NA,
                        "runtime" = NA)
  
  # tuning
  for(i in 1:length(num_knots_candi)){
    # fit
    start_time <- Sys.time()
    fit <- fit_hal(x,
                   y, 
                   return_x_basis = TRUE,
                   family = type,
                   num_knots = num_knots_generator(
                     max_degree = ifelse(ncol(x) >= 20, 2, 3),
                     smoothness_orders = 1,
                     base_num_knots_0 = 500,
                     base_num_knots_1 = num_knots_candi[i]
                   )
    )
    time_elapse <- as.numeric(Sys.time() - start_time, units="mins")
    
    # non-zero coef
    nonzero_col <- which(fit$coefs[-1] != 0)
    
    # performance
    preds <- predict(fit, new_data = x)
    
    df_tune[i, 1] <- num_knots_candi[i]
    df_tune[i, 2] <- mean((preds - y)^2)
    df_tune[i, 3] <- length(nonzero_col)
    df_tune[i, 4] <- time_elapse
  }
  return(df_tune)
}

# tic()
# temp_tune <- tune_knots(df=data_set, 
#                         yname='haz', 
#                         xname=c(w,a), 
#                         num_knots_candi = c(150, 200, 250, 300, 350),
#                         type = "gaussian")
# toc()
