
rm(list=ls())

library(dplyr)
library(tidyr)




result <- dplyr::bind_rows( .id = 'studyid')

result = result %>%
  dplyr::mutate(gcomp_upper= gcomp + 1.96*tmle_se,
                gcomp_lower= gcomp - 1.96*tmle_se,)

table_results_data <- result %>%
  dplyr::mutate(tmle_proportion = ss <= tmle_upper & ss >= tmle_lower,
                cvtmle_proportion = ss <= cvtmle_upper & ss >= cvtmle_lower,
                ctmle_proportion = ss <= ctmle_upper & ss >= ctmle_lower,
                iptw_proportion = ss <= iptw_upper & ss >= iptw_lower,
                cviptw_proportion = ss <= cviptw_upper & ss >= cviptw_lower,
                aiptw_proportion = ss <= aiptw_upper & ss >= aiptw_lower,
                cvaiptw_proportion = ss <= cvaiptw_upper & ss >= cvaiptw_lower,
                gcomp_proportion = ss <= gcomp_upper & ss >= gcomp_lower,
                
                
                tmle_widthCI = tmle_upper-tmle_lower,
                cvtmle_widthCI = cvtmle_upper-cvtmle_lower,
                ctmle_widthCI = ctmle_upper-ctmle_lower,
                iptw_widthCI = iptw_upper-iptw_lower,
                cviptw_widthCI = cviptw_upper-cviptw_lower,
                aiptw_widthCI = aiptw_upper-aiptw_lower,
                cvaiptw_widthCI = cvaiptw_upper-cvaiptw_lower,
                gcomp_widthCI = gcomp_upper-gcomp_lower
  ) %>%
  dplyr::group_by(studyid, ss) %>%
  summarize(tmle_coverage = mean(tmle_proportion),
            cvtmle_coverage = mean(cvtmle_proportion),
            ctmle_coverage = mean(ctmle_proportion),
            iptw_coverage = mean(iptw_proportion),
            cviptw_coverage = mean(cviptw_proportion),
            aiptw_coverage = mean(aiptw_proportion),
            cvaiptw_coverage = mean(cvaiptw_proportion),
            gcomp_coverage = mean(gcomp_proportion),
            
            tmle_bias = mean(tmle) - mean(ss),
            cvtmle_bias = mean(cvtmle) - mean(ss),
            ctmle_bias = mean(ctmle) - mean(ss),
            iptw_bias = mean(iptw) - mean(ss),
            cviptw_bias = mean(cviptw) - mean(ss),
            aiptw_bias = mean(aiptw) - mean(ss),
            cvaiptw_bias = mean(cvaiptw) - mean(ss),
            gcomp_bias = mean(gcomp) - mean(ss),
            
            
            tmle_var = var(tmle),
            cvtmle_var = var(cvtmle),
            ctmle_var = var(ctmle),
            iptw_var = var(iptw),
            cviptw_var = var(cviptw),
            aiptw_var = var(aiptw),
            cvaiptw_var = var(cvaiptw),
            gcomp_var = var(gcomp),
            
            tmle_mse = tmle_bias^2 + var(tmle),
            cvtmle_mse = cvtmle_bias^2 + var(cvtmle),
            ctmle_mse = ctmle_bias^2 + var(ctmle),
            iptw_mse = iptw_bias^2 + var(iptw),
            cviptw_mse = cviptw_bias^2 + var(cviptw),
            aiptw_mse = aiptw_bias^2 + var(aiptw),
            cvaiptw_mse = cvaiptw_bias^2 + var(cvaiptw),
            gcomp_mse = gcomp_bias^2 + var(gcomp),
            
            # April-2021 rMSE
            tmle_rmse = tmle_mse/iptw_mse,
            cvtmle_rmse = cvtmle_mse/iptw_mse,
            ctmle_rmse = ctmle_mse/iptw_mse,
            iptw_rmse = iptw_mse/iptw_mse,
            cviptw_rmse = cviptw_mse/iptw_mse,
            aiptw_rmse = aiptw_mse/iptw_mse,
            cvaiptw_rmse = cvaiptw_mse/iptw_mse,
            gcomp_rmse = gcomp_mse/iptw_mse,
            
            # 2020-02-01 coverage of oracle CI
            tmle_oracle = mean(ss <= tmle + 1.96*sd(tmle) & ss >= tmle - 1.96*sd(tmle)),
            cvtmle_oracle = mean(ss <= cvtmle + 1.96*sd(cvtmle) & ss >= tmle - 1.96*sd(cvtmle)),
            ctmle_oracle = mean(ss <= ctmle + 1.96*sd(ctmle) & ss >= ctmle - 1.96*sd(ctmle)),
            iptw_oracle = mean(ss <= iptw + 1.96*sd(iptw) & ss >= iptw - 1.96*sd(iptw)),
            cviptw_oracle = mean(ss <= cviptw + 1.96*sd(cviptw) & ss >= cviptw - 1.96*sd(cviptw)),
            aiptw_oracle = mean(ss <= aiptw + 1.96*sd(aiptw) & ss >= aiptw - 1.96*sd(aiptw)),
            cvaiptw_oracle = mean(ss <= cvaiptw + 1.96*sd(cvaiptw) & ss >= cvaiptw - 1.96*sd(cvaiptw)),
            
            gcomp_oracle = mean(ss <= gcomp + 1.96*sd(gcomp) & ss >= gcomp - 1.96*sd(gcomp)),

            
            tmle_meanwidthCI = mean(tmle_widthCI),
            cvtmle_meanwidthCI = mean(cvtmle_widthCI),
            ctmle_meanwidthCI = mean(ctmle_widthCI),
            iptw_meanwidthCI = mean(iptw_widthCI),
            cviptw_meanwidthCI = mean(cviptw_widthCI),
            aiptw_meanwidthCI = mean(aiptw_widthCI),
            cvaiptw_meanwidthCI = mean(cvaiptw_widthCI),
            gcomp_meanwidthCI = mean(gcomp_widthCI),
            
            drop_cov = 1- mean(as.numeric(drop_cov==0))
  )


