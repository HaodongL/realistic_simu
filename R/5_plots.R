# -----------------------------------------------------------------------------
# Plots and figures
# -----------------------------------------------------------------------------
# April-2021 new plots
# Figure 1, major performance metrics

rm(list=ls())
table_results_data <- read.csv(paste0('~/Repo/...', "HAL_500res_07212021", '.csv'))

library(ggplot2)
library(ggpubr)

theme_nice <- theme(panel.background = element_blank(), 
                    panel.grid = element_line(size = 0.3, linetype = 'dashed',
                                              colour = "grey"))

temp_res <- na.omit(table_results_data)
temp_tbl_p1 <- gather(temp_res, "method", "coverage", tmle_coverage:cvaiptw_coverage) %>% select(coverage)
temp_tbl_p2 <- gather(temp_res, "method", "bias", tmle_bias:cvaiptw_bias) %>% select(bias)
temp_tbl_p3 <- gather(temp_res, "method", "meanCIwidth", tmle_meanwidthCI:cvaiptw_meanwidthCI) %>% select(meanCIwidth)
temp_tbl_p4 <- gather(temp_res, "method", "rmse", tmle_rmse:cvaiptw_rmse) %>% select(rmse)


temp_tbl_long <- data.frame(studyid = rep(temp_res$studyid[1:10],7),
                            true_psi = rep(temp_res$ss[1:10],7),
                            method = rep(c("TMLE", "CVTMLE", "CTMLE", 
                                           "IPTW", "CVIPTW", "AIPTW", "CVAIPTW"), each = 10))

temp_tbl_long$studyid <- as.factor(rep(c(7,4,6,9,10,8,1,5,2,3),7))

temp_tbl_long <- cbind(temp_tbl_long, 
                       temp_tbl_p1,
                       temp_tbl_p2,
                       temp_tbl_p3,
                       temp_tbl_p4)

temp_tbl_long <- temp_tbl_long %>% 
  group_by(method) %>% 
  mutate(ave_rmse = mean(rmse),
         ave_coverage = mean(coverage),
         ave_bias = mean(bias),
         ave_meanCIwidth = mean(meanCIwidth),
         med_rmse = median(rmse),
         med_coverage = median(coverage),
         med_bias = median(bias),
         med_meanCIwidth = median(meanCIwidth))

temp_tbl_long <- temp_tbl_long %>% arrange(method, studyid)

library(RColorBrewer)
color.pal <- brewer.pal(n = 10, name = "Paired")

# rMSE
plot1 <- 
  temp_tbl_long %>% 
  mutate(rmse = ifelse(rmse >= 14, 6.3, rmse)) %>% 
  ggplot() + 
  geom_point(aes(x= rmse, y =method, color = studyid, shape= method), size= 4) +
  scale_shape_manual(values=rep(7,7), guide=FALSE) + 
  scale_color_manual(values=color.pal, name = "StudyID") + 
  geom_point(aes(x= med_rmse, y = method), 
             shape = 21, color = "black", fill = "black", size= 2)+
  coord_cartesian(xlim = c(0,6)) +
  geom_vline(aes(xintercept = 1), colour="#000000",size=0.5)+
  scale_y_discrete(limits=rev) + 
  labs(y=NULL,
       x="rMSE") +
  theme_nice

# Bias
plot2 <- 
  temp_tbl_long %>% 
  ggplot() + 
  geom_point(aes(x= bias, y =method, color = studyid, shape= method), size= 4) +
  scale_shape_manual(values=rep(7,7), guide=FALSE) + 
  scale_color_manual(values=color.pal, name = "StudyID") + 
  geom_point(aes(x= med_bias, y = method), 
             shape = 21, color = "black", fill = "black", size= 2) +
  coord_cartesian(xlim = c(-0.35,0.1)) + 
  geom_vline(aes(xintercept = 0), colour="#000000",size=0.5) +
  scale_y_discrete(limits=rev) + 
  labs(y=NULL,
       x="Bias") +
  theme_nice

# Coverage
plot3 <- 
  temp_tbl_long %>% 
  ggplot() + 
  geom_point(aes(x= coverage, y =method, color = studyid, shape= method), size= 4)+
  scale_shape_manual(values=rep(7,7), guide=FALSE) + 
  scale_color_manual(values=color.pal, name = "StudyID") + 
  geom_point(aes(x= med_coverage, y = method), 
             shape = 21, color = "black", fill = "black", size= 2)+
  # coord_cartesian(xlim = c(0,1)) + 
  geom_vline(aes(xintercept = 0.95), colour="#000000",size=0.5) +
  scale_y_discrete(limits=rev) + 
  labs(y=NULL,
       x="Coverage of 95% CI",
       legen) +
  theme_nice

# meanCIwidth
plot4 <- 
  temp_tbl_long %>% 
  ggplot() + 
  geom_point(aes(x= meanCIwidth, y =method, color = studyid, shape= method), size= 4) +
  scale_shape_manual(values=rep(7,7), guide=FALSE) + 
  scale_color_manual(values=color.pal, name = "StudyID") + 
  geom_point(aes(x= med_meanCIwidth, y = method), 
             shape = 21, color = "black", fill = "black", size= 2)+
  scale_y_discrete(limits=rev) + 
  theme_nice + 
  labs(y=NULL,
       x="Average width of 95% CI") +
  theme(legend.key.width=unit(1,"cm"))

plot_dots <-
  ggarrange(plot1,
            plot2,
            plot3,
            plot4,
            legend = "top", common.legend = TRUE)

plotpath = paste0("~/Repo/...")
ggsave('plot_dots.pdf', width = 14, height = 8,
       plot = plot_dots, path = plotpath)


