rm(list=ls())

source(here::here("0-config.R"))

d <- d %>%
  filter(tr %in% c("Nutrition + WSH", "Control"))

# individual cytokines for th2/il10, th1 (IL-12 + IFN), th2 (IL-4 + IL-5 + IL-13) at year 1 
# outcomes y1 development, year 2 easq

Xvars <- c("t2_ln_il12", "t2_ln_il4", "t2_ln_il5", "t2_ln_il13", 
           "t2_ratio_il12_il10", "t2_ratio_ifn_il10", 
           "t2_ratio_il4_il10", "t2_ratio_il5_il10", "t2_ratio_il13_il10")            
Yvars <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2", "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
           "z_cdi_say_t3", "z_cdi_und_t3") 

#Fit models
models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    models <- bind_rows(models, res)
  }
}

Xvars <- c("t3_ln_il12", "t3_ln_il4", "t3_ln_il5", "t3_ln_il13", 
           "t3_ratio_il12_il10", "t3_ratio_ifn_il10", 
           "t3_ratio_il4_il10", "t3_ratio_il5_il10", "t3_ratio_il13_il10")            
Yvars <- c("z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
           "z_cdi_say_t3", "z_cdi_und_t3") 

#Fit models
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    models <- bind_rows(models, res)
  }
}

#Get primary contrasts
results <- NULL
for(i in 1:nrow(models)){
  res <- data.frame(X=models$X[i], Y=models$Y[i])
  preds <- predict_gam_diff(fit=models$fit[i][[1]], d=models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  results <-  bind_rows(results , preds$res)
}

#Make list of plots
plot_list <- NULL
plot_data <- NULL
for(i in 1:nrow(models)){
  res <- data.frame(X=models$X[i], Y=models$Y[i])
  simul_plot <- gam_simul_CI(models$fit[i][[1]], models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  plot_list[[i]] <-  simul_plot$p
  plot_data <-  rbind(plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}

#Save models
saveRDS(models, here("models/individual_unadj_models.RDS"))

#Save results
saveRDS(results, here("results/unadjusted/individual_unadj_res.RDS"))

saveRDS(plot_data, here("figure-data/individual_unadj_spline_data.RDS"))



#Set list of adjustment variables
#Make vectors of adjustment variable names
Wvars<-c("sex","birthord", "momage","momheight","momedu", 
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "HHwealth_scaled", "tr",
         "ari7d_t2", "diar7d_t2", "nose7d_t2", 
         "fci_t2", "cesd_sum_t2", "life_viol_any_t3")

Wvars[!(Wvars %in% colnames(d))]

#Add in time varying covariates:

#NOTES
#Does monsoon_ut2 need to be replaced with monsoon_ht2 for growth measures? (and agemth_ut2 with agedays_ht2?)
Wvars22<-c("ageday_bt2", "agedays_motor", "month_bt2", "month_motor", "laz_t1_cat", "waz_t1_cat") 
Wvars33<-c("ageday_bt3", "month_bt3",
           "ari7d_t3", "diar7d_t3", "nose7d_t3", "fci_t3", "laz_t2_cat", "waz_t2_cat",
           "cesd_sum_ee_t3", "pss_sum_mom_t3") 
Wvars23<-c("ageday_bt2", "month_bt2",
           "ari7d_t3", "diar7d_t3", "nose7d_t3", "fci_t3", "laz_t2", "waz_t2",
           "cesd_sum_ee_t3", "pss_sum_mom_t3")

W2_immune.W2_dev <- c(Wvars, Wvars22) %>% unique(.)
W3_immune.W3_dev <- c(Wvars, Wvars33) %>% unique(.)
W2_immune.W3_dev <- c(Wvars, Wvars23) %>% unique(.)

W2_immune.W2_dev[!(W2_immune.W2_dev %in% colnames(d))]
W3_immune.W3_dev[!(W3_immune.W3_dev %in% colnames(d))]
W2_immune.W3_dev[!(W2_immune.W3_dev %in% colnames(d))]

add_t3_covariates <- function(j, W){
  if(grepl("easq", j)){return (c(W, "agedays_easq", "month_easq"))}
  else if(grepl("cdi", j) & grepl("t3", j)){return (c(W, "agedays_cdi_t3", "month_cdi_t3"))}
}

pick_covariates<-function(j){
  if(j %in% c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2")){Wset=W2_immune.W2_dev}
  else{Wset=add_t3_covariates(j, W2_immune.W3_dev)}
  return(Wset)
}


# individual cytokines for th2/il10, th1 (IL-12 + IFN), th2 (IL-4 + IL-5 + IL-13) at year 1 
# outcomes y1 development, year 2 easq

Xvars <- c("t2_ln_il12", "t2_ln_il4", "t2_ln_il5", "t2_ln_il13", 
           "t2_ratio_il12_il10", "t2_ratio_ifn_il10", 
           "t2_ratio_il4_il10", "t2_ratio_il5_il10", "t2_ratio_il13_il10")            
Yvars <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2", "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
           "z_cdi_say_t3", "z_cdi_und_t3") 

#Fit models
models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=pick_covariates(j))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    models <- bind_rows(models, res)
  }
}

Xvars <- c("t3_ln_il12", "t3_ln_il4", "t3_ln_il5", "t3_ln_il13", 
           "t3_ratio_il12_il10", "t3_ratio_ifn_il10", 
           "t3_ratio_il4_il10", "t3_ratio_il5_il10", "t3_ratio_il13_il10")            
Yvars <- c("z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
           "z_cdi_say_t3", "z_cdi_und_t3") 

#Fit models
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=W3_immune.W3_dev)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    models <- bind_rows(models, res)
  }
}

#Get primary contrasts
results <- NULL
for(i in 1:nrow(models)){
  res <- data.frame(X=models$X[i], Y=models$Y[i])
  preds <- predict_gam_diff(fit=models$fit[i][[1]], d=models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  results <-  bind_rows(results , preds$res)
}

#Make list of plots
plot_list <- NULL
plot_data <- NULL
for(i in 1:nrow(models)){
  res <- data.frame(X=models$X[i], Y=models$Y[i])
  simul_plot <- gam_simul_CI(models$fit[i][[1]], models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  plot_list[[i]] <-  simul_plot$p
  plot_data <-  rbind(plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}

#Save models
saveRDS(models, here("models/individual_adj_models.RDS"))

#Save results
saveRDS(results, here("results/adjusted/individual_adj_res.RDS"))

saveRDS(plot_data, here("figure-data/individual_adj_spline_data.RDS"))




#### hazard ratio for WHO motor milestones
Xvars <- c("t2_ln_il12", "t2_ln_il4", "t2_ln_il5", "t2_ln_il13", 
           "t2_ratio_il12_il10", "t2_ratio_ifn_il10", 
           "t2_ratio_il4_il10", "t2_ratio_il5_il10", "t2_ratio_il13_il10")            
Yvars <- grep("who_", colnames(d), value=T)

HR_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_HR_GAM(d=d, X=i, Y=j, age="agedays_motor", W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    HR_models <- bind_rows(HR_models, res)
  }
}

HR_res <- NULL
for(i in 1:nrow(HR_models)){
  res <- data.frame(X=HR_models$X[i], Y=HR_models$Y[i])
  preds <- predict_gam_HR(fit=HR_models$fit[i][[1]], d=HR_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  HR_res <-  bind_rows(HR_res , preds$res)
}

#Make list of plots
HR_plot_list <- NULL
HR_plot_data <- NULL
for(i in 1:nrow(HR_models)){
  res <- data.frame(X=HR_models$X[i], Y=HR_models$Y[i])
  simul_plot <- gam_simul_CI(HR_models$fit[i][[1]], HR_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  HR_plot_list[[i]] <-  simul_plot$p
  HR_plot_data <-  rbind(HR_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


#Save models
saveRDS(HR_models, here("models/HR_ind_models.RDS"))

#Save results
saveRDS(HR_res, here("results/unadjusted/HR_ind_res.RDS"))

#Save plot data
saveRDS(HR_plot_data, here("figure-data/HR_ind_unadj_spline_data.RDS"))



Xvars <- c("t2_ln_il12", "t2_ln_il4", "t2_ln_il5", "t2_ln_il13", 
           "t2_ratio_il12_il10", "t2_ratio_ifn_il10", 
           "t2_ratio_il4_il10", "t2_ratio_il5_il10", "t2_ratio_il13_il10")            
Yvars <- grep("who_", colnames(d), value=T)

HR_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_HR_GAM(d=d, X=i, Y=j, age="agedays_motor", W=W2_immune.W2_dev)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    HR_models <- bind_rows(HR_models, res)
  }
}

HR_res <- NULL
for(i in 1:nrow(HR_models)){
  res <- data.frame(X=HR_models$X[i], Y=HR_models$Y[i])
  preds <- predict_gam_HR(fit=HR_models$fit[i][[1]], d=HR_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  HR_res <-  bind_rows(HR_res , preds$res)
}

#Make list of plots
HR_plot_list <- NULL
HR_plot_data <- NULL
for(i in 1:nrow(HR_models)){
  res <- data.frame(X=HR_models$X[i], Y=HR_models$Y[i])
  simul_plot <- gam_simul_CI(HR_models$fit[i][[1]], HR_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  HR_plot_list[[i]] <-  simul_plot$p
  HR_plot_data <-  rbind(HR_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


#Save models
saveRDS(HR_models, here("models/HR_ind_adj_models.RDS"))

#Save results
saveRDS(HR_res, here("results/adjusted/HR_ind_res.RDS"))

#Save plot data
saveRDS(HR_plot_data, here("figure-data/HR_ind_adj_spline_data.RDS"))




Xvars <- c("t2_ratio_th1_th2")
Yvars <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2", "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq",
           "z_cdi_say_t3", "z_cdi_und_t3")

# eff_models <- NULL
# for(i in Xvars){
#   for(j in Yvars){
#     print(i)
#     print(j)
#     W=pick_covariates(j)
#     res_adj <- fit_RE_gam(d=d, X=i, Y=j, V="tr", W=W[W!="tr"])
#     res <- data.frame(X=i, Y=j, V="tr", int.p =res_adj$int.p, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
#     eff_models <- bind_rows(eff_models, res)
#   }
# }
# 
# eff_results <- NULL
# for(i in 1:nrow(eff_models)){
#   preds <- predict_gam_emm(fit=eff_models$fit[i][[1]], d=eff_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=eff_models$X[i], Yvar=eff_models$Y[i])
#   gamm_diff_res <- data.frame(V="tr", preds$res) %>% mutate(int.Pval = c(NA, eff_models$int.p[[i]]))
#   eff_results <-  bind_rows(eff_results, gamm_diff_res)
# }
# 
# view(eff_results)
