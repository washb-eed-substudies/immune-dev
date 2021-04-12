rm(list=ls())

source(here::here("0-config.R"))

dfull <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-development-analysis-dataset.rds"))
d <- dfull %>% filter(tr %in% c("Nutrition + WSH", "Control")) 
  
#Set list of adjustment variables
#Make vectors of adjustment variable names
Wvars<-c("sex","birthord", "momage","momheight","momedu", 
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "HHwealth", "tr",
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

#Loop over exposure-outcome pairs

#### Hypothesis 1: immune status associated with concurrent child development ####
# all immune ratios at Y1 v. development outcomes at Y1
Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "t2_ln_ifn")            
Yvars <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2") 

#Fit models
H1_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=W2_immune.W2_dev)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_models <- bind_rows(H1_adj_models, res)
  }
}

# all immune outcomes at y2 and growth outcomes at y2
Xvars <- c("t3_ratio_pro_il10", "t3_ratio_il2_il10", "t3_ratio_gmc_il10", "t3_ratio_th1_il10", "t3_ratio_th2_il10",     
           "t3_ratio_th17_il10", "t3_ratio_th1_th2", "t3_ratio_th1_th17", "t3_ln_agp", "t3_ln_crp", "t3_ln_ifn")            
Yvars <- c("z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
           "z_cdi_say_t3", "z_cdi_und_t3") 

for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    if(i %in% c("t3_ln_agp", "t3_ln_crp")){
      dfunc <- dfull
    }else{dfunc <- d}
    res_adj <- fit_RE_gam(d=dfunc, X=i, Y=j,  W=add_t3_covariates(j, W3_immune.W3_dev))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_models <- bind_rows(H1_adj_models, res)
  }
}

#Get primary contrasts
H1_adj_res <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H1_adj_models$fit[i][[1]], d=H1_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H1_adj_res <-  bind_rows(H1_adj_res , preds$res)
}

#Make list of plots
H1_adj_plot_list <- NULL
H1_adj_plot_data <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H1_adj_models$fit[i][[1]], H1_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_adj_plot_list[[i]] <-  simul_plot$p
  H1_adj_plot_data <-  rbind(H1_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred %>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
saveRDS(H1_adj_models, here("models/H1_adj_models.RDS"))

#Save results
saveRDS(H1_adj_res, here("results/adjusted/H1_adj_res.RDS"))


#Save plots
#saveRDS(H1_adj_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H1_adj_splines.RDS"))

#Save plot data
saveRDS(H1_adj_plot_data, here("figure-data/H1_adj_spline_data.RDS"))



#### Hypothesis 2: immune status and subsequent development ####
# all immune outcomes at y1 v. development at y2
Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "t2_ln_ifn")            
Yvars <- c("z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
           "z_cdi_say_t3", "z_cdi_und_t3") 

#Fit models
H2_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=add_t3_covariates(j, W2_immune.W3_dev))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H2_adj_models <- bind_rows(H2_adj_models, res)
  }
}

#Get primary contrasts
H2_adj_res <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H2_adj_models$fit[i][[1]], d=H2_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H2_adj_res <-  bind_rows(H2_adj_res , preds$res)
}

#Make list of plots
H2_plot_list <- NULL
H2_plot_data <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H2_adj_models$fit[i][[1]], H2_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H2_plot_list[[i]] <-  simul_plot$p
  H2_plot_data <-  rbind(H2_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
saveRDS(H2_adj_models, here("models/H2_adj_models.RDS"))

#Save results
saveRDS(H2_adj_res, here("results/adjusted/H2_adj_res.RDS"))


#Save plots
#saveRDS(H2_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H2_adj_splines.RDS"))

#Save plot data
saveRDS(H2_plot_data, here("figure-data/H2_adj_spline_data.RDS"))



#### Hypothesis 3: sum score and development ####
# sum score and concurrent and subsequent development
Xvars <- c("sumscore_t2_Z")            
Yvars <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2",
           "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
           "z_cdi_say_t3", "z_cdi_und_t3") 

pick_covariates<-function(j){
  if(j %in% c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2")){Wset=W2_immune.W2_dev}
  else{Wset=add_t3_covariates(j, W2_immune.W3_dev)}
  return(Wset)
}

#Fit models
H3_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset=pick_covariates(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H3_adj_models <- bind_rows(H3_adj_models, res)
  }
}

Xvars <- c("sumscore_t3_Z")            
Yvars <- c("z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
           "z_cdi_say_t3", "z_cdi_und_t3")

for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=add_t3_covariates(j, W3_immune.W3_dev))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H3_adj_models <- bind_rows(H3_adj_models, res)
  }
}

#Get primary contrasts
H3_adj_res <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H3_adj_models$fit[i][[1]], d=H3_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H3_adj_res <-  bind_rows(H3_adj_res , preds$res)
}

#Make list of plots
H3_plot_list <- NULL
H3_plot_data <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H3_adj_models$fit[i][[1]], H3_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H3_plot_list[[i]] <-  simul_plot$p
  H3_plot_data <-  rbind(H3_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
saveRDS(H3_adj_models, here("models/H3_adj_models.RDS"))

#Save results
saveRDS(H3_adj_res, here("results/adjusted/H3_adj_res.RDS"))


#Save plots
#saveRDS(H3_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H3_adj_splines.RDS"))

#Save plot data
saveRDS(H3_plot_data, here("figure-data/H3_adj_spline_data.RDS"))


#### Hypothesis 4 ####
# IGF and concurrent and subsequent development
Xvars <- c("t2_ln_igf")            
Yvars <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2",
           "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
           "z_cdi_say_t3", "z_cdi_und_t3")

#Fit models
H4_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset<-pick_covariates(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H4_adj_models <- bind_rows(H4_adj_models, res)
  }
}

Xvars <- c("t3_ln_igf")            
Yvars <- c("z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
           "z_cdi_say_t3", "z_cdi_und_t3") # add "z_cdi_und_t2", "z_cdi_say_t2" later

for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=add_t3_covariates(j, W3_immune.W3_dev))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H4_adj_models <- bind_rows(H4_adj_models, res)
  }
}

#Get primary contrasts
H4_adj_res <- NULL
for(i in 1:nrow(H4_adj_models)){
  res <- data.frame(X=H4_adj_models$X[i], Y=H4_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H4_adj_models$fit[i][[1]], d=H4_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H4_adj_res <-  bind_rows(H4_adj_res , preds$res)
}

#Make list of plots
H4_adj_plot_list <- NULL
H4_adj_plot_data <- NULL
for(i in 1:nrow(H4_adj_models)){
  res <- data.frame(X=H4_adj_models$X[i], Y=H4_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H4_adj_models$fit[i][[1]], H4_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H4_adj_plot_list[[i]] <-  simul_plot$p
  H4_adj_plot_data <-  rbind(H4_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
saveRDS(H4_adj_models, here("models/H4_adj_models.RDS"))

#Save results
saveRDS(H4_adj_res, here("results/adjusted/H4_adj_res.RDS"))


#Save plots
#saveRDS(H4_adj_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H4_adj_adj_splines.RDS"))

#Save plot data
saveRDS(H4_adj_plot_data, here("figure-data/H4_adj_spline_data.RDS"))


#### hazard ratio for WHO motor milestones

Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "t2_ln_ifn", "sumscore_t2_Z",
           "t2_ln_igf") 
Yvars <- grep("who_", colnames(d), value=T)

HR_models <- NULL
for(i in Xvars){
  for(j in Yvars[4:6]){
    print(i)
    print(j)
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
HR_adj_plot_list <- NULL
HR_adj_plot_data <- NULL
for(i in 1:nrow(HR_models)){
  res <- data.frame(X=HR_models$X[i], Y=HR_models$Y[i])
  simul_plot <- gam_simul_CI(HR_models$fit[i][[1]], HR_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  HR_adj_plot_list[[i]] <-  simul_plot$p
  HR_adj_plot_data <-  rbind(HR_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
saveRDS(HR_models, here("models/HR_adj_models.RDS"))

#Save results
saveRDS(HR_res, here("results/adjusted/HR_adj_res.RDS"))

#Save plot data
saveRDS(HR_adj_plot_data, here("figure-data/HR_adj_spline_data.RDS"))

