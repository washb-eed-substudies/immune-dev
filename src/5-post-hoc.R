rm(list=ls())

source(here::here("0-config.R"))

d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-development-analysis-dataset.rds"))

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

pick_covariates<-function(j){
  if(j %in% c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2")){Wset=W2_immune.W2_dev}
  else{Wset=add_t3_covariates(j, W2_immune.W3_dev)}
  return(Wset)
}



# individual cytokines for th2/il10, th1 (IL-12 + IFN), th2 (IL-4 + IL-5 + IL-13) at year 1 
# outcomes y1 development, year 2 easq

Xvars <- c("t2_ln_il12", "t2_ln_il4", "t2_ln_il5", "t2_ln_il13", 
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

#Get primary contrasts
results <- NULL
for(i in 1:nrow(models)){
  res <- data.frame(X=models$X[i], Y=models$Y[i])
  preds <- predict_gam_diff(fit=models$fit[i][[1]], d=models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  results <-  bind_rows(results , preds$res)
}


Xvars <- c("t2_ratio_th1_th2")            
Yvars <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2", "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
           "z_cdi_say_t3", "z_cdi_und_t3") 

eff_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j, V="tr", W=pick_covariates(j))
    res <- data.frame(X=i, Y=j, V="tr", int.p =res_adj$int.p, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    eff_models <- bind_rows(eff_models, res)
  }
}

eff_results <- NULL
for(i in 1:nrow(eff_models)){
  preds <- predict_gam_emm(fit=eff_models$fit[i][[1]], d=eff_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=eff_models$X[i], Yvar=eff_models$Y[i])
  gamm_diff_res <- data.frame(V="tr", preds$res)
  eff_results <-  bind_rows(eff_results, preds$res)
}
