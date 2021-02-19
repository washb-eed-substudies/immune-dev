rm(list=ls())



source(here::here("0-config.R"))

# contains immune-growth covariates and exposures
d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))
names(d)

# merge in development outcomes
dev <- readRDS(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-development.RDS"))

dfull <- left_join(d, dev, "childid")



# check missingness
#Set list of adjustment variables
#Make vectors of adjustment variable names
Wvars<-c("sex","birthord", "momage","momheight","momedu", 
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "HHwealth",
         "ari7d_t2", "diar7d_t2", "nose7d_t2", 
         "fci_t2", "cesd_sum_t2", "life_viol_any_t3", "tr")

generate_miss_tbl <- function(Wvars, d){
  W <- d %>% select(all_of(Wvars))  
  miss <- data.frame(name = names(W), missing = colSums(is.na(W))/nrow(W), row.names = c(1:ncol(W)))
  for (i in 1:nrow(miss)) {
    miss$class[i] <- class(W[,which(colnames(W) == miss[i, 1])])
  }
  miss 
}

generate_miss_tbl(Wvars, dfull)

Wvars22<-c("ageday_bt2", "agedays_motor", "month_bt2", "month_motor", "laz_t1", "waz_t1") 
Wvars33<-c("ageday_bt3", "month_bt3",
           "ari7d_t3", "diar7d_t3", "nose7d_t3", "fci_t3", "laz_t2", "waz_t2",
           "cesd_sum_ee_t3", "pss_sum_mom_t3") 
Wvars23<-c("ageday_bt2", "month_bt2",
           "ari7d_t3", "diar7d_t3", "nose7d_t3", "fci_t3", "laz_t2", "waz_t2",
           "cesd_sum_ee_t3", "pss_sum_mom_t3")

W2_immune.W2_dev <- c(Wvars, Wvars22) %>% unique(.)
W3_immune.W3_dev <- c(Wvars, Wvars33) %>% unique(.)
W2_immune.W3_dev <- c(Wvars, Wvars23) %>% unique(.)

add_t3_covariates <- function(j, W){
  if(grepl("easq", j)){return (c(W, "agedays_easq", "month_easq"))}
  else if(grepl("cdi", j) & grepl("t3", j)){return (c(W, "agedays_cdi_t3", "month_cdi_t3"))}
}

Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "t2_ln_ifn", "sumscore_t2_Z", "t2_ln_igf")            
Yvars <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2") 

for (i in Xvars){
  for (j in Yvars){
    print(i)
    print(j)
    W = W2_immune.W2_dev
    d_sub <- subset(dfull, !is.na(dfull[,i]) & !is.na(dfull[,j]))
    print(generate_miss_tbl(W, d_sub))
  }
}

# all exposures at y2 and growth outcomes at y2
Xvars <- c("t3_ratio_pro_il10", "t3_ratio_il2_il10", "t3_ratio_gmc_il10", "t3_ratio_th1_il10", "t3_ratio_th2_il10",     
           "t3_ratio_th17_il10", "t3_ratio_th1_th2", "t3_ratio_th1_th17", "t3_ln_ifn", "sumscore_t3_Z", "t3_ln_igf")            
Yvars <- c("z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
           "z_cdi_say_t3", "z_cdi_und_t3") 

for (i in Xvars){
  for (j in Yvars){
    print(i)
    print(j)
    W = add_t3_covariates(j, W3_immune.W3_dev)
    d_sub <- subset(dfull, !is.na(dfull[,i]) & !is.na(dfull[,j]))
    print(generate_miss_tbl(W, d_sub))
  }
}

# all exposures at y2 and growth outcomes at y2
Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "t2_ln_ifn", "sumscore_t2_Z", "t2_ln_igf")            
Yvars <- c("z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
           "z_cdi_say_t3", "z_cdi_und_t3") 

for (i in Xvars){
  for (j in Yvars){
    print(i)
    print(j)
    W = add_t3_covariates(j, W2_immune.W3_dev)
    d_sub <- subset(dfull, !is.na(dfull[,i]) & !is.na(dfull[,j]))
    print(generate_miss_tbl(W, d_sub))
  }
}

saveRDS(dfull, paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-development-analysis-dataset.rds"))
