rm(list=ls())

library('flextable')
library('officer')
source(here::here("0-config.R"))

# load enrollment characteristics and results
# d <- read.csv(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-stress-growth-covariates-stresslab-anthro.csv"))
H1 <- readRDS(here('results/unadjusted/H1_res.RDS'))
H2 <- readRDS(here('results/unadjusted/H2_res.RDS'))
H3 <- readRDS(here('results/unadjusted/H3_res.RDS'))
H4 <- readRDS(here('results/unadjusted/H4_res.RDS'))
H1adj <- readRDS(here('results/adjusted/H1_adj_res.RDS'))
H2adj <- readRDS(here('results/adjusted/H2_adj_res.RDS'))
H3adj <- readRDS(here('results/adjusted/H3_adj_res.RDS'))
H4adj <- readRDS(here('results/adjusted/H4_adj_res.RDS'))

unadj <- rbind(H1, H2, H3, H4)
adj <- rbind(H1adj, H2adj, H3adj, H4adj)


#### Functions for growth tables ####
growth_tbl <- function(name, expo_var, out_var, exposure, outcome, results, results_adj){
  ### name: string name of group of exposures
  ### expo_var: vector of string exposures to include in table
  ### out_var: vector of string outcomes to include in table
  ### exposure: vector of string exposure variable names
  ### outcome: vector of string outcome variable names
  ### results: data frame with unadjusted results
  ### results_adj: data fram with adjusted results
  
  ### this function produces a table that can be saved as a csv
  
  tbl <- data.table(name = character(), "Outcome" = character(), "N" = character(), "25th Percentile" = character(), "75th Percentile" = character(),
                    " Outcome, 75th Percentile v. 25th Percentile" = character(), " " = character(), " " = character(), " " = character(), 
                    " " = character(), " " = character(), " " = character(), " " = character())
  tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", "Unadjusted", " ", " ", " ", "Fully adjusted", " ", " ", " "))
  tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", 
                         "Predicted Outcome at 25th Percentile", "Predicted Outcome at 75th Percentile", "Coefficient (95% CI)", "P-value", 
                         "Predicted Outcome at 25th Percentile", "Predicted Outcome at 75th Percentile", "Coefficient (95% CI)", "P-value"))
  for (i in 1:length(exposure)) {
    for (j in 1:length(outcome)) {
      exp <- exposure[i]
      out <- outcome[j]
      filtered_res <- results[results$Y==out & results$X==exp,]
      filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
      unadj <- paste(round(filtered_res$`point.diff`, 2), " (", round(filtered_res$`lb.diff`, 2), ", ", round(filtered_res$`ub.diff`, 2), ")", sep="")
      adj <- paste(round(filtered_adj$`point.diff`, 2), " (", round(filtered_adj$`lb.diff`, 2), ", ", round(filtered_adj$`ub.diff`, 2), ")", sep="")
      if (j==1){
        tbl <- rbind(tbl, list(expo_var[i], out_var[j], filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$corrected.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$corrected.Pval, 2)))
      }else {
        tbl <- rbind(tbl, list("", out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$corrected.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$corrected.Pval, 2)))
      }
    }
    if (i != length(exposure)) {
      tbl <- rbind(tbl, list("","","","","","","","","","","","",""))
    }
  }
  tbl
}

growth_tbl_flex <- function(name, expo_var, out_var, exposure, outcome, results, results_adj){
  ### name: string name of group of exposures
  ### expo_var: vector of string exposures to include in table
  ### out_var: vector of string outcomes to include in table
  ### exposure: vector of string exposure variable names
  ### outcome: vector of string outcome variable names
  ### results: data frame with unadjusted results
  ### results_adj: data fram with adjusted results
  
  ### this function produces a table that can be saved as an image or 
  ### directly to a word document!
  
  # build table
  tbl <- data.table(matrix(nrow=0, ncol=13))
  for (i in 1:length(exposure)) {
    for (j in 1:length(outcome)) {
      exp <- exposure[i]
      out <- outcome[j]
      filtered_res <- results[results$Y==out & results$X==exp,]
      filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
      unadj <- paste(round(filtered_res$`point.diff`, 2), " (", round(filtered_res$`lb.diff`, 2), ", ", round(filtered_res$`ub.diff`, 2), ")", sep="")
      adj <- paste(round(filtered_adj$`point.diff`, 2), " (", round(filtered_adj$`lb.diff`, 2), ", ", round(filtered_adj$`ub.diff`, 2), ")", sep="")
      if (j==1){
        tbl <- rbind(tbl, list(expo_var[i], out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$corrected.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$corrected.Pval, 2)))
      }else {
        tbl <- rbind(tbl, list(" ", out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$corrected.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$corrected.Pval, 2)))
      }
    }
    if (i != length(exposure)) {
      tbl <- rbind(tbl, list("","","","","","","","", "","","","",""))
    }
  }
  
  # format for export
  flextbl<-flextable(tbl, col_keys=names(tbl))
  flextbl <- set_header_labels(flextbl,
                               values = list("V1" = " ", "V2" = " ", "V3" = " ", "V4" = " ", "V5" = " ",
                                             "V6" = "Predicted Outcome at 25th Percentile", "V7" = "Predicted Outcome at 75th Percentile", "V8" = "Coefficient (95% CI)", "V9" = "P-value",
                                             "V10" = "Predicted Outcome at 25th Percentile", "V11" = "Predicted Outcome at 75th Percentile", "V12" = "Coefficient (95% CI)", "V13" = "P-value"))
  flextbl <- add_header_row(flextbl, values = c("","","","","", "Unadjusted", "Fully adjusted"), colwidths=c(1,1,1,1,1,4,4))
  # flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- add_header_row(flextbl, values = c(name, "Outcome","N","25th Percentile","75th Percentile", "Outcome, 75th Percentile v. 25th Percentile"), colwidths=c(1,1,1,1,1,8))
  # flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- hline(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- hline_bottom(flextbl, part="body", border=fp_border(color="black"))
  flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- align(flextbl, align = "center", part = "all")
  flextbl <- align(flextbl, j = c(1, 2), align = "left", part="all")
  flextbl <- autofit(flextbl, part = "all")
  flextbl <- fit_to_width(flextbl, max_width=8)
  
  flextbl
}



#### MAIN TABLES ####
# #### Table 1 ####
# # Characteristics of participants
# nperc <- function(vector){
#   n <- sum(vector==1, na.rm=T)
#   perc <- round(n/sum(!is.na(vector))*100)
#   paste(n, " (", perc, "%)", sep="")
# }
# 
# mediqr <- function(vector){
#   quantiles <- round(quantile(vector, na.rm=T), 2)
#   paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
# }
# 
# n_med_col <- c(nperc(d$sex), mediqr(d$t2_f2_8ip), mediqr(d$t2_f2_23d), mediqr(d$t2_f2_VI), mediqr(d$t2_f2_12i),
#                mediqr(d$t3_cort_slope), mediqr(d$t3_residual_cort), mediqr(d$t3_saa_slope), mediqr(d$t3_residual_saa),
#                mediqr(d$t3_map), mediqr(d$t3_hr_mean), mediqr(d$t3_gcr_mean), mediqr(d$t3_gcr_cpg12),
#                mediqr(d$laz_t2), mediqr(d$waz_t2), mediqr(d$whz_t2), mediqr(d$hcz_t2),
#                mediqr(d$laz_t3), mediqr(d$waz_t3), mediqr(d$whz_t3), mediqr(d$hcz_t3),
#                nperc(d$diar7d_t2), nperc(d$diar7d_t3), mediqr(d$momage), mediqr(d$momheight), 
#                mediqr(d$momeduy), mediqr(d$cesd_sum_t2), mediqr(d$cesd_sum_ee_t3), mediqr(d$pss_sum_mom_t3), 
#                nperc(d$life_viol_any_t3))
# 
# tbl1 <- data.table("C1" = c("Child","","","","","","","","","","","","","","","","","","","","","","","Mother","","","","","",""),
#                    "C2" = c("", "Urinary F2-isoprostanes (Year 1)","","","", "Salivary cortisol reactivity (Year 2)","", "sAA reactivity (Year 2)","",
#                            "SAM biomarkers (Year 2)","", "Glucocorticoid receptor","", "Anthropometry (14 months, Year 1)","","","",
#                            "Anthropometry (28 months, Year 2)","","","", "Diarrhea (14 months, Year 1)", "Diarrhea (28 months, Year 2)","",
#                            "Anthropometry at enrollment", "Education", "Depression at Year 1", "Depression at Year 2", "Perceived stress at Year 2", 
#                            "Intimate partner violence"),
#                    "C3" = c("Female", "iPF(2a)-III", "2,3-dinor-iPF(2a)-III", "iPF(2a-VI", "8,12-iso-iPF(2a)-VI", 
#                            "Change in slope between pre- and post-stressor cortisol", "Cortisol residualized gain score", 
#                            "Change in slope between pre- and post-stressor sAA change", "sAA residualized gain score",
#                            "Mean arterial pressure", "Resting heart rate", "NR3C1 exon 1F promoter methylation", "NGFI-A transcription factor binding site methylation",
#                            "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
#                            "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
#                            "Caregiver-reported 7-day recall", "Caregiver-reported 7-day recall", "Age (years)", "Height (cm)", "Schooling completed (years)",
#                            "CES-D score", "CES-D score", "Perceived Stress Scale score", "Any lifetime exposure"),
#                    "C4" = n_med_col)
# 
# tbl1flex <- flextable(tbl1, col_keys=names(tbl1))
# tbl1flex <- set_header_labels(tbl1flex,
#                         values = list("C1" = "", "C2" = "", "C3" = "", "C4" = "n (%) or median (IQR)"))
# tbl1flex <- hline_top(tbl1flex, part="header", border=fp_border(color="black", width = 1))
# tbl1flex <- hline_bottom(tbl1flex, part="all", border=fp_border(color="black", width = 1))
# tbl1flex <- autofit(tbl1flex, part = "all")
# tbl1flex <- align(tbl1flex, j = c(1, 2, 3), align = "left", part="all")
# tbl1flex <- align(tbl1flex, j = 4, align = "center", part="all")
# tbl1flex <- fit_to_width(tbl1flex, max_width=8)
# names(tbl1)<- c("","","","n (%) or median (IQR)")


#### Table 2 ####
# concurrent y1 immune markers and y1 who motor milestones

exposure <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
              "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "t2_ln_ifn", "sumscore_t2_Z", "t2_ln_igf") 
outcome <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2")
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10",     
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17", "Ln AGP", "Ln CRP", "Ln IFN-y", "Sum score of 13 cytokines", "Ln IGF-1") 
out_var <- c("Sum of 2nd, 4th, 5th, and 6th WHO motor milestones", "CDI comprehension Z-score", "CDI expressive language Z-score")

tbl2 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, unadj, adj)
tbl2flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, unadj, adj)

#### Table 3 ####
# concurrent y2 immune markers and y2 easq scores

exposure <- c("t3_ratio_pro_il10", "t3_ratio_il2_il10", "t3_ratio_gmc_il10", "t3_ratio_th1_il10", "t3_ratio_th2_il10",     
              "t3_ratio_th17_il10", "t3_ratio_th1_th2", "t3_ratio_th1_th17", "t3_ln_ifn", "sumscore_t3_Z", "t3_ln_igf")   
outcome <- c("z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq")  
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10",     
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17", "Ln IFN-y", "Sum score of 13 cytokines", "Ln IGF-1") 
out_var <- c("EASQ communication Z-score", "EASQ gross motor Z-score", "EASQ personal social Z-score", "EASQ combined Z-core") 
tbl3 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, unadj, adj)
tbl3flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, unadj, adj)

#### Table 4 ####
# concurrent y2 immune markers and y2 cdi scores

exposure <- c("t3_ratio_pro_il10", "t3_ratio_il2_il10", "t3_ratio_gmc_il10", "t3_ratio_th1_il10", "t3_ratio_th2_il10",     
              "t3_ratio_th17_il10", "t3_ratio_th1_th2", "t3_ratio_th1_th17", "t3_ln_ifn", "sumscore_t3_Z", "t3_ln_igf")   
outcome <- c("z_cdi_und_t3", "z_cdi_say_t3")  
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10",     
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17", "Ln IFN-y", "Sum score of 13 cytokines", "Ln IGF-1") 
out_var <- c("CDI comprehension Z-score", "CDI expressive language Z-score") 
tbl4 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, unadj, adj)
tbl4flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, unadj, adj)


#### Table 5 ####
# y1 immune markers and subsequent y2 easq scores

exposure <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
              "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "t2_ln_ifn", "sumscore_t2_Z", "t2_ln_igf") 
outcome <- c("z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq")  
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10",     
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17", "Ln AGP", "Ln CRP", "Ln IFN-y", "Sum score of 13 cytokines", "Ln IGF-1") 
out_var <- c("EASQ communication Z-score", "EASQ gross motor Z-score", "EASQ personal social Z-score", "EASQ combined Z-core") 

tbl5 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, unadj, adj)
tbl5flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, unadj, adj)

#### Table 6 ####
# y1 immune markers and subsequent y2 cdi scores

exposure <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
              "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "t2_ln_ifn", "sumscore_t2_Z", "t2_ln_igf") 
outcome <- c("z_cdi_und_t3", "z_cdi_say_t3")  
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10",     
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17", "Ln AGP", "Ln CRP", "Ln IFN-y", "Sum score of 13 cytokines", "Ln IGF-1") 
out_var <- c("CDI comprehension Z-score", "CDI expressive language Z-score") 

tbl6 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, unadj, adj)
tbl6flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, unadj, adj)



#### SAVE TABLES ####

# write.csv(tbl1, file=here("tables/main/stress-growth-table1.csv"))
write.csv(tbl2, here('tables/immune-dev-table2.csv'))
write.csv(tbl3, here('tables/immune-dev-table3.csv'))
write.csv(tbl4, here('tables/immune-dev-table4.csv'))
write.csv(tbl5, here('tables/immune-dev-table5.csv'))
write.csv(tbl6, here('tables/immune-dev-table6.csv'))

save_as_docx("Table 2: Association between Immune Status and Child Development Outcomes at Year 1" = tbl2flex, 
             "Table 3: Association between Immune Status and EASQ Scores at Year 2" = tbl3flex, 
             "Table 4: Association between Immune Status and CDI Scores at Year 2" = tbl4flex, 
             "Table 5: Association between Immune Status at Year 1 and EASQ Scores at Year 2" = tbl5flex, 
             "Table 6: Association between Immune Status at Year 1 and CDI Scores at Year 2" = tbl6flex, 
             path='C:/Users/Sophia/Documents/WASH/WASH Immune and Child Development/immune-dev main tables v2.docx')
