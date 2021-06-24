rm(list=ls())

source(here::here("0-config.R"))
source(here::here("table-functions.R"))

# load enrollment characteristics and results
d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-development-analysis-dataset.rds"))
H1 <- readRDS(here('results/unadjusted/H1_res.RDS'))
H2 <- readRDS(here('results/unadjusted/H2_res.RDS'))
H3 <- readRDS(here('results/unadjusted/H3_res.RDS'))
H4 <- readRDS(here('results/unadjusted/H4_res.RDS'))
HR <- readRDS(here('results/unadjusted/HR_res.RDS'))

H1adj <- readRDS(here('results/adjusted/H1_adj_res.RDS'))
H2adj <- readRDS(here('results/adjusted/H2_adj_res.RDS'))
H3adj <- readRDS(here('results/adjusted/H3_adj_res.RDS'))
H4adj <- readRDS(here('results/adjusted/H4_adj_res.RDS'))
HRadj <- readRDS(here('results/adjusted/HR_adj_res.RDS'))
indadj <- readRDS(here('results/adjusted/individual_adj_res.RDS'))

unadj <- rbind(H1, H2, H3, H4)
adj <- rbind(H1adj, H2adj, H3adj, H4adj, indadj)

#### MAIN TABLES ####
# #### Table 1 ####
filtering <- function(row){
  any(!is.na(row))
}

exp_t2 <- c(grep("t2_ln", names(d), value=TRUE), "sumscore_t2_Z")
out_t2 <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2")
exp_t3 <- c(grep("t3_ln", names(d), value=TRUE), "sumscore_t3_Z")
out_t3 <- c("z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq", 
            "z_cdi_say_t3", "z_cdi_und_t3")
t2 <- d[apply(select(d, all_of(exp_t2)), 1, filtering),] # only has rows where we have exposure data at t2
t2 <- t2[apply(select(t2, all_of(c(out_t2, out_t3))), 1, filtering),] # only has rows where we have both some exposure data and some outcome data
t3 <- d[apply(select(d, all_of(exp_t3)), 1, filtering),] # only has rows where we have exposure data at t3
t3 <- t3[apply(select(t3, all_of(out_t3)), 1, filtering),] # only has rows where we have both some exposure data and some outcome data
d <- t2 %>% full_join(t3, by=names(d)) # has all children included in this analysis

characteristics <- function(d, child_char = NULL, child_char_names = NULL, mom_char = NULL, mom_char_names = NULL) {
  nperc <- function(vector){
    n <- sum(vector==1, na.rm=T)
    perc <- round(n/sum(!is.na(vector))*100)
    paste(n, " (", perc, "%)", sep="")
  }
  
  mediqr <- function(vector){
    quantiles <- round(quantile(vector, na.rm=T), 2)
    paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
  }
  
  child <- c('sex',child_char,'laz_t1','waz_t1','whz_t1','hcz_t1',
             'laz_t2','waz_t2','whz_t2','hcz_t2',
             'laz_t3','waz_t3','whz_t3','hcz_t3','diar7d_t2','diar7d_t3')
  
  mom <- c('momage', 'momheight', 'momeduy', mom_char, 'cesd_sum_t2', 'cesd_sum_ee_t3', 'pss_sum_mom_t3', 'life_viol_any_t3')
  
  n_med_col <- NULL
  for (var in c(child, mom)) {
    if (var %in% c('sex', 'diar7d_t2', 'diar7d_t3', 'life_viol_any_t3') | is.factor(d[[var]])) {
      if (var == 'sex') {
        n <- sum(d$sex=='female', na.rm=T)
        perc <- round(n/sum(!is.na(d$sex))*100)
        n_med_col <- c(n_med_col, paste(n, " (", perc, "%)", sep=""))
      }else {
        d[[var]] <- na_if(d[[var]], "Missing")
        n_med_col <- c(n_med_col, nperc(d[[var]]))
      }
    }else {
      n_med_col <- c(n_med_col, mediqr(d[[var]]))
    }
  }
  
  tbl1 <- data.table("C1" = c("Child", rep("", length(child)-1),"Mother", rep("",length(mom)-1)),
                     "C2" = c("", rep("", length(child_char)), "Anthropometry (3 months, Year 1)","","","",
                              "Anthropometry (14 months, Year 1)","","","",
                              "Anthropometry (28 months, Year 2)","","","", "Diarrhea (14 months, Year 1)", "Diarrhea (28 months, Year 2)","",
                              "Anthropometry at enrollment", "Education", rep("", length(mom_char)), "Depression at Year 1", "Depression at Year 2", "Perceived stress at Year 2", 
                              "Intimate partner violence"),
                     "C3" = c("Female", child_char_names,
                              "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
                              "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
                              "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
                              "Caregiver-reported 7-day recall", "Caregiver-reported 7-day recall", "Age (years)", "Height (cm)", "Schooling completed (years)",
                              mom_char_names, "CES-D score", "CES-D score", "Perceived Stress Scale score", "Any lifetime exposure"),
                     "C4" = n_med_col)
  
  tbl1flex <- flextable(tbl1, col_keys=names(tbl1))
  tbl1flex <- set_header_labels(tbl1flex,
                                values = list("C1" = "", "C2" = "", "C3" = "", "C4" = "n (%) or median (IQR)"))
  tbl1flex <- hline_top(tbl1flex, part="header", border=fp_border(color="black", width = 1))
  tbl1flex <- hline_bottom(tbl1flex, part="all", border=fp_border(color="black", width = 1))
  tbl1flex <- autofit(tbl1flex, part = "all")
  tbl1flex <- align(tbl1flex, j = c(1, 2, 3), align = "left", part="all")
  tbl1flex <- align(tbl1flex, j = 4, align = "center", part="all")
  tbl1flex <- fit_to_width(tbl1flex, max_width=8)
  tbl1flex
}

enroll <- characteristics(d)


#### Table 2 ####
# concurrent y1 immune markers and y1 who motor milestones

exposure <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", 
              "t2_ratio_th1_il10", "t2_ratio_il12_il10", "t2_ratio_ifn_il10", "t2_ln_il12", "t2_ln_ifn", 
              "t2_ratio_th2_il10", "t2_ratio_il4_il10", "t2_ratio_il5_il10", "t2_ratio_il13_il10", "t2_ln_il4", "t2_ln_il5", "t2_ln_il13",
              "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", 
              "sumscore_t2_Z", "t2_ln_agp", "t2_ln_crp") 
outcome <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2")
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", 
              "Ln Th1/IL-10", "Ln IL-12/IL-10", "Ln IFN-y/IL-10", "Ln IL-12", "Ln IFN-y",
              "Ln Th2/IL-10", "Ln IL-4/IL-10", "Ln IL-5/IL-10", "Ln IL-13/IL-10",
              "Ln IL-4", "Ln IL-5", "Ln IL-13",
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17",
              "Sum score of 13 cytokines", "Ln AGP", "Ln CRP") 
out_var <- c("Sum of 2nd, 4th, 5th, and 6th WHO motor milestones", "CDI comprehension Z-score", "CDI expressive language Z-score")

tbl2 <- growth_tbl("Immune Status at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, T)
tbl2flex <- growth_tbl_flex("Immune Status at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, T, 1.3, 1.5)
tbl2supp <- growth_tbl("Immune Status at Year 1", expo_var, out_var, exposure, outcome, unadj, adj)
tbl2flexsupp <- growth_tbl_flex("Immune Status at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, F, .9, 1.1)


#### Table 3 ####
# concurrent y1 immune markers and y1 who motor milestones hazard ratios

exposure <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
              "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_ifn", "sumscore_t2_Z", "t2_ln_agp", "t2_ln_crp") 
outcome <- c("who_crawl", "who_stand_supp", 
             "who_walk_supp", "who_stand_nosupp", "who_walk_nosup")
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10",     
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17", "Ln IFN-y", "Sum score of 13 cytokines", "Ln AGP", "Ln CRP") 
out_var <- c("Hands-and-knees crawling", "Standing with assistance",
             "Walking with assistance", "Standing alone", "Walking alone")

tbl3 <- hr_tbl("Immune Status at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj, T)
tbl3flex <- hr_tbl_flex("Immune Status at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj, T, 1.3, 1.3)
tbl3supp <- hr_tbl("Immune Status at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj)
tbl3flexsupp <- hr_tbl_flex("Immune Status at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj, F, 1.2, 1.3)

#### Table 4 ####
# concurrent y2 immune markers and y2 easq scores

exposure <- c("t3_ratio_pro_il10", "t3_ratio_il2_il10", "t3_ratio_gmc_il10", "t3_ratio_th1_il10", "t3_ratio_th2_il10",     
              "t3_ratio_th17_il10", "t3_ratio_th1_th2", "t3_ratio_th1_th17", "t3_ln_ifn", "sumscore_t3_Z", "t3_ln_agp", "t3_ln_crp")   
outcome <- c("z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq")  
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10",     
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17","Ln IFN-y", "Sum score of 13 cytokines", "Ln AGP", "Ln CRP") 
out_var <- c("EASQ communication Z-score", "EASQ gross motor Z-score", "EASQ personal social Z-score", "EASQ combined Z-score") 

tbl4 <- growth_tbl("Immune Status at Year 2", expo_var, out_var, exposure, outcome, unadj, adj, T)
tbl4flex <- growth_tbl_flex("Immune Status at Year 2", expo_var, out_var, exposure, outcome, unadj, adj, T, 1.3, 1.5)
tbl4supp <- growth_tbl("Immune Status at Year 2", expo_var, out_var, exposure, outcome, unadj, adj)
tbl4flexsupp <- growth_tbl_flex("Immune Status at Year 2", expo_var, out_var, exposure, outcome, unadj, adj, F, .9,1.1)

#### Table 5 ####
# concurrent y2 immune markers and y2 cdi scores

exposure <- c("t3_ratio_pro_il10", "t3_ratio_il2_il10", "t3_ratio_gmc_il10", "t3_ratio_th1_il10", "t3_ratio_th2_il10",     
              "t3_ratio_th17_il10", "t3_ratio_th1_th2", "t3_ratio_th1_th17", "t3_ln_ifn", "sumscore_t3_Z", "t3_ln_agp", "t3_ln_crp")   
outcome <- c("z_cdi_und_t3", "z_cdi_say_t3")  
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10",     
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17", "Ln IFN-y", "Sum score of 13 cytokines", "Ln AGP", "Ln CRP") 
out_var <- c("CDI comprehension Z-score", "CDI expressive language Z-score") 

tbl5 <- growth_tbl("Immune Status at Year 2", expo_var, out_var, exposure, outcome, unadj, adj, T)
tbl5flex <- growth_tbl_flex("Immune Status at Year 2", expo_var, out_var, exposure, outcome, unadj, adj, T, 1.3, 1.6)
tbl5supp <- growth_tbl("Immune Status at Year 2", expo_var, out_var, exposure, outcome, unadj, adj)
tbl5flexsupp <- growth_tbl_flex("Immune Status at Year 2", expo_var, out_var, exposure, outcome, unadj, adj, F, .9, 1.1)

#### Table 6 ####
# y1 immune markers and subsequent y2 easq scores

exposure <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", 
              "t2_ratio_th1_il10", "t2_ratio_il12_il10", "t2_ratio_ifn_il10", "t2_ln_il12", "t2_ln_ifn", 
              "t2_ratio_th2_il10", "t2_ratio_il4_il10", "t2_ratio_il5_il10", "t2_ratio_il13_il10", "t2_ln_il4", "t2_ln_il5", "t2_ln_il13",
              "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", 
              "sumscore_t2_Z", "t2_ln_agp", "t2_ln_crp") 
outcome <- c("z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq")  
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", 
              "Ln Th1/IL-10", "Ln IL-12/IL-10", "Ln IFN-y/IL-10", "Ln IL-12", "Ln IFN-y",
              "Ln Th2/IL-10", "Ln IL-4/IL-10", "Ln IL-5/IL-10", "Ln IL-13/IL-10",
              "Ln IL-4", "Ln IL-5", "Ln IL-13",
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17",
              "Sum score of 13 cytokines", "Ln AGP", "Ln CRP") 
out_var <- c("EASQ communication Z-score", "EASQ gross motor Z-score", "EASQ personal social Z-score", "EASQ combined Z-score") 

tbl6 <- growth_tbl("Immune Status at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, T)
tbl6flex <- growth_tbl_flex("Immune Status at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, T, 1.3, 1.5)
tbl6supp <- growth_tbl("Immune Status at Year 1", expo_var, out_var, exposure, outcome, unadj, adj)
tbl6flexsupp <- growth_tbl_flex("Immune Status at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, F, .8, 1)

#### Table 7 ####
# y1 immune markers and subsequent y2 cdi scores

exposure <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", 
              "t2_ratio_th1_il10", "t2_ratio_il12_il10", "t2_ratio_ifn_il10", "t2_ln_il12", "t2_ln_ifn", 
              "t2_ratio_th2_il10", "t2_ratio_il4_il10", "t2_ratio_il5_il10", "t2_ratio_il13_il10", "t2_ln_il4", "t2_ln_il5", "t2_ln_il13",
              "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", 
              "sumscore_t2_Z", "t2_ln_agp", "t2_ln_crp") 
outcome <- c("z_cdi_und_t3", "z_cdi_say_t3")  
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", 
              "Ln Th1/IL-10", "Ln IL-12/IL-10", "Ln IFN-y/IL-10", "Ln IL-12", "Ln IFN-y",
              "Ln Th2/IL-10", "Ln IL-4/IL-10", "Ln IL-5/IL-10", "Ln IL-13/IL-10",
              "Ln IL-4", "Ln IL-5", "Ln IL-13",
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17",
              "Sum score of 13 cytokines", "Ln AGP", "Ln CRP") 
out_var <- c("CDI comprehension Z-score", "CDI expressive language Z-score") 

tbl7 <- growth_tbl("Immune Status at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, T)
tbl7flex <- growth_tbl_flex("Immune Status at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, T, 1.3, 1.6)
tbl7supp <- growth_tbl("Immune Status at Year 1", expo_var, out_var, exposure, outcome, unadj, adj)
tbl7flexsupp <- growth_tbl_flex("Immune Status at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, F, .9, 1.1)

#### Table 8 ####
# IGF and all development outcomes

exposure <- c("t2_ln_igf", "t3_ln_igf") 
outcome <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2", "z_cdi_und_t3", "z_cdi_say_t3",
             "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq")  
expo_var <- c("Ln IGF-1 Year 1", "Ln IGF-1 Year 2") 
out_var <- c("Sum of 2nd, 4th, 5th, and 6th WHO motor milestones", "CDI comprehension Z-score Year 1", "CDI expressive language Z-score Year 1",
             "CDI comprehension Z-score Year 2", "CDI expressive language Z-score Year 2", 
             "EASQ communication Z-score Year 2", "EASQ gross motor Z-score Year 2", "EASQ personal social Z-score Year 2", "EASQ combined Z-score Year 2") 

tbl8 <- growth_tbl("IGF-1", expo_var, out_var, exposure, outcome, unadj, adj, T)
tbl8flex <- growth_tbl_flex("IGF-1", expo_var, out_var, exposure, outcome, unadj, adj, T, .5, 1.6)
tbl8supp <- growth_tbl("IGF-1", expo_var, out_var, exposure, outcome, unadj, adj)
tbl8flexsupp <- growth_tbl_flex("IGF-1", expo_var, out_var, exposure, outcome, unadj, adj, F, .5, 1)

#### Table 9 ####
# IGF and hazard ratios

exposure <- c("t2_ln_igf") 
outcome <- c("who_crawl", "who_stand_supp", 
             "who_walk_supp", "who_stand_nosupp", "who_walk_nosup")
expo_var <- c("Ln IGF-1 Year 1") 
out_var <- c("Hands-and-knees crawling", "Standing with assistance",
             "Walking with assistance", "Standing alone", "Walking alone")

tbl9 <- hr_tbl("IGF-1", expo_var, out_var, exposure, outcome, HR, HRadj, T)
tbl9flex <- hr_tbl_flex("IGF-1", expo_var, out_var, exposure, outcome, HR, HRadj, T, .8, 1.5)
tbl9supp <- hr_tbl("IGF-1", expo_var, out_var, exposure, outcome, HR, HRadj)
tbl9flexsupp <- hr_tbl_flex("IGF-1", expo_var, out_var, exposure, outcome, HR, HRadj, F, .5, .9)



#### SAVE TABLES ####

#write.csv(tbl1, file=here('tables/main/immune-dev-table1.csv'))
write.csv(tbl2, here('tables/main/immune-dev-table2.csv'))
write.csv(tbl3, here('tables/main/immune-dev-table3.csv'))
write.csv(tbl4, here('tables/main/immune-dev-table4.csv'))
write.csv(tbl5, here('tables/main/immune-dev-table5.csv'))
write.csv(tbl6, here('tables/main/immune-dev-table6.csv'))
write.csv(tbl7, here('tables/main/immune-dev-table7.csv'))
write.csv(tbl8, here('tables/main/immune-dev-table8.csv'))
write.csv(tbl9, here('tables/main/immune-dev-table9.csv'))

write.csv(tbl2supp, here('tables/supplementary/immune-dev-table1.csv'))
write.csv(tbl3supp, here('tables/supplementary/immune-dev-table2.csv'))
write.csv(tbl4supp, here('tables/supplementary/immune-dev-table3.csv'))
write.csv(tbl5supp, here('tables/supplementary/immune-dev-table4.csv'))
write.csv(tbl6supp, here('tables/supplementary/immune-dev-table5.csv'))
write.csv(tbl7supp, here('tables/supplementary/immune-dev-table6.csv'))
write.csv(tbl8supp, here('tables/supplementary/immune-dev-table7.csv'))
write.csv(tbl9supp, here('tables/supplementary/immune-dev-table8.csv'))


save_as_docx("Table 1: Characteristics of Participants" = enroll, 
             "Table 2: Association between Immune Status and Child Development Outcomes at Year 1" = tbl2flex, 
             "Table 3: Hazard Ratio for Motor Milestone Attainment for Immune Status at Year 1" = tbl3flex, 
             "Table 4: Association between Immune Status and EASQ Scores at Year 2" = tbl4flex, 
             "Table 5: Association between Immune Status and CDI Scores at Year 2" = tbl5flex, 
             "Table 6: Association between Immune Status at Year 1 and EASQ Scores at Year 2" = tbl6flex, 
             "Table 7: Association between Immune Status at Year 1 and CDI Scores at Year 2" = tbl7flex, 
             "Table 8: Association between IGF-1 and Child Development" = tbl8flex, 
             "Table 9: Hazard Ratio for Motor Milestone Attainment for IGF-1" = tbl9flex, 
             path='C:/Users/Sophia/Documents/WASH/WASH Immune and Child Development/immune-dev main tables v10.docx', 
             pr_section = sect_properties)

save_as_docx("Table 1: Association between Immune Status and Child Development Outcomes at Year 1" = tbl2flexsupp, 
             "Table 2: Hazard Ratio for Motor Milestone Attainment for Immune Status at Year 1" = tbl3flexsupp, 
             "Table 3: Association between Immune Status and EASQ Scores at Year 2" = tbl4flexsupp, 
             "Table 4: Association between Immune Status and CDI Scores at Year 2" = tbl5flexsupp, 
             "Table 5: Association between Immune Status at Year 1 and EASQ Scores at Year 2" = tbl6flexsupp, 
             "Table 6: Association between Immune Status at Year 1 and CDI Scores at Year 2" = tbl7flexsupp, 
             "Table 7: Association between IGF-1 and Child Development" = tbl8flexsupp, 
             "Table 8: Hazard Ratio for Motor Milestone Attainment for IGF-1" = tbl9flexsupp, 
             path='C:/Users/Sophia/Documents/WASH/WASH Immune and Child Development/immune-dev supplementary v9.docx',
             pr_section = sect_properties)
