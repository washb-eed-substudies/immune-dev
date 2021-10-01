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
indunadj <- readRDS(here('results/unadjusted/individual_unadj_res.RDS'))

H1adj <- readRDS(here('results/adjusted/H1_adj_res.RDS'))
H2adj <- readRDS(here('results/adjusted/H2_adj_res.RDS'))
H3adj <- readRDS(here('results/adjusted/H3_adj_res.RDS'))
H4adj <- readRDS(here('results/adjusted/H4_adj_res.RDS'))
HRadj <- readRDS(here('results/adjusted/HR_adj_res.RDS'))
indadj <- readRDS(here('results/adjusted/individual_adj_res.RDS'))

unadj <- rbind(H1, H2, H3, H4, indunadj)
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
# y1 individual cytokines and development (main and supplementary)
exposure <- c("t2_ln_il12", "t2_ln_ifn", "t2_ln_il4", "t2_ln_il5", "t2_ln_il13", "sumscore_t2_Z") 
outcome <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2", "z_cdi_und_t3", "z_cdi_say_t3",
             "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq")
expo_var <- c("Ln IL-12 Year 1", "Ln IFN-y Year 1", "Ln IL-4 Year 1", "Ln IL-5 Year 1", "Ln IL-13 Year 1", "Sum score of 13 cytokines Year 1") 
out_var <- c("Sum of 2nd, 4th, 5th, and 6th WHO motor milestones Year 1", "CDI comprehension Z-score Year 1", "CDI expressive language Z-score Year 1", 
             "CDI comprehension Z-score Year 2", "CDI expressive language Z-score Year 2",
             "EASQ communication Z-score Year 2", "EASQ gross motor Z-score Year 2", "EASQ personal social Z-score Year 2", "EASQ combined Z-score Year 2")

title_y1ind <- "Association between Individual Cytokines at Year 1 and Child Development Outcomes"
y1ind <- growth_tbl("Individual Cytokines at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, T)
y1indflex <- growth_tbl_flex("Individual Cytokines at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, T, 1.1, 1.5)
y1indsupp <- growth_tbl("Individual Cytokines at Year 1", expo_var, out_var, exposure, outcome, unadj, adj,)
y1indflexsupp <- growth_tbl_flex("Individual Cytokines at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, F, 1, 1.3)

# add hazard ratios table 
exposure <- c("t2_ln_il12", "t2_ln_ifn", "t2_ln_il4", "t2_ln_il5", "t2_ln_il13", "sumscore_t2_Z") 
outcome <- c("who_crawl", "who_stand_supp", 
             "who_walk_supp", "who_stand_nosupp", "who_walk_nosup")
expo_var <- c("Ln IFN-y Year 1", "Sum score of 13 cytokines Year 1") 
out_var <- c("Hands-and-knees crawling", "Standing with assistance",
             "Walking with assistance", "Standing alone", "Walking alone")

title_y1indhr <- "Hazard Ratio for Motor Milestone Attainment for Individual Cytokines at Year 1"
y1indhr <- hr_tbl("Individual Cytokines at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj, T)
y1indhrflex <- hr_tbl_flex("Individual Cytokines at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj, T, 1.3, 1.3)
y1indhrsupp <- hr_tbl("Individual Cytokines at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj)
y1indhrflexsupp <- hr_tbl_flex("Individual Cytokines at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj, F, 1.2, 1.3)

# y2 individual cytokines and concurrent development (supplementary table only)
exposure <- c("t3_ln_il12", "t3_ln_ifn", "t3_ln_il4", "t3_ln_il5", "t3_ln_il13", "sumscore_t3_Z") 
outcome <- c("z_cdi_und_t3", "z_cdi_say_t3",
             "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq")
expo_var <- c("Ln IL-12 Year 2", "Ln IFN-y Year 2", "Ln IL-4 Year 2", "Ln IL-5 Year 2", "Ln IL-13 Year 2", "Sum score of 13 cytokines Year 2") 
out_var <- c("CDI comprehension Z-score Year 2", "CDI expressive language Z-score Year 2",
             "EASQ communication Z-score Year 2", "EASQ gross motor Z-score Year 2", "EASQ personal social Z-score Year 2", "EASQ combined Z-score Year 2")

title_y2ind <- "Association between Individual Cytokines at Year 2 and Child Development Outcomes"
y2indsupp <- growth_tbl("Individual Cytokines at Year 2", expo_var, out_var, exposure, outcome, unadj, adj,)
y2indflexsupp <- growth_tbl_flex("Individual Cytokines at Year 2", expo_var, out_var, exposure, outcome, unadj, adj, F, 1, 1.3)

# y1 cytokine ratios and development (main and supplementary)
exposure <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", 
              "t2_ratio_th1_il10", "t2_ratio_il12_il10", "t2_ratio_ifn_il10", 
              "t2_ratio_th2_il10", "t2_ratio_il4_il10", "t2_ratio_il5_il10", "t2_ratio_il13_il10",
              "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17") 
outcome <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2", "z_cdi_und_t3", "z_cdi_say_t3",
             "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq")
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10 Year 1", "Ln IL-2/IL-10 Year 1", "Ln GM-CSF/IL-10 Year 1", 
              "Ln Th1/IL-10 Year 1", "Ln IL-12/IL-10 Year 1", "Ln IFN-y/IL-10 Year 1",
              "Ln Th2/IL-10 Year 1", "Ln IL-4/IL-10 Year 1", "Ln IL-5/IL-10 Year 1", "Ln IL-13/IL-10 Year 1",
              "Ln Th17/IL-10 Year 1", "Ln Th1/Th2 Year 1", "Ln Th1/Th17 Year 1") 
out_var <- c("Sum of 2nd, 4th, 5th, and 6th WHO motor milestones Year 1", "CDI comprehension Z-score Year 1", "CDI expressive language Z-score Year 1",
             "CDI comprehension Z-score Year 2", "CDI expressive language Z-score Year 2",
             "EASQ communication Z-score Year 2", "EASQ gross motor Z-score Year 2", "EASQ personal social Z-score Year 2", "EASQ combined Z-score Year 2")

title_y1rat <- "Association between Cytokine Ratios at Year 1 and Child Development Outcomes"
y1rat <- growth_tbl("Cytokine Ratios at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, T)
y1ratflex <- growth_tbl_flex("Cytokine Ratios at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, T, 1.1, 1.5)
y1ratsupp <- growth_tbl("Cytokine Ratios at Year 1", expo_var, out_var, exposure, outcome, unadj, adj,)
y1ratflexsupp <- growth_tbl_flex("Cytokine Ratios at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, F, 1, 1.3)

# concurrent y1 immune markers and y1 who motor milestones hazard ratios
exposure <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
              "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17") 
outcome <- c("who_crawl", "who_stand_supp", 
             "who_walk_supp", "who_stand_nosupp", "who_walk_nosup")
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10 Year 1", "Ln IL-2/IL-10 Year 1", "Ln GM-CSF/IL-10 Year 1", 
              "Ln Th1/IL-10 Year 1", "Ln Th2/IL-10 Year 1",     
              "Ln Th17/IL-10 Year 1", "Ln Th1/Th2 Year 1", "Ln Th1/Th17 Year 1") 
out_var <- c("Hands-and-knees crawling", "Standing with assistance",
             "Walking with assistance", "Standing alone", "Walking alone")

title_y1rathr <- "Hazard Ratio for Motor Milestone Attainment for Cytokine Ratios at Year 1"
y1rathr <- hr_tbl("Cytokine Ratios at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj, T)
y1rathrflex <- hr_tbl_flex("Cytokine Ratios at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj, T, 1.3, 1.3)
y1rathrsupp <- hr_tbl("Cytokine Ratios at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj)
y1rathrflexsupp <- hr_tbl_flex("Cytokine Ratios at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj, F, 1.2, 1.3)

# y2 cytokine ratios and concurrent development (supplementary only)
exposure <- c("t3_ratio_pro_il10", "t3_ratio_il2_il10", "t3_ratio_gmc_il10", 
              "t3_ratio_th1_il10", "t3_ratio_il12_il10", "t3_ratio_ifn_il10", 
              "t3_ratio_th2_il10", "t3_ratio_il4_il10", "t3_ratio_il5_il10", "t3_ratio_il13_il10",
              "t3_ratio_th17_il10", "t3_ratio_th1_th2", "t3_ratio_th1_th17") 
outcome <- c("z_cdi_und_t3", "z_cdi_say_t3",
             "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq")
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10 Year 2", "Ln IL-2/IL-10 Year 2", "Ln GM-CSF/IL-10 Year 2", 
              "Ln Th1/IL-10 Year 2", "Ln IL-12/IL-10 Year 2", "Ln IFN-y/IL-10 Year 2",
              "Ln Th2/IL-10 Year 2", "Ln IL-4/IL-10 Year 2", "Ln IL-5/IL-10 Year 2", "Ln IL-13/IL-10 Year 2",
              "Ln Th17/IL-10 Year 2", "Ln Th1/Th2 Year 2", "Ln Th1/Th17 Year 2") 
out_var <- c("CDI comprehension Z-score Year 2", "CDI expressive language Z-score Year 2",
             "EASQ communication Z-score Year 2", "EASQ gross motor Z-score Year 2", "EASQ personal social Z-score Year 2", "EASQ combined Z-score Year 2")

title_y2rat <- "Association between Cytokine Ratios at Year 2 and Child Development Outcomes"
y2ratsupp <- growth_tbl("Cytokine Ratios at Year 2", expo_var, out_var, exposure, outcome, unadj, adj,)
y2ratflexsupp <- growth_tbl_flex("Cytokine Ratios at Year 2", expo_var, out_var, exposure, outcome, unadj, adj, F, .9, 1.1)

# y1 crp/agp and development (main and supplementary)
exposure <- c("t2_ln_crp", "t2_ln_agp") 
outcome <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2", "z_cdi_und_t3", "z_cdi_say_t3",
             "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq")
expo_var <- c("Ln CRP Year 1", "Ln AGP Year 1") 
out_var <- c("Sum of 2nd, 4th, 5th, and 6th WHO motor milestones Year 1", "CDI comprehension Z-score Year 1", "CDI expressive language Z-score Year 1",
             "CDI comprehension Z-score Year 2", "CDI expressive language Z-score Year 2",
             "EASQ communication Z-score Year 2", "EASQ gross motor Z-score Year 2", "EASQ personal social Z-score Year 2", "EASQ combined Z-score Year 2")

title_y1crpagp <- "Association between CRP and AGP at Year 1 and Child Development Outcomes"
y1crpagp <- growth_tbl("CRP and AGP at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, T)
y1crpagpflex <- growth_tbl_flex("CRP and AGP at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, T, 1, 1.5)
y1crpagpsupp <- growth_tbl("CRP and AGP at Year 1", expo_var, out_var, exposure, outcome, unadj, adj,)
y1crpagpflexsupp <- growth_tbl_flex("CRP and AGP at Year 1", expo_var, out_var, exposure, outcome, unadj, adj, F, 1, 1.3)

# crp/agp who motor milestones hazard ratios
outcome <- c("who_crawl", "who_stand_supp", 
             "who_walk_supp", "who_stand_nosupp", "who_walk_nosup")
out_var <- c("Hands-and-knees crawling", "Standing with assistance",
             "Walking with assistance", "Standing alone", "Walking alone")

title_y1crpagphr <- "Hazard Ratio for Motor Milestone Attainment for CRP and AGP at Year 1"
y1crpagphr <- hr_tbl("Cytokine Ratios at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj, T)
y1crpagphrflex <- hr_tbl_flex("Cytokine Ratios at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj, T, 1.3, 1.3)
y1crpagphrsupp <- hr_tbl("Cytokine Ratios at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj)
y1crpagphrflexsupp <- hr_tbl_flex("Cytokine Ratios at Year 1", expo_var, out_var, exposure, outcome, HR, HRadj, F, 1.2, 1.3)

# y2 crp/agp and concurrent development (supplementary only)
exposure <- c("t3_ln_crp", "t3_ln_agp") 
outcome <- c("z_cdi_und_t3", "z_cdi_say_t3",
             "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq")
expo_var <- c("Ln CRP Year 2", "Ln AGP Year 2") 
out_var <- c("CDI comprehension Z-score Year 2", "CDI expressive language Z-score Year 2",
             "EASQ communication Z-score Year 2", "EASQ gross motor Z-score Year 2", "EASQ personal social Z-score Year 2", "EASQ combined Z-score Year 2")

title_y2crpagp <- "Association between CRP and AGP at Year 2 and Child Development Outcomes"
y2crpagpsupp <- growth_tbl("CRP and AGP at Year 2", expo_var, out_var, exposure, outcome, unadj, adj,)
y2crpagpflexsupp <- growth_tbl_flex("CRP and AGP at Year 2", expo_var, out_var, exposure, outcome, unadj, adj, F, .9, 1.1)

#### Table 8 ####
# IGF and all development outcomes
exposure <- c("t2_ln_igf", "t3_ln_igf") 
outcome <- c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2", "z_cdi_und_t3", "z_cdi_say_t3",
             "z_comm_easq", "z_motor_easq", "z_personal_easq", "z_combined_easq")  
expo_var <- c("Ln IGF-1 Year 1", "Ln IGF-1 Year 2") 
out_var <- c("Sum of 2nd, 4th, 5th, and 6th WHO motor milestones", "CDI comprehension Z-score Year 1", "CDI expressive language Z-score Year 1",
             "CDI comprehension Z-score Year 2", "CDI expressive language Z-score Year 2", 
             "EASQ communication Z-score Year 2", "EASQ gross motor Z-score Year 2", "EASQ personal social Z-score Year 2", "EASQ combined Z-score Year 2") 

title_igf <- "Association between IGF-1 and Child Development Outcomes"
igf <- growth_tbl("IGF-1", expo_var, out_var, exposure, outcome, unadj, adj, T)
igfflex <- growth_tbl_flex("IGF-1", expo_var, out_var, exposure, outcome, unadj, adj, T, .5, 1)
igfsupp <- growth_tbl("IGF-1", expo_var, out_var, exposure, outcome, unadj, adj)
igfflexsupp <- growth_tbl_flex("IGF-1", expo_var, out_var, exposure, outcome, unadj, adj, F, .5, 1)

#### Table 9 ####
# IGF and hazard ratios

exposure <- c("t2_ln_igf") 
outcome <- c("who_crawl", "who_stand_supp", 
             "who_walk_supp", "who_stand_nosupp", "who_walk_nosup")
expo_var <- c("Ln IGF-1 Year 1") 
out_var <- c("Hands-and-knees crawling", "Standing with assistance",
             "Walking with assistance", "Standing alone", "Walking alone")

title_igfhr <- "Hazard Ratio for Motor Milestone Attainment for IGF-1 at Year 1"
igfhrsupp <- hr_tbl("IGF-1", expo_var, out_var, exposure, outcome, HR, HRadj)
igfhrflexsupp <- hr_tbl_flex("IGF-1", expo_var, out_var, exposure, outcome, HR, HRadj, F, .5, .9)


ls()

#### SAVE TABLES ####

write.csv(y1crpagp, here('tables/main/immune-dev-y1-crp-agp.csv'))
write.csv(y1crpagphr, here('tables/main/immune-dev-y1-crp-agp-hr.csv'))
write.csv(y1ind, here('tables/main/immune-dev-y1-ind-cyt.csv'))
write.csv(y1indhr, here('tables/main/immune-dev-y1-ind-cyt-hr.csv'))
write.csv(y1rat, here('tables/main/immune-dev-y1-ratios.csv'))
write.csv(y1rathr, here('tables/main/immune-dev-y1-ratios-hr.csv'))
write.csv(igf, here('tables/main/immune-dev-igf.csv'))

write.csv(y1crpagpsupp, here('tables/supplementary/immune-dev-y1-crp-agp.csv'))
write.csv(y1crpagphrsupp, here('tables/supplementary/immune-dev-y1-crp-agp-hr.csv'))
write.csv(y1indsupp, here('tables/supplementary/immune-dev-y1-ind-cyt.csv'))
write.csv(y1indhrsupp, here('tables/supplementary/immune-dev-y1-ind-cyt-hr.csv'))
write.csv(y1ratsupp, here('tables/supplementary/immune-dev-y1-ratios.csv'))
write.csv(y1rathrsupp, here('tables/supplementary/immune-dev-y1-ratios-hr.csv'))
write.csv(y2crpagpsupp, here('tables/supplementary/immune-dev-y2-crp-agp.csv'))
write.csv(y2indsupp, here('tables/supplementary/immune-dev-y2-ind-cyt.csv'))
write.csv(y2ratsupp, here('tables/supplementary/immune-dev-y2-ratios.csv'))
write.csv(igfsupp, here('tables/supplementary/immune-dev-igf.csv'))
write.csv(igfhrsupp, here('tables/supplementary/immune-dev-igfhr.csv'))


save_as_docx("Table 1: Characteristics of Participants" = enroll, 
             "Table 2" = y1indhrflex, 
             "Table 3" = y1indflex, 
             "Table 4" = y1rathrflex, 
             "Table 5" = y1ratflex, 
             "Table 6" = y1crpagphrflex, 
             "Table 7" = y1crpagpflex, 
             "Table 8" = igfflex, 
             path='/Users/sophiatan/Documents/WASH/immune-dev main tables v12.docx', 
             pr_section = sect_properties)

save_as_docx("Table 1" = y1indhrflexsupp, 
             "Table 2" = y1indflexsupp, 
             "Table 3" = y2indflexsupp, 
             "Table 4" = y1rathrflexsupp, 
             "Table 5" = y1ratflexsupp, 
             "Table 6" = y2ratflexsupp, 
             "Table 7" = y1crpagphrflexsupp, 
             "Table 8" = y1crpagpflexsupp, 
             "Table 9" = y2crpagpflexsupp, 
             "Table 10" = igfhrflexsupp, 
             "Table 11" = igfflexsupp, 
             path='/Users/sophiatan/Documents/WASH/immune-dev supplementary v12.docx',
             pr_section = sect_properties)
