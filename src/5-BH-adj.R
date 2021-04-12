#### Adjust all pvalues with BH procedure ####
rm(list=ls())

source(here::here("0-config.R"))

# load all results
H1 <- readRDS(here('results/unadjusted/H1_res.RDS')) %>% select(X, Y, Pval)
H2 <- readRDS(here('results/unadjusted/H2_res.RDS')) %>% select(X, Y, Pval)
H3 <- readRDS(here('results/unadjusted/H3_res.RDS')) %>% select(X, Y, Pval)
H4 <- readRDS(here('results/unadjusted/H4_res.RDS')) %>% select(X, Y, Pval)
HR <- readRDS(here('results/unadjusted/HR_res.RDS')) %>% select(X, Y, Pval)

H1_adj <- readRDS(here('results/adjusted/H1_adj_res.RDS')) %>% select(X, Y, Pval)
H2_adj <- readRDS(here('results/adjusted/H2_adj_res.RDS')) %>% select(X, Y, Pval)
H3_adj <- readRDS(here('results/adjusted/H3_adj_res.RDS')) %>% select(X, Y, Pval)
H4_adj <- readRDS(here('results/adjusted/H4_adj_res.RDS')) %>% select(X, Y, Pval)
HR_adj <- readRDS(here('results/adjusted/HR_adj_res.RDS')) %>% select(X, Y, Pval)
ind_adj <- readRDS(here('results/adjusted/individual_adj_res.RDS')) %>% select(X, Y, Pval)

H1_res <- H1 %>% select(X, Y, Pval)
H2_res <- H2 %>% select(X, Y, Pval)
H3_res <- H3 %>% select(X, Y, Pval)
H4_res <- H4 %>% select(X, Y, Pval)
HR_res <- HR %>% select(X, Y, Pval)

H1_adj_res <- H1_adj %>% select(X, Y, Pval)
H2_adj_res <- H2_adj %>% select(X, Y, Pval)
H3_adj_res <- H3_adj %>% select(X, Y, Pval)
H4_adj_res <- H4_adj %>% select(X, Y, Pval)
HR_adj_res <- HR_adj %>% select(X, Y, Pval)
ind_adj_res <- ind_adj %>% select(X, Y, Pval)


H1_res$H = 1
H2_res$H = 2
H3_res$H = 3
H4_res$H = 4
HR_res$H = 5

H1_adj_res$H = 1
H2_adj_res$H = 2
H3_adj_res$H = 3
H4_adj_res$H = 4
HR_adj_res$H = 5
ind_adj_res$H = 6

H1_res$G = 1
H2_res$G = 1
H3_res$G = 1
H4_res$G = 2
HR_res$G = if_else(grepl("igf", HR_res$X), 2, 1)

H1_adj_res$G = 1
H2_adj_res$G = 1
H3_adj_res$G = 1
H4_adj_res$G = 2
HR_adj_res$G = if_else(grepl("igf", HR_adj_res$X), 2, 1)
ind_adj_res$G = 2



H1_res$time <- ifelse(H1_res$Y %in% c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2"), "t2C", "t3C")
H2_res$time <- "t3S"
H3_res$time <- ifelse(H3_res$Y %in% c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2"), "t2C", 
                      ifelse(grepl("t2", H3_res$X), "t3S", "t3C"))
H4_res$time <- ifelse(H4_res$Y %in% c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2"), "t2C", 
                      ifelse(grepl("t2", H4_res$X), "t3S", "t3C"))
HR_res$time <- "t2C"

H1_adj_res$time <- ifelse(H1_adj_res$Y %in% c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2"), "t2C", "t3C")
H2_adj_res$time <- "t3S"
H3_adj_res$time <- ifelse(H3_adj_res$Y %in% c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2"), "t2C", 
                      ifelse(grepl("t2", H3_adj_res$X), "t3S", "t3C"))
H4_adj_res$time <- ifelse(H4_adj_res$Y %in% c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2"), "t2C", 
                      ifelse(grepl("t2", H4_adj_res$X), "t3S", "t3C"))
HR_adj_res$time <- "t2C"
ind_adj_res$time <- ifelse(ind_adj_res$Y %in% c("sum_who", "z_cdi_und_t2", "z_cdi_say_t2"), "t2C", "t3S")


full_res <- rbind(H1_res, H2_res, H3_res, H4_res, HR_res)
full_adj_res <- rbind(H1_adj_res, H2_adj_res, H3_adj_res, H4_adj_res, HR_adj_res, ind_adj_res)

full_res <- full_res %>% group_by(G, time) %>% 
  mutate(BH.Pval=p.adjust(Pval, method="BH")) %>%
  ungroup() %>%
  as.data.frame()

full_adj_res <- full_adj_res %>% group_by(G, time) %>% 
  mutate(BH.Pval=p.adjust(Pval, method="BH")) %>%
  ungroup() %>%
  as.data.frame()

H1$BH.Pval = (full_res %>% filter(H==1))$BH.Pval
H2$BH.Pval = (full_res %>% filter(H==2))$BH.Pval
H3$BH.Pval = (full_res %>% filter(H==3))$BH.Pval
H4$BH.Pval = (full_res %>% filter(H==4))$BH.Pval
HR$BH.Pval = (full_res %>% filter(H==5))$BH.Pval

H1_adj$BH.Pval = (full_adj_res %>% filter(H==1))$BH.Pval
H2_adj$BH.Pval = (full_adj_res %>% filter(H==2))$BH.Pval
H3_adj$BH.Pval = (full_adj_res %>% filter(H==3))$BH.Pval
H4_adj$BH.Pval = (full_adj_res %>% filter(H==4))$BH.Pval
HR_adj$BH.Pval = (full_adj_res %>% filter(H==5))$BH.Pval
ind_adj$BH.Pval = (full_adj_res %>% filter(H==6))$BH.Pval

saveRDS(H1, here("results/unadjusted/H1_res.RDS"))
saveRDS(H2, here("results/unadjusted/H2_res.RDS"))
saveRDS(H3, here("results/unadjusted/H3_res.RDS"))
saveRDS(H4, here("results/unadjusted/H4_res.RDS"))
saveRDS(HR, here("results/unadjusted/HR_res.RDS"))

saveRDS(H1_adj, here("results/adjusted/H1_adj_res.RDS"))
saveRDS(H2_adj, here("results/adjusted/H2_adj_res.RDS"))
saveRDS(H3_adj, here("results/adjusted/H3_adj_res.RDS"))
saveRDS(H4_adj, here("results/adjusted/H4_adj_res.RDS"))
saveRDS(HR_adj, here("results/adjusted/HR_adj_res.RDS"))
saveRDS(ind_adj, here("results/adjusted/individual_adj_res.RDS"))
