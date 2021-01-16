rm(list=ls())

source(here::here("0-config.R"))

# contains immune-growth covariates and exposures
d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))

# merge in development outcomes
dev <- readRDS(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-development.RDS"))

d_full <- merge(d, dev, "childid")

saveRDS(d_full, paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-development-analysis-dataset.rds"))
