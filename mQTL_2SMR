library(devtools)
install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
devtools::install_github("MRCIEU/MRInstruments")
library(MRInstruments)
library(dplyr)
library(R.utils)
library(data.table)

# List available GWASs
ao <- available_outcomes()


for (i in 1:nrow(cpg_list)) {
  ## get mqtl for exposure
  tryCatch(exposure_dat <- format_aries_mqtl(subset(aries_mqtl, cpg==cpg_list[i,1] & age == "Birth")))

  # Extract outcome data
  outcome_dat <- extract_outcome_data(exposure_dat$SNP, outcomes1$id)
  # Set units to SD
  outcome_dat <- subset(outcome_dat, id.outcome %in% outcomes1$id)
  outcomes <- subset(outcomes1, id %in% outcome_dat$id.outcome)
  index <- ! grepl("SD", outcomes$unit, ignore.case=TRUE) & ! grepl("log odds", outcomes$unit, ignore.case=TRUE) & ! is.na(outcomes$sd)
  temp <- outcomes[index,] %>% select(id, sd)
  outcome_dat <- merge(outcome_dat, temp, by.x="id.outcome", by.y="id", all.x=TRUE)
  index <- !is.na(outcome_dat$sd)
  outcome_dat$beta.outcome[index] <- outcome_dat$beta.outcome[index] / outcome_dat$sd[index]
  outcome_dat$se.outcome[index] <- outcome_dat$se.outcome[index] / outcome_dat$sd[index]
    
  # Harmonise data
  tryCatch(dat <- harmonise_data(exposure_dat, outcome_dat), error=function(e) NULL)
    
  # Perform MR
  tryCatch(res <- mr(dat), error=function(e) NULL)
  
  tryCatch(write.table(cbind(i,res), "mQTL_2SMR.txt", append = T, col.names = FALSE, row.names=F), error=function(e) NULL)    
  }
}
