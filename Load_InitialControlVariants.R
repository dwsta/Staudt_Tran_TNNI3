source("Utility/contractility_kinetical_loader.R")
source("Utility/QC_functions.R")
library(plater)

###########
#load platemaps
###########
platemap1 <- read_plate("Data/InitialControlVariants/Platemaps/platemap96_PT230515.csv") %>% rename(well = Wells)
platemap2 <- read_plate("Data/InitialControlVariants/Platemaps/platemap96_PT230519.csv") %>% rename(well = Wells)
platemap3 <- read_plate("Data/InitialControlVariants/Platemaps/platemap96_PT230714.csv") %>% rename(well = Wells)

#combine and label with readout day
platemaps <-bind_rows(platemap1 %>% mutate(plate.ID = "20230523"),
                      platemap2 %>% mutate(plate.ID = "20230526"),
                      platemap3 %>% mutate(plate.ID = "20230722"))

###########
#load contractility
###########

contractility.files <- c("Data/InitialControlVariants/Contractility/contractility_traction_20230525T175922_PT230515.csv",
                         "Data/InitialControlVariants/Contractility/contractility_traction_20230528T180150_PT230519.csv",
                         "Data/InitialControlVariants/Contractility/contractility_traction_20230724T075529_PT230714.csv")

contractility.list <- loadContractility(contractility.files,loadMetrics = T,loadTraces = T)
contractility.metrics <- contractility.list$metrics
contractility.traces <- contractility.list$traces
rm(contractility.list)

contractility.metrics <- contractility.metrics %>%
  parse.well() %>%
  make.plate.names.nochannel() %>%
  mutate(pw90_corr = mean_pw90/mean_peak_duration,
         rise_corr = mean_rise_time/mean_peak_duration,
         fall_corr = mean_fall_time/mean_peak_duration,
         rise_rate = mean_amplitude/mean_rise_time,
         fall_rate = mean_amplitude/mean_fall_time,
         rate_ratio = mean_rise_time/mean_fall_time,
         rise_fall_ratio = mean_rise_time/mean_fall_time,
         plate.ID = str_extract(plate.name,"_(\\d{8})",group = 1))

contractility.traces <- contractility.traces %>%
  parse.well() %>%
  make.plate.names.nochannel() %>%
  mutate(plate.ID = str_extract(plate.name,"_(\\d{8})",group = 1)) %>%
  filter(!grepl("iso",plate.name)) #get rid of iso treated data in this figure


###########
#load calcium
###########

calcium.files <- c("Data/InitialControlVariants/Calcium/calcium_whole_roi_PT230515_PT230519.csv",
                   "Data/InitialControlVariants/Calcium/calcium_whole_roi_PT230714.csv")

calcium.list <- loadContractility(calcium.files,loadMetrics = T,loadTraces = T)
calcium.metrics <- calcium.list$metrics
calcium.traces <- calcium.list$traces
rm(calcium.list)

calcium.metrics <- calcium.metrics %>%
  parse.well() %>%
  make.plate.names.nochannel() %>%
  mutate(pw90_corr = mean_pw90/mean_peak_duration,
         rise_corr = mean_rise_time/mean_peak_duration,
         fall_corr = mean_fall_time/mean_peak_duration,
         rise_rate = mean_amplitude/mean_rise_time,
         fall_rate = mean_amplitude/mean_fall_time,
         rate_ratio = mean_rise_time/mean_fall_time,
         rise_fall_ratio = mean_rise_time/mean_fall_time,
         deltaF_over_baseline = mean_amplitude/mean_baseline,
         plate.ID = str_extract(plate.name,"_(\\d{8})",group = 1))

calcium.traces <- calcium.traces %>%
  parse.well() %>%
  #parse.subposition() %>%
  make.plate.names.nochannel() %>%
  mutate(plate.ID = str_extract(plate.name,"_(\\d{8})",group = 1)) 

#extract interleaved channel from calcium to merge shortly
matched.calcium.metrics <- calcium.metrics %>% filter(grepl("Cal520",alias))
matched.calcium.traces <- calcium.traces %>% filter(grepl("Cal520",alias))

#############
# merge the two datasets
#############

combined.metrics <- left_join(contractility.metrics %>% select(-alias), #clean the data a little
                              matched.calcium.metrics %>% select(-alias,-well.letter,-well.number,-plate.name),
                              by = c("plate.ID","well"),
                              suffix = c(".contractility",".calcium")) %>%
  mutate(pw90_ratio = mean_pw90.contractility/mean_pw90.calcium) %>%
  mutate(plate_well_sub = paste(plate.ID,well,sep = "_"))


combined.traces <- left_join(contractility.traces %>% select(-alias),
                             matched.calcium.traces %>% select(-alias,-well.letter,-well.number,-plate.name),
                             by = c("plate.ID","well","time"),
                             suffix = c(".contractility",".calcium")) %>%
  mutate(plate_well_sub = paste(plate.ID,well,sep = "_"))


#############
# QC #
#############

# -1 means needs to be read
# 0 means do not keep (flatlines or close to)
# 1 means keep (mostly rhythmic)
# 2 means "maybe" (arrhythmias but some activity)
# -2 means error

pws <- runQC_Traces(combined.traces,"Data/InitialControlVariants/QC/QC.csv",QC_reset = FALSE)

combined.metrics <- combined.metrics %>% left_join(pws)

combined.traces <- combined.traces %>% left_join(pws)

#well metrics will take median of each well
#then join to everything
well.metrics <- combined.metrics %>% 
  group_by(plate.ID,well,well.letter,well.number) %>%
  summarise(across(where(is.numeric),
            ~median(.x,na.rm = TRUE))) %>%
  ungroup() %>%
  left_join(platemaps)
saveRDS(well.metrics,"Data/InitialControlVariants/CompiledData/per.well.metrics.rds")

# #normalize to wild type

WT_values <- well.metrics %>%
  filter(Virus == "TNNI3_WT") %>%
  group_by(plate.ID) %>%
  summarise(across(where(is.numeric),~median(.x,na.rm = TRUE),.names = "plateWT_median.{col}")) %>%
  ungroup()

well.metrics.plateNorm <- well.metrics %>% left_join(WT_values) %>%
  mutate(across(where(is.numeric) & !starts_with("plateWT_median"),
                ~ .x / get(paste0("plateWT_median.", cur_column())))) %>%
  select(-starts_with("plateWT_median"))

saveRDS(well.metrics.plateNorm,"Data/InitialControlVariants/CompiledData/per.well.metrics.WT.normalized.rds")

#label the combined traces
combined.traces <- combined.traces %>% 
  left_join(platemaps) 

saveRDS(combined.traces,"Data/InitialControlVariants/CompiledData/combined.traces.rds")
