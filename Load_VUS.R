source("Utility/contractility_kinetical_loader.R")
source("Utility/QC_functions.R")
library(plater)


###########
#load platemaps
###########
platemap <- read_plate("Data/VUS/Platemaps/platemap_VUS.csv") %>% rename(well = Wells)

###########
#load contractility
###########

contractility.files <- "Data/VUS/Contractility/contractility_traction_20241009T114413.csv"

contractility.list <- loadContractility(contractility.files,loadMetrics = T,loadTraces = T)
contractility.metrics <- contractility.list$metrics
contractility.traces <- contractility.list$traces
rm(contractility.list)

contractility.metrics <- contractility.metrics %>%
  parse.well() %>%
  parse.subposition() %>%
  make.plate.names.nochannel() %>%
  mutate(pw90_corr = mean_pw90/mean_peak_duration,
         rise_corr = mean_rise_time/mean_peak_duration,
         fall_corr = mean_fall_time/mean_peak_duration,
         rise_rate = mean_amplitude/mean_rise_time,
         fall_rate = mean_amplitude/mean_fall_time,
         rate_ratio = mean_rise_time/mean_fall_time,
         rise_fall_ratio = mean_rise_time/mean_fall_time) %>%
  mutate(plate.ID = str_extract(plate.name,"plate(\\d{8}[:alpha:]?)",group = 1))

contractility.traces <- contractility.traces %>%
  parse.well() %>%
  parse.subposition() %>%
  make.plate.names.nochannel() %>%
  mutate(plate.ID = str_extract(plate.name,"plate(\\d{8}[:alpha:]?)",group = 1))


###########
#load calcium
###########

calcium.files <- "Data/VUS/Calcium/calcium_whole_roi_20241011T100324.csv"

calcium.list <- loadContractility(calcium.files,loadMetrics = T,loadTraces = T)
calcium.metrics <- calcium.list$metrics
calcium.traces <- calcium.list$traces
rm(calcium.list)

calcium.metrics <- calcium.metrics %>%
  parse.well() %>%
  parse.subposition() %>%
  make.plate.names.nochannel() %>%
  mutate(pw90_corr = mean_pw90/mean_peak_duration,
         rise_corr = mean_rise_time/mean_peak_duration,
         fall_corr = mean_fall_time/mean_peak_duration,
         rise_rate = mean_amplitude/mean_rise_time,
         fall_rate = mean_amplitude/mean_fall_time,
         rate_ratio = mean_rise_time/mean_fall_time,
         rise_fall_ratio = mean_rise_time/mean_fall_time,
         deltaF_over_baseline = mean_amplitude/mean_baseline) %>%
  mutate(plate.ID = str_extract(plate.name,"plate(\\d{8}[:alpha:]?)",group = 1)) %>%
  mutate(channel = ifelse(grepl("_CyanEx",alias),"CyanEx","FITC"))

calcium.traces <- calcium.traces %>%
  parse.well() %>%
  parse.subposition() %>%
  make.plate.names.nochannel() %>%
  mutate(plate.ID = str_extract(plate.name,"plate(\\d{8}[:alpha:]?)",group = 1)) %>%
  mutate(channel = ifelse(grepl("_CyanEx",alias),"CyanEx","FITC"))


#extract interleaved channel from calcium to merge shortly
matched.calcium.metrics <- calcium.metrics %>% filter(channel == "CyanEx")
matched.calcium.traces <- calcium.traces %>% filter(channel == "CyanEx")

#############
# merge the two datasets
#############

combined.metrics <- left_join(contractility.metrics %>% select(-alias), #clean the data a little
                              matched.calcium.metrics %>% select(-alias,-well.letter,-well.number,-subposition.row,-subposition.col),
                              by = c("plate.ID","well","subposition"),
                              suffix = c(".contractility",".calcium")) %>%
  mutate(pw90_ratio = mean_pw90.contractility/mean_pw90.calcium) %>%
  mutate(plate_well_sub = paste(plate.ID,well,subposition,sep = "_"))
  

combined.traces <- left_join(contractility.traces %>% select(-alias),
                             matched.calcium.traces %>% select(-alias,-well.letter,-well.number,-subposition.row,-subposition.col),
                             by = c("plate.ID","well","subposition","time"),
                             suffix = c(".contractility",".calcium")) %>%
  mutate(plate_well_sub = paste(plate.ID,well,subposition,sep = "_"))
  

#############
# QC #
#############

# -1 means needs to be read
# 0 means do not keep (flatlines or close to)
# 1 means keep (mostly rhythmic)
# 2 means "maybe" (arrhythmias but some activity)
# -2 means error

pws <- runQC_Traces(combined.traces,"Data/VUS/QC/QC.csv",QC_reset = FALSE)

combined.metrics <- combined.metrics %>% left_join(pws)

combined.traces <- combined.traces %>% left_join(pws)

#well metrics will take median of each well
#then join to everything
well.metrics <- combined.metrics %>% 
  filter(keep_trace == 1) %>%
  group_by(plate.ID,well,well.letter,well.number) %>%
  summarise(across(where(is.numeric),
            ~median(.x,na.rm = TRUE))) %>%
  ungroup() %>%
  left_join(platemap)
saveRDS(well.metrics,"Data/VUS/CompiledData/per.well.metrics.rds")

# #normalize to wild type

WT_values <- well.metrics %>%
  filter(Virus == "WT") %>%
  group_by(plate.ID) %>%
  summarise(across(where(is.numeric),~median(.x,na.rm = TRUE),.names = "plateWT_median.{col}")) %>%
  ungroup()

well.metrics.plateNorm <- well.metrics %>% left_join(WT_values) %>%
  mutate(across(where(is.numeric) & !starts_with("plateWT_median"),
                ~ .x / get(paste0("plateWT_median.", cur_column())))) %>%
  select(-starts_with("plateWT_median"))

saveRDS(well.metrics.plateNorm,"Data/VUS/CompiledData/per.well.metrics.WT.normalized.rds")

#label the combined traces
combined.traces <- combined.traces %>% 
  left_join(platemap) 

saveRDS(combined.traces,"Data/VUS/CompiledData/combined.traces.rds")
