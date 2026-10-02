source("Utility/contractility_kinetical_loader.R")
source("Utility/QC_functions.R")
library(plater)


###########
#load platemaps
###########
platemap1 <- read_plate("Data/PathogenicScreen/Platemaps/platemap96_PT_20240618.csv") %>% rename(well = Wells)
platemap2 <- read_plate("Data/PathogenicScreen/Platemaps/platemap96_PT_20240621.csv") %>% rename(well = Wells)
platemap3 <- read_plate("Data/PathogenicScreen/Platemaps/platemap96_PT_20240624.csv") %>% rename(well = Wells)

#combine and label with readout day
platemaps <-bind_rows(platemap1 %>% mutate(plate.ID = "20240618"),
                      platemap2 %>% mutate(plate.ID = "20240621"),
                      platemap3 %>% mutate(plate.ID = "20240624"))


###########
#load contractility
###########

contractility.files <- c("Data/PathogenicScreen/Contractility/contractility_traction_20240730T110043.csv",
                         "Data/PathogenicScreen/Contractility/contractility_traction_20240720T123707.csv")

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
         rise_fall_ratio = mean_rise_time/mean_fall_time,
         proportion.valley = mean_valley/mean_peak_value,
         plate.ID = str_extract(alias,"PT_(\\d{8})",group = 1),
         plate.letter = str_extract(alias,"plate([A-Z])",group = 1))

contractility.traces <- contractility.traces %>%
  parse.well() %>%
  parse.subposition() %>%
  make.plate.names.nochannel() %>%
  mutate(plate.ID = str_extract(alias,"PT_(\\d{8})",group = 1),
         plate.letter = str_extract(alias,"plate([A-Z])",group = 1))

# remove duplicated row in plateC from repeat 2
# (cal520_bluebeads_PT_20240621_10s_20240621180454_plateC)
# this was rerun due to no beating
contractility.metrics <- contractility.metrics %>%
  filter(!(plate.name == "cal520_bluebeads_PT_20240621_10s_20240621180454_plateC" & well.letter == "B"))

contractility.traces <- contractility.traces %>%
  filter(!(plate.name == "cal520_bluebeads_PT_20240621_10s_20240621180454_plateC" & well.letter == "B"))

###########
#load calcium
###########

calcium.files <- c("Data/PathogenicScreen/Calcium/calcium_whole_roi_20240704T103710.csv")

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
         rise_fall_ratio = mean_rise_time/mean_fall_time,
         deltaF_over_F0 = mean_amplitude/mean_valley,
         proportion.valley = mean_valley/mean_peak_value, 
         plate.ID = str_extract(alias,"PT_(\\d{8})",group = 1),
         channel = str_extract(alias,"c_\\d{4}_([[:alpha:]]+)",group = 1),
         plate.letter = str_extract(alias,"plate([A-Z])",group = 1))

calcium.traces <- calcium.traces %>%
  parse.well() %>%
  parse.subposition() %>%
  make.plate.names.nochannel() %>%
  mutate(plate.ID = str_extract(alias,"PT_(\\d{8})",group = 1),
         channel = str_extract(alias,"c_\\d{4}_([[:alpha:]]+)",group = 1),
         plate.letter = str_extract(alias,"plate([A-Z])",group = 1))

# remove duplicated row in plateC from repeat 2
# (cal520_bluebeads_PT_20240621_10s_20240621180454_plateC)
# this was rerun due to no beating

calcium.metrics <- calcium.metrics %>%
  filter(!(plate.name == "cal520_bluebeads_PT_20240621_10s_20240621180454_plateC" & well.letter == "B"))

calcium.traces <- calcium.traces %>%
  filter(!(plate.name == "cal520_bluebeads_PT_20240621_10s_20240621180454_plateC" & well.letter == "B"))


#extract interleaved channel from calcium to merge shortly
matched.calcium.metrics <- calcium.metrics %>% filter(channel == "CyanEx")
matched.calcium.traces <- calcium.traces %>% filter(channel == "CyanEx")

F0.df <- calcium.traces %>% filter(channel == "FITC") %>%
  group_by(plate.ID,plate.letter,well,subposition) %>%
  summarise(F0 = min(signal))

#####
#load mcherry fluorescence
#####

mcherry.data <- read.csv("Data/PathogenicScreen/mcherry/mcherry.fluorescence.csv",stringsAsFactors = FALSE) %>%
  mutate(plate.ID = as.character(plate.ID))

#############
# merge the two datasets
#############

combined.metrics <- left_join(contractility.metrics %>% select(-alias), #clean the data a little
                              matched.calcium.metrics %>% select(-alias,-well.letter,-well.number,-subposition.row,-subposition.col),
                              by = c("plate.ID","plate.letter","well","subposition"),
                              suffix = c(".contractility",".calcium")) %>%
  mutate(pw90_ratio = mean_pw90.contractility/mean_pw90.calcium,
         rise_ratio = mean_rise_time.contractility/mean_rise_time.calcium,
         fall_ratio = mean_fall_time.contractility/mean_fall_time.calcium,
         ratio_of_fall_to_rise_ratios = fall_ratio/rise_ratio) %>%
  mutate(plate_well_sub = paste(plate.ID,plate.letter,well,subposition,sep = "_")) %>%
  left_join(F0.df) %>%
  left_join(mcherry.data)


combined.traces <- left_join(contractility.traces %>% select(-alias),
                             matched.calcium.traces %>% select(-alias,-well.letter,-well.number,-subposition.row,-subposition.col),
                             by = c("plate.ID","plate.letter","well","subposition","time"),
                             suffix = c(".contractility",".calcium")) %>%
  mutate(plate_well_sub = paste(plate.ID,plate.letter,well,subposition,sep = "_"))


#############
# QC #
#############

# -1 means needs to be read
# 0 means do not keep (flatlines or close to)
# 1 means keep (mostly rhythmic)
# 2 means "maybe" (arrhythmias but some activity)
# -2 means error

pws <- runQC_Traces(combined.traces,"Data/PathogenicScreen/QC/QC.csv",QC_reset = FALSE)

combined.metrics <- combined.metrics %>% left_join(pws)

combined.traces <- combined.traces %>% left_join(pws)

#additional, secondary QC by another author
ppqtran.qc <- read.csv("Data/PathogenicScreen/QC/QC_PPQTRAN.csv",stringsAsFactors = FALSE) %>%
  rename(keep_trace_ppqtran = keep_trace,
         alias = plate_well_sub) %>%
  mutate(plate.ID = str_extract(alias,"PT_(\\d{8})",group = 1),
         plate.letter = str_extract(alias,"_plate([A-Z])",group = 1),
         well = str_extract(alias,"_([A-Z]\\d{2})_",group = 1),
         subposition = str_extract(alias,"_(r\\d{2}c\\d{2})",group = 1),
         plate_well_sub = paste(plate.ID,plate.letter,well,subposition,sep = "_")) %>%
  select(plate_well_sub,keep_trace_ppqtran)
  
combined.metrics <- combined.metrics %>% left_join(ppqtran.qc)

combined.traces <- combined.traces %>% left_join(ppqtran.qc)

#well metrics will take median of each well
#then join to everything
well.metrics <- combined.metrics %>% 
  filter(keep_trace == 1,keep_trace_ppqtran == 1) %>%
  group_by(plate.ID,plate.letter,well,well.letter,well.number) %>%
  summarise(across(where(is.numeric),
            ~median(.x,na.rm = TRUE))) %>%
  ungroup() %>%
  left_join(platemaps)

saveRDS(well.metrics %>% select(-mean.mcherry.fluorescence),"Data/PathogenicScreen/CompiledData/per.well.metrics.rds") # main analysis file
saveRDS(well.metrics,"Data/PathogenicScreen/CompiledData/per.well.metrics.with.mcherry.rds") # extra with mcherry data


# #normalize to wild type

WT_values <- well.metrics %>%
  filter(Virus.ID == "3WT.7") %>%
  group_by(plate.ID,plate.letter) %>%
  summarise(across(where(is.numeric),~median(.x,na.rm = TRUE),.names = "plateWT_median.{col}")) %>%
  ungroup()

well.metrics.plateNorm <- well.metrics %>% left_join(WT_values) %>%
  mutate(across(where(is.numeric) & !starts_with("plateWT_median"),
                ~ .x / get(paste0("plateWT_median.", cur_column())))) %>%
  select(-starts_with("plateWT_median"))

saveRDS(well.metrics.plateNorm %>% select(-mean.mcherry.fluorescence),"Data/PathogenicScreen/CompiledData/per.well.metrics.WT.normalized.rds")
saveRDS(well.metrics.plateNorm,"Data/PathogenicScreen/CompiledData/per.well.metrics.WT.normalized.with.mcherry.rds")

#label the combined traces
combined.traces <- combined.traces %>% 
  left_join(platemaps) 

saveRDS(combined.traces,"Data/PathogenicScreen/CompiledData/combined.traces.rds")
