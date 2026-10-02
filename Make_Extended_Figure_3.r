#mcherry data for paper
library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggpmisc)

normalized.metrics.mcherry <- readRDS("Data/PathogenicScreen/CompiledData/per.well.metrics.WT.normalized.with.mcherry.rds")

#load lookup table
virus.ID.lookup <- read.csv("Data/PathogenicScreen/Lookup/TNNI3screen_pathvariants_lookuptable.csv")

#join the lookup table to the normalized data
normalized.metrics.mcherry <- normalized.metrics.mcherry %>%
  left_join(virus.ID.lookup %>% select(-X)) %>%
  mutate(Clinical.significance = ifelse(Clinical.significance == "Other","Pathogenic",Clinical.significance)) # rename one variant

#select out only TNNI3 variants (and not phophomimetics)
normalized.metrics.TNNI3.mcherry <- normalized.metrics.mcherry %>% filter(Gene == "TNNI3",!(Clinical.significance %in% c(NA,"Phosphomimetic")))

mean_amplitude_mcherry_plot <- 
  ggplot(normalized.metrics.TNNI3.mcherry,
       aes(x = mean.mcherry.fluorescence,y = mean_amplitude.contractility))+
  geom_point()+
  stat_poly_line()+
  stat_poly_eq(use_label("rr"))+
  xlab("mCherry Fluorescence (Normalized)")+
  ylab("Contractile Amplitude (Normalized)")+
  theme(plot.tag.position = c(0.02,0.98),plot.tag = element_text(face = "bold",size = 20))

mean_valley_mcherry_plot <- 
  ggplot(normalized.metrics.TNNI3.mcherry,
         aes(x = mean.mcherry.fluorescence,y = mean_valley.contractility))+
  geom_point()+
  stat_poly_line()+
  stat_poly_eq(use_label("rr"))+
  xlab("mCherry Fluorescence (Normalized)")+
  ylab("Diastolic Tension (Normalized)")+
  theme(plot.tag.position = c(0.02,0.98),plot.tag = element_text(face = "bold",size = 20))


mean_peak_value_mcherry_plot <- 
  ggplot(normalized.metrics.TNNI3.mcherry,
         aes(x = mean.mcherry.fluorescence,y = mean_peak_value.contractility))+
  geom_point()+
  stat_poly_line()+
  stat_poly_eq(use_label("rr"))+
  xlab("mCherry Fluorescence (Normalized)")+
  ylab("Maximum Tension (Normalized)")+
  theme(plot.tag.position = c(0.02,0.98),plot.tag = element_text(face = "bold",size = 20))

pw90_ratio_mcherry_plot <- 
  ggplot(normalized.metrics.TNNI3.mcherry,
         aes(x = mean.mcherry.fluorescence,y = pw90_ratio))+
  geom_point()+
  stat_poly_line()+
  stat_poly_eq(use_label("rr"))+
  xlab("mCherry Fluorescence (Normalized)")+
  ylab("Contractile Duration/Calcium Duration (Normalized)")+
  theme(plot.tag.position = c(0.02,0.98),plot.tag = element_text(face = "bold",size = 20))

fall_ratio_mcherry_plot <- 
  ggplot(normalized.metrics.TNNI3.mcherry,
         aes(x = mean.mcherry.fluorescence,y = fall_ratio))+
  geom_point()+
  stat_poly_line()+
  stat_poly_eq(use_label("rr"))+
  xlab("mCherry Fluorescence (Normalized)")+
  ylab("Relaxation Time/Calcium Fall Time (Normalized)")+
  theme(plot.tag.position = c(0.02,0.98),plot.tag = element_text(face = "bold",size = 20))

proportion.valley.contractility_mcherry_plot <- 
  ggplot(normalized.metrics.TNNI3.mcherry,
         aes(x = mean.mcherry.fluorescence,y = proportion.valley.contractility))+
  geom_point()+
  stat_poly_line()+
  stat_poly_eq(use_label("rr"))+
  xlab("mCherry Fluorescence (Normalized)")+
  ylab("Diastolic Tension/Maximum Tension (Normalized)")+
  theme(plot.tag.position = c(0.02,0.98),plot.tag = element_text(face = "bold",size = 20))

fall_rate.contractility_mcherry_plot <- 
  ggplot(normalized.metrics.TNNI3.mcherry,
         aes(x = mean.mcherry.fluorescence,y = fall_rate.contractility))+
  geom_point()+
  stat_poly_line()+
  stat_poly_eq(use_label("rr"))+
  xlab("mCherry Fluorescence (Normalized)")+
  ylab("Relaxation Rate (Normalized)")+
  theme(plot.tag.position = c(0.02,0.98),plot.tag = element_text(face = "bold",size = 20))

rise_rate.calcium_mcherry_plot <- 
  ggplot(normalized.metrics.TNNI3.mcherry,
         aes(x = mean.mcherry.fluorescence,y = rise_rate.calcium))+
  geom_point()+
  stat_poly_line()+
  stat_poly_eq(use_label("rr"))+
  xlab("mCherry Fluorescence (Normalized)")+
  ylab("Calcium Upstroke Time (Normalized)")+
  theme(plot.tag.position = c(0.02,0.98),plot.tag = element_text(face = "bold",size = 20))

mcherry_plots <-
  ggarrange(mean_amplitude_mcherry_plot+labs(tag = "a"),mean_valley_mcherry_plot+labs(tag = "b"),mean_peak_value_mcherry_plot+labs(tag = "c"),
          pw90_ratio_mcherry_plot+labs(tag = "d"),fall_ratio_mcherry_plot+labs(tag = "e"),
          proportion.valley.contractility_mcherry_plot+labs(tag = "f"),fall_rate.contractility_mcherry_plot+labs(tag = "g"),
          rise_rate.calcium_mcherry_plot+labs(tag = "h"))

ggsave("Figures/Extended_Data_Figure_3.pdf",mcherry_plots,width = 13,height = 13)

mcherry.data <- normalized.metrics.TNNI3.mcherry %>%
  select(plate.ID,plate.letter,well,
         mean_amplitude.contractility,
         mean_valley.contractility,
         mean_peak_value.contractility,
         pw90_ratio,
         fall_ratio,
         proportion.valley.contractility,
         fall_rate.contractility,
         rise_rate.calcium,
         mean.mcherry.fluorescence,
         Variant
  )


write.csv(mcherry.data,"Files_For_Source_Data/Extended_Figure3_data.csv",row.names = FALSE)
