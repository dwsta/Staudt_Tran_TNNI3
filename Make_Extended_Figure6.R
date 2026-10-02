source("Utility/contractility_kinetical_loader.R")
source("Utility/QC_functions.R")
library(data.table)
library(plater)
library(ggprism)
library(ggpubr)
library(ggbeeswarm)
library(ggpmisc)
library(rstatix)
library(lmerTest)
library(patchwork)

platemap <- read_plate("Data/Cal520ToxicityTest/20260707_platemap.csv") %>% rename(well = Wells)

contractility.file <- "Data/Cal520ToxicityTest/contractility_traction_metrics_20260706T160544.csv"
contractility.list <- loadContractility(contractility.file,loadMetrics = T,loadTraces = T)

contractility.metrics <- contractility.list$metrics
contractility.traces <- contractility.list$traces
rm(contractility.list)

contractility.metrics <- contractility.metrics %>%
  parse.well() %>%
  parse.subposition() %>%
  make.plate.names() %>%
  left_join(platemap,by = c("well")) %>%
  mutate(plate_well_sub = paste0(well,subposition))

contractility.traces <- contractility.traces %>%
  parse.well() %>%
  parse.subposition() %>%
  make.plate.names() %>%
  left_join(platemap,by = c("well")) %>%
  mutate(plate_well_sub = paste0(well,subposition))
  

unique(contractility.metrics$plate.name)

#also number the plates
contractility.metrics <- contractility.metrics %>% mutate(corrected_pw90 = mean_pw90/mean_peak_duration,
                                                          corrected_fall_time = mean_fall_time/mean_peak_duration,
                                                          corrected_rise_time = mean_rise_time/mean_peak_duration)


############

pws <- runQC_Traces_singlechannel(contractility.traces,"Data/Cal520ToxicityTest/QC.csv",QC_reset = FALSE)

###########


#####
# load focus QC files
#####

movie.qc <- read.csv("Data/Cal520ToxicityTest/focus_annotation_bluebead_movie.csv",stringsAsFactors = FALSE) %>%
  parse.well(traceID = Label) %>%
  parse.subposition(traceID = Label) %>%
  mutate(keep_movie = Mean != 0) %>%
  select(well,subposition,keep_movie)

focus.qc <- read.csv("Data/Cal520ToxicityTest/focus_annotation_reference_post.csv",stringsAsFactors = FALSE) %>%
  parse.well(traceID = Label) %>%
  parse.subposition(traceID = Label) %>%
  mutate(keep_focus = Mean != 0) %>%
  select(well,subposition,keep_focus)

summarized.metrics <- contractility.metrics %>%
  left_join(pws) %>%
  left_join(movie.qc) %>%
  left_join(focus.qc) %>%
  filter(keep_trace == 1,keep_movie == TRUE,keep_focus == TRUE) %>%
  group_by(well,Virus,Dye) %>%
  summarise(across(where(is.numeric),function(x) {median(x,na.rm = TRUE)}),
            keep_count = n()) %>%
  ungroup() %>%
  filter(keep_count >=3) %>%
  mutate(UID = factor(paste(Virus,Dye),levels = c("WT None",
                                                  "WT Cal520",
                                                  "R192H None",
                                                  "R192H Cal520",
                                                  "K36Q None",
                                                  "K36Q Cal520")))

#plots

#####
stat_test_valley <- 
  t_test(summarized.metrics,
         mean_valley~UID) %>%
  add_y_position()

stat_test_amplitude <- 
  t_test(summarized.metrics,
         mean_amplitude~UID)%>%
  add_y_position()

stat_test_peak <- 
  t_test(summarized.metrics,
         mean_peak_value~UID)%>%
  add_y_position()


stat_test <- bind_rows(stat_test_valley,
                       stat_test_amplitude,
                       stat_test_peak) %>%
  select(-p.adj,-p.adj.signif) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance()
  

stat_test %>% filter()


stat_test_out <- stat_test %>% rename(variable = .y.) %>% select(variable,group1,group2,p,p.adj)
write.csv(stat_test_out,"Files_For_Source_Data/Extended_Figure6_pvalue.csv",row.names = FALSE)



#####

valley.plot.pairs <- 
  ggplot(summarized.metrics,
       aes(x = UID,y = mean_valley,color = Dye)) + 
  geom_violin()+
  geom_beeswarm()+
  theme_pubr()+
  theme()+
  #stat_compare_means(method = "t.test",comparisons = list(c(1,2),c(3,4),c(5,6)))+
  stat_pvalue_manual(data = stat_test %>% filter(.y. == "mean_valley") %>%
                       filter(grepl("K36Q",group1) & grepl("K36Q",group2) |
                                grepl("R192H",group1) & grepl("R192H",group2) |
                                grepl("WT",group1) & grepl("WT",group2)))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+ylab("Diastolic Tension")+
  labs(tag = "A") + 
  theme(plot.tag = element_text(size = 20, face = "bold")) 

valley.plot.pairs


amplitude.plot.pairs <-
  ggplot(summarized.metrics,
         aes(x = UID,y = mean_amplitude,color = Dye)) + 
  geom_violin()+
  geom_beeswarm()+
  theme_pubr()+
  theme()+
  stat_pvalue_manual(data = stat_test %>% filter(.y. == "mean_amplitude") %>%
                       filter(grepl("K36Q",group1) & grepl("K36Q",group2) |
                                grepl("R192H",group1) & grepl("R192H",group2) |
                                grepl("WT",group1) & grepl("WT",group2)))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+ylab("Generated Force (Amplitude)")+
  labs(tag = "D") + 
  theme(plot.tag = element_text(size = 20, face = "bold")) 

peak.plot.pairs <-
  ggplot(summarized.metrics,
         aes(x = UID,y = mean_peak_value,color = Dye)) + 
  geom_violin()+
  geom_beeswarm()+
  theme_pubr()+
  theme()+
  stat_pvalue_manual(data = stat_test %>% filter(.y. == "mean_peak_value") %>%
                       filter(grepl("K36Q",group1) & grepl("K36Q",group2) |
                                grepl("R192H",group1) & grepl("R192H",group2) |
                                grepl("WT",group1) & grepl("WT",group2)))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+ylab("Maximum Force (Diastolic + Systolic)")+
  labs(tag = "G") + 
  theme(plot.tag = element_text(size = 20, face = "bold")) 


pair.plot <- valley.plot.pairs / amplitude.plot.pairs / peak.plot.pairs
pair.plot

#########

valley.plot.dye <- 
  ggplot(summarized.metrics %>% filter(Dye == "Cal520"),
         aes(x = UID,y = mean_valley,color = Dye)) + 
  geom_violin()+
  geom_beeswarm()+
  theme_pubr()+
  theme()+
  stat_pvalue_manual(data = stat_test %>% filter(.y. == "mean_valley",grepl("Cal520",group1),grepl("Cal520",group2),grepl("WT",group1) | grepl("WT",group2)))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+ylab("Diastolic Tension")+
  labs(tag = "C") + 
  theme(plot.tag = element_text(size = 20, face = "bold")) 

amplitude.plot.dye <-
  ggplot(summarized.metrics %>% filter(Dye == "Cal520"),
         aes(x = UID,y = mean_amplitude,color = Dye)) + 
  geom_violin()+
  geom_beeswarm()+
  theme_pubr()+
  theme()+
  stat_pvalue_manual(data = stat_test %>% filter(.y. == "mean_amplitude",grepl("Cal520",group1),grepl("Cal520",group2),grepl("WT",group1) | grepl("WT",group2)))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+ylab("Generated Force (Amplitude)")+
  labs(tag = "F") + 
  theme(plot.tag = element_text(size = 20, face = "bold")) 

peak.plot.dye <-
  ggplot(summarized.metrics %>% filter(Dye == "Cal520"),
         aes(x = UID,y = mean_peak_value,color = Dye)) + 
  geom_violin()+
  geom_beeswarm()+
  theme_pubr()+
  theme()+
  stat_pvalue_manual(data = stat_test %>% filter(.y. == "mean_peak_value",grepl("Cal520",group1),grepl("Cal520",group2),grepl("WT",group1) | grepl("WT",group2)))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+ylab("Maximum Force (Diastolic + Systolic)")+
  labs(tag = "I") + 
  theme(plot.tag = element_text(size = 20, face = "bold")) 


###########

valley.plot.none <- 
  ggplot(summarized.metrics %>% filter(Dye == "None"),
         aes(x = UID,y = mean_valley,color = Dye)) + 
  geom_violin()+
  geom_beeswarm()+
  theme_pubr()+
  theme()+
  stat_pvalue_manual(data = stat_test %>% filter(.y. == "mean_valley",grepl("None",group1),grepl("None",group2),grepl("WT",group1) | grepl("WT",group2)))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+ylab("Diastolic Tension")+
  labs(tag = "B") + 
  theme(plot.tag = element_text(size = 20, face = "bold")) 

amplitude.plot.none <-
  ggplot(summarized.metrics %>% filter(Dye == "None"),
         aes(x = UID,y = mean_amplitude,color = Dye)) + 
  geom_violin()+
  geom_beeswarm()+
  theme_pubr()+
  theme()+
  stat_pvalue_manual(data = stat_test %>% filter(.y. == "mean_amplitude",grepl("None",group1),grepl("None",group2),grepl("WT",group1) | grepl("WT",group2)))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+ylab("Generated Force (Amplitude)")+
  labs(tag = "E") + 
  theme(plot.tag = element_text(size = 20, face = "bold")) 

peak.plot.none <-
  ggplot(summarized.metrics %>% filter(Dye == "None"),
         aes(x = UID,y = mean_peak_value,color = Dye)) + 
  geom_violin()+
  geom_beeswarm()+
  theme_pubr()+
  theme()+
  stat_pvalue_manual(data = stat_test %>% filter(.y. == "mean_peak_value",grepl("None",group1),grepl("None",group2),grepl("WT",group1) | grepl("WT",group2)))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+ylab("Maximum Force (Diastolic + Systolic)")+
  labs(tag = "H") + 
  theme(plot.tag = element_text(size = 20, face = "bold")) 



#(valley.plot.pairs+valley.plot.none+valley.plot.dye)/(amplitude.plot.pairs+amplitude.plot.none+amplitude.plot.dye)/(peak.plot.pairs+peak.plot.none+peak.plot.dye)


combined <-
  valley.plot.pairs    + valley.plot.none    + valley.plot.dye +
  amplitude.plot.pairs + amplitude.plot.none + amplitude.plot.dye +
  peak.plot.pairs      + peak.plot.none      + peak.plot.dye +
  plot_layout(ncol = 3, widths = c(2, 1, 1))

combined

ggsave("Figures/Extended_Data_Figure6.pdf",combined,width = 13,height = 18)

data.for.output <- summarized.metrics %>%
  select(UID,mean_valley,mean_amplitude,mean_peak_value)

write.csv(data.for.output,"Files_For_Source_Data/Extended_Figure6_data.csv",row.names = FALSE)
