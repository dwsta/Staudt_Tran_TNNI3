source("Utility/contractility_kinetical_loader.R")
library(ggpubr)
library(rstatix)
library(scales)

#paths to data
basepath = "Data/"

well.metrics <- readRDS("Data/InitialControlVariants/CompiledData/per.well.metrics.rds")
normalized.metrics <- readRDS("Data/InitialControlVariants/CompiledData/per.well.metrics.WT.normalized.rds")
combined.traces <- readRDS("Data/InitialControlVariants/CompiledData/combined.traces.rds")


well.metrics <- well.metrics %>% filter(keep_trace == 1) %>%
  mutate(Virus_ordered = factor(Virus,levels = c("TNNI3_WT","TNNI3_R192H","TNNI3_K36Q")))

normalized.metrics <- normalized.metrics %>% filter(keep_trace == 1 | Virus == "GFP") %>%
  mutate(Virus_ordered = factor(str_remove(Virus,"TNNI3_"),levels = c("WT","R192H","K36Q","GFP")))

#############
#Calculate Statistics Using Mixed Model with lmer
############

selected_metrics <- c("mean_valley.contractility",
                      "mean_amplitude.contractility",
                      "pw90_corr.contractility",
                      "pw90_ratio",
                      "pw90_corr.calcium",
                      "mean_peak_value.contractility")

normalized.metrics.long <- normalized.metrics %>%
  filter(keep_trace == 1 | Virus == "GFP") %>%
  pivot_longer(cols = all_of(selected_metrics),
               names_to = "variable",
               values_to = "value") %>%
  select(Virus_ordered,variable,value,plate.ID,well) 

variants.to.analyze <- normalized.metrics %>%
  filter(Virus_ordered != "WT") %>%
  distinct(Virus_ordered) %>% pull()

p_values <- c()
group1 <- c()
group2 <- c()
metric.names <- c()
counter <- 1

for (m in selected_metrics) {
  for (v in variants.to.analyze) {
    df <- normalized.metrics.long %>%
      filter(Virus_ordered %in% c(v,"WT"),
             variable == m)
    
    stat.v <- lmerTest::lmer(formula("value ~ Virus_ordered + (1 | plate.ID)"),
                             data = df)
    p_values[counter] <- summary(stat.v)$coefficients[2,"Pr(>|t|)"]
    group1[counter] <- "WT"
    group2[counter] <- v
    metric.names[counter] <- m
    counter <- counter + 1
  }
}

stat_calc <- data.frame(group1 = group1,group2 = group2,p = p_values,variable = metric.names) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance() %>%
  rename(p.adj.signif.old = p.adj.signif) %>%
  mutate(p.adj.signif = ifelse(p.adj.signif.old == "**","†",
                               ifelse(p.adj.signif.old == "***","#",
                                      ifelse(p.adj.signif.old == "****","§",p.adj.signif.old))))

#save statistics
write.csv(stat_calc,"Files_For_Source_Data/Figure1_stats.csv")

#Panel 1d: Representative Traces

representative_traces <- ggplot(combined.traces %>% filter(plate.ID == "20230526",
                                                           well %in% c("C04","C06","C08","C10"),
                                                           time < 5) %>%
                                  mutate(Virus_ordered = factor(str_remove(Virus,"TNNI3_"),
                                                                levels = c("WT","R192H","K36Q","GFP"))),
                                aes(x = time,y = signal.contractility,color = Virus_ordered)) + 
  geom_line(size = 0.7) + 
  facet_grid(cols = vars(Virus_ordered))+
  theme_pubr()+
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text.x = element_text(size = 18,face = "bold"),
        legend.position = "none")+
  ylim(0,270)+
  xlab("Time (s)") + ylab("Tension (Pa)")
representative_traces

ggsave("Figures/Figure1d.pdf",
       plot = representative_traces,
       width = 10,height = 6)

##########
#Figure Panels 1e,1f,1g,1i,1j,1k
##########



#function to plot each variable
plot_variable_panel <- function(variable.name,ylabel,includeGFP = TRUE) {
  p <- ggplot(normalized.metrics %>% filter(keep_trace == 1 | (Virus == "GFP" & includeGFP == TRUE)),
              aes(x = Virus_ordered,y = .data[[variable.name]])) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(color = plate.ID))+
    stat_pvalue_manual(stat_calc %>% filter(variable == variable.name,group1=="WT",(group2 != "GFP" | includeGFP == TRUE)) %>%
                         left_join(get_y_position(normalized.metrics,as.formula(paste0(variable.name,"~Virus_ordered")))))+
    ylab(ylabel)+
    xlab("")+
    theme_classic()+
    scale_y_continuous(labels = label_number(accuracy = 0.1))+
    theme(axis.text.x = element_text(size = 18,angle = 90),
          axis.text.y = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.position = "none")
  return(p)
}

w <- 3
h <- 6
ggsave("Figures/Figure1e.pdf",
       plot_variable_panel(variable.name="mean_peak_value.contractility",ylabel="Maximum Tension (Normalized)"),
       width = w,height = h)
ggsave("Figures/Figure1f.pdf",
       plot_variable_panel(variable.name="mean_amplitude.contractility",ylabel="Contractile Amplitude (Normalized)"),
       width = w,height = h)
ggsave("Figures/Figure1g.pdf",
       plot_variable_panel(variable.name="mean_valley.contractility",ylabel="Diastolic Tension (Normalized)"),
       width = w,height = h)

ggsave("Figures/Figure1i.pdf",
       plot_variable_panel(variable.name="pw90_corr.contractility",ylabel="Contractile Duration (Normalized)",includeGFP = FALSE),
       width = w,height = h)
ggsave("Figures/Figure1j.pdf",
       plot_variable_panel(variable.name="pw90_corr.calcium",ylabel="Calcium Duration (Normalized)",includeGFP = FALSE),
       width = w,height = h)
ggsave("Figures/Figure1k.pdf",
       plot_variable_panel(variable.name="pw90_ratio",ylabel="Contractile Duration/Calcium Duration (Normalized)",includeGFP = FALSE),
       width = w,height = h)


#Output the above figure data for data reporting summary
figure.data <- normalized.metrics %>% filter(keep_trace == 1 | (Virus == "GFP")) %>%
  select(plate.ID,well,Virus_ordered,
         mean_amplitude.contractility,
         mean_valley.contractility,
         pw90_corr.contractility,
         pw90_ratio,
         pw90_corr.calcium,
         mean_peak_value.contractility)

write.csv(figure.data,"Files_For_Source_Data/Figure1_sourcedata.csv",row.names = FALSE)
