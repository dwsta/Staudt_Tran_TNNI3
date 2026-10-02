source("Utility/contractility_kinetical_loader.R")
library(ggpubr)
library(rstatix)
library(scales)

virus.variant.lookup <- read.csv("Data/VUS/Lookup/VUS.variant.lookup.csv",stringsAsFactors = FALSE)

well.metrics <- readRDS("Data/VUS/CompiledData/per.well.metrics.rds") %>% 
  left_join(virus.variant.lookup)
normalized.metrics <- readRDS("Data/VUS/CompiledData/per.well.metrics.WT.normalized.rds") %>%
  left_join(virus.variant.lookup)
combined.traces <- readRDS("Data/VUS/CompiledData/combined.traces.rds") %>%
  left_join(virus.variant.lookup)

ordering.pw90_ratio.df <- normalized.metrics %>%
  select(Variant,pw90_ratio) %>% 
  group_by(Variant) %>%
  summarise(pw90_ratio = median(pw90_ratio,na.rm=TRUE)) %>%
  arrange(desc(pw90_ratio))


well.metrics <- well.metrics %>% #filter(keep_trace == 1) %>%
  #mutate(Variant_ordered = factor(Variant,levels = c("TNNI3_WT","TNNI3_R192H","TNNI3_K36Q")))
  mutate(Variant_ordered = factor(Variant,levels = ordering.pw90_ratio.df$Variant))

normalized.metrics <- normalized.metrics %>% #filter(keep_trace == 1) %>%
  #mutate(Variant_ordered = factor(str_remove(Variant,"TNNI3_"),levels = c("WT","R192H","K36Q")))
  mutate(Variant_ordered = factor(Variant,levels = ordering.pw90_ratio.df$Variant))

combined.traces <- readRDS("Data/VUS/CompiledData/combined.traces.rds") %>%
  left_join(virus.variant.lookup) %>%
  mutate(Variant_ordered = factor(Variant,levels = ordering.pw90_ratio.df$Variant))

#statistics

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
  select(Variant_ordered,variable,value,plate.ID,well) 

variants.to.analyze <- normalized.metrics %>%
  filter(Variant_ordered != "WT") %>%
  distinct(Variant_ordered) %>% pull()

p_values <- c()
group1 <- c()
group2 <- c()
metric.names <- c()
counter <- 1

for (m in selected_metrics) {
  for (v in variants.to.analyze) {
    df <- normalized.metrics.long %>%
      filter(Variant_ordered %in% c(v,"WT"),
             variable == m)
    
    stat.v <- lmerTest::lmer(formula("value ~ Variant_ordered + (1 | plate.ID)"),
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

write.csv(stat_calc,"Files_For_Source_Data/Figure3_pvalues.csv",row.names = FALSE)

#Manual Subsets

#FROM pw90_ratio data NEW
HCM.subset <- c("RR145-146dup",
                "R69S",
                "WT")
DCM.subset <- c("V104M",
                "Q48P",
                "WT")

nothing.subset <- c("I56M",
                    "R103H",
                    "A116G",
                    "R103C",
                    "A161V",
                    "WT")



#function to plot each variable
plot_variable_panel <- function(variable.name,ylabel,virus.subset) {
  p <- ggplot(normalized.metrics %>% filter(keep_trace == 1,Variant_ordered %in% virus.subset),
         aes(x = Variant_ordered,y = .data[[variable.name]],
             color = ifelse(Variant_ordered %in% nothing.subset,"WT/Benign",
                            ifelse(Variant_ordered %in% HCM.subset,"HCM",
                                   ifelse(Variant_ordered %in% DCM.subset,"DCM","NA"))))) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter()+
    #geom_jitter(aes(color = plate.ID))+
    # stat_pvalue_manual(stat_calc %>% filter(variable == variable.name) %>%
    #                      left_join(get_y_position(normalized.metrics,as.formula(paste0(variable.name,"~Variant_ordered")),ref.group = "WT")),
    #                    remove.bracket = TRUE)+
    stat_pvalue_manual(stat_calc %>% filter(variable == variable.name,group2 %in% virus.subset),
                         y.position = min(get_y_position(normalized.metrics,as.formula(paste0(variable.name,"~Variant_ordered")),ref.group = "WT") %>% pull(y.position)),
                       remove.bracket = TRUE)+
    geom_hline(yintercept = 1,color = "red",lty =2)+
    ylab(ylabel)+
    xlab("")+
    theme_classic()+
    scale_y_continuous(labels = label_number(accuracy = 0.1))+
    theme(axis.text.x = element_text(size = 12,angle = 90),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.position = "none")
  return(p)
}

w <- 5
h <- 5
# 
# plot_variable_panel(variable.name="pw90_ratio",ylabel="Contractile Duration/Calcium Duration (Normalized)")
# # 
# get_y_position(normalized.metrics,as.formula(paste0("pw90_ratio","~Variant_ordered")),ref.group = "WT")

all_variants <- unique(normalized.metrics$Variant)

ggsave("Figures/Figure3a.pdf",
       plot_variable_panel(variable.name="pw90_ratio",ylabel="Contractile Duration/Calcium Duration (Normalized)",all_variants),
       width = w,height = h)



#####
#data for figure 3a
#####

figure3a.data <- normalized.metrics %>% filter(keep_trace == 1,Variant_ordered %in% all_variants) %>%
  select(plate.ID,well,Variant_ordered,pw90_ratio) %>%
  mutate(potential.disease.category = ifelse(Variant_ordered %in% nothing.subset,"WT/Benign",
                                             ifelse(Variant_ordered %in% HCM.subset,"HCM",
                                                    ifelse(Variant_ordered %in% DCM.subset,"DCM","NA")))) %>%
  rename(Variant = Variant_ordered)

write.csv(figure3a.data,"Files_For_Source_Data/Figure3a_data.csv",row.names = FALSE)

#####
# New Figures
#####

w<-2
h<-5

wtrace <- 7
htrace <- 7

normalized.metrics <- normalized.metrics %>%
  mutate(Variant_ordered.stable = Variant_ordered,
         Variant_reversed = factor(Variant_ordered,levels = rev(levels(Variant_ordered))))

combined.traces <- combined.traces %>%
  mutate(Variant_ordered.stable = Variant_ordered,
         Variant_reversed = factor(Variant_ordered,levels = rev(levels(Variant_ordered))))

#########
#Plots for A143C/Q48P
#########

# Plot with dual y-axes
combined.q48.trace <- ggplot(combined.traces %>% filter(plate.ID == "20240925A",
                                                        subposition == "r06c04",
                                                        Virus %in% c("WT","A143C"),
                                                        well.letter == "E",
                                                        time<10) %>%
                               mutate(signal.calcium.adjusted = signal.calcium-80) %>%
                               pivot_longer(cols = c(signal.contractility, signal.calcium.adjusted), 
                                            names_to = "variable", 
                                            values_to = "signal"), 
                             aes(x = time, y = signal, color = Variant_reversed, linetype = variable)) + 
  geom_line(size = 0.7) + 
  facet_grid(cols = vars(Variant_ordered)) +
  theme_pubr() +
  theme(legend.position = "none")+
  xlab("Time (s)") + 
  ylab("Tension (Pa)") +
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_y_continuous(
    sec.axis = sec_axis(~ (.+80),name = "Calcium Fluorescence (AU)")) +
  theme()

combined.q48.trace

ggsave("Figures/Figure3j.pdf",width = wtrace,height = htrace,
       combined.q48.trace)

ggsave("Figures/Figure3k.pdf",width = w,height = h,
       plot_variable_panel("mean_amplitude.contractility","Systolic Amplitude (Normalized)",c("WT","Q48P")))



#RR145-146dup


dup.example.plot <- ggplot(combined.traces %>% filter(plate.ID == "20240925A",
                                                      subposition == "r05c05",
                                                      Variant %in% c("WT","RR145-146dup"),
                                                      well.letter == "C",
                                                      time<10),
                           aes(x = time,y = signal.contractility,color = Variant_ordered.stable)) + 
  geom_line(size = 0.7) + 
  facet_grid(cols = vars(Variant_reversed))+
  theme_pubr()+
  theme(legend.position = "none")+
  xlab("Time (s)") + ylab("Tension (Pa)")
dup.example.plot
ggsave("Figures/Figure3d.pdf",width = wtrace,height = htrace,
       dup.example.plot)

normalized.metrics <- normalized.metrics %>%
  mutate(Variant_ordered = Variant_reversed)

ggsave("Figures/Figure3e.pdf",width = w,height = h,
       plot_variable_panel("pw90_ratio","Contractile Duration/Calcium Duration (Normalized)",c("WT","RR145-146dup")))
ggsave("Figures/Figure3f.pdf",width = w,height = h,
       plot_variable_panel("mean_peak_value.contractility","Maximum Tension (Normalized)",c("WT","RR145-146dup")))
ggsave("Figures/Figure3g.pdf",width = w,height = h,
       plot_variable_panel("mean_valley.contractility","Diastolic Tension (Normalized)",c("WT","RR145-146dup")))
ggsave("Figures/Figure3h.pdf",width = w,height = h,
       plot_variable_panel("mean_amplitude.contractility","Contractile Amplitude (Normalized)",c("WT","RR145-146dup")))

normalized.metrics <- normalized.metrics %>%
  mutate(Variant_ordered = Variant_ordered.stable)






#############
# Path vs VUS
#############
normalized.metrics.path <- readRDS("Data/PathogenicScreen/CompiledData/per.well.metrics.WT.normalized.rds")
path.virus.variant.lookup <- read.csv("Data/PathogenicScreen/Lookup/TNNI3screen_pathvariants_lookuptable.csv",
                                      stringsAsFactors = FALSE)


normalized.metrics.path.labeled <- normalized.metrics.path %>%
  left_join(path.virus.variant.lookup) %>%
  filter(Gene == "TNNI3",Clinical.significance!="Phosphomimetic") %>%
  mutate(Binned.Clinical.significance = ifelse(Clinical.significance %in% c("Pathogenic","Pathogenic/Likely pathogenic","Other"),
                                               "Pathogenic/Likely pathogenic","Wildtype/Benign/Synonymous"))


VUS.summarized.metrics <- normalized.metrics %>%
  filter(Variant != "WT") %>%
  group_by(Variant) %>%
  summarise(pw90_ratio = median(pw90_ratio,na.rm = TRUE)) %>% ungroup() %>%
  ungroup() %>%
  mutate(Binned.Clinical.significance = "VUS",Clinical.significance = "VUS")

path.summarized.metrics <- normalized.metrics.path.labeled %>%
  group_by(Variant,Binned.Clinical.significance,Clinical.significance) %>%
  summarise(pw90_ratio = median(pw90_ratio,na.rm = TRUE)) %>% ungroup() %>%
  ungroup() %>%
  mutate(Binned.Clinical.significance.split = ifelse(Binned.Clinical.significance == "Pathogenic/Likely pathogenic" & pw90_ratio>1,
                                                     "Pathogenic/Likely pathogenic HCM/RCM",
                                                     ifelse(Binned.Clinical.significance == "Pathogenic/Likely pathogenic" & pw90_ratio<1,
                                                            "Pathogenic/Likely pathogenic DCM",Binned.Clinical.significance)))


all.summarized.metrics <- bind_rows(path.summarized.metrics,VUS.summarized.metrics) %>%
  mutate(Binned.Clinical.significance.split = ifelse(Binned.Clinical.significance == "VUS","VUS",Binned.Clinical.significance.split)) %>%
  mutate(Binned.Clinical.significance.split = factor(Binned.Clinical.significance.split,
                                                        levels = c("Wildtype/Benign/Synonymous",
                                                                   "Pathogenic/Likely pathogenic DCM",
                                                                   "Pathogenic/Likely pathogenic HCM/RCM",
                                                                   "VUS")))

lowest.hcm.pw90_ratio <- all.summarized.metrics %>% filter(Variant == "G203S") %>% pull(pw90_ratio)
highest.dcm.pw90_ratio <- all.summarized.metrics %>% filter(Variant == "N185K") %>% pull(pw90_ratio)

VUS.path.plot <- ggplot(all.summarized.metrics,
       aes(x = Binned.Clinical.significance.split,
           y = pw90_ratio)) + 
  geom_hline(yintercept = lowest.hcm.pw90_ratio,color = "grey",lty = 2)+#G203S pw90_ratio (lowest HCM path)
  geom_hline(yintercept = highest.dcm.pw90_ratio,color = "grey",lty = 2)+#N185K pw90_ratio (highest DCM path)
  geom_violin()+
  geom_jitter(aes(color = Binned.Clinical.significance),width = 0.1,size = 1)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
        legend.position = "none")+
  xlab("") + ylab("Contractility Duration/Calcium Duration (Normalized)")

VUS.path.plot

w <- 3
h <- 7

ggsave("Figures/Figure3b.pdf",
       VUS.path.plot,
       width = w,height = h)


figure.3eh.data <- normalized.metrics %>% filter(keep_trace == 1,Variant_ordered %in% c("WT","RR145-146dup")) %>%
  select(plate.ID,well,Variant_ordered,pw90_ratio,mean_amplitude.contractility,mean_valley.contractility,mean_peak_value.contractility)

write.csv(figure.3eh.data,"Files_For_Source_Data/Figure3efgh_data.csv",row.names = FALSE)

figure.3k.data <- normalized.metrics %>% filter(keep_trace == 1,Variant_ordered %in% c("WT","Q48P")) %>%
  select(plate.ID,well,Variant_ordered,pw90_ratio,mean_amplitude.contractility)

write.csv(figure.3k.data,"Files_For_Source_Data/Figure3k_data.csv",row.names = FALSE)
