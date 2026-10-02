source("Utility/contractility_kinetical_loader.R")
library(ggpubr)
library(rstatix)
library(scales)
library(plotROC)
library(ggbeeswarm)
#for mixed model
library(lme4)
library(lmerTest)

#load data
well.metrics <- readRDS("Data/PathogenicScreen/CompiledData/per.well.metrics.rds")
normalized.metrics <- readRDS("Data/PathogenicScreen/CompiledData/per.well.metrics.WT.normalized.rds")
combined.traces <- readRDS("Data/PathogenicScreen/CompiledData/combined.traces.rds")

#load lookup table
virus.ID.lookup <- read.csv("Data/PathogenicScreen/Lookup/TNNI3screen_pathvariants_lookuptable.csv")

#join the lookup table to the normalized data
normalized.metrics <- normalized.metrics %>%
  left_join(virus.ID.lookup %>% select(-X)) %>%
  mutate(Clinical.significance = ifelse(Clinical.significance == "Other","Pathogenic",Clinical.significance)) # rename one variant

#select out only TNNI3 variants (and not phophomimetics)
normalized.metrics.TNNI3 <- normalized.metrics %>% filter(Gene == "TNNI3",!(Clinical.significance %in% c(NA,"Phosphomimetic")))

#create a dataframe to order the median pw90 values
ordering.pw90_ratio.df <- normalized.metrics.TNNI3 %>%
  select(Variant,Clinical.significance,pw90_ratio) %>% 
  group_by(Variant,Clinical.significance) %>%
  summarise(pw90_ratio = median(pw90_ratio,na.rm=TRUE)) %>%
  mutate(Clinical.significance.ordered = factor(Clinical.significance,
                                                levels = c("Wildtype","Synonymous","Benign/Likely benign",
                                                           "Pathogenic/Likely pathogenic","Pathogenic"))) %>%
  arrange(desc(pw90_ratio))

#########
#Calculate Statistics with mixed model incorporating batch and plate
#########
variants.to.analyze <- normalized.metrics.TNNI3 %>%
  filter(Variant != "WT") %>%
  distinct(Variant) %>% pull()

p_values <- c()
group1 <- c()
group2 <- c()

for (i in 1:length(variants.to.analyze)) {
  v <- variants.to.analyze[i]
  df <- normalized.metrics.TNNI3 %>%
    filter(Variant %in% c(v,"WT"))
  
  stat.v <- lmerTest::lmer(pw90_ratio ~ Variant + (1 | plate.ID) + (1 | plate.ID:plate.letter),
                                  data = df)
  
  p_values[i] <- summary(stat.v)$coefficients[2,"Pr(>|t|)"]
  group1[i] <- "WT"
  group2[i] <- v
}

p.df <- data.frame(group1 = group1,group2 = group2,p = p_values) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance() %>%
    mutate(p.adj.signif2 = ifelse(p.adj.signif == "**","†",
                                  ifelse(p.adj.signif == "***","#",
                                         ifelse(p.adj.signif == "****","§",p.adj.signif))))


write.csv(p.df,"Files_For_Source_Data/Figure2c_pvalues.csv")

normalized.metrics.TNNI3.DCM.extra <- normalized.metrics.TNNI3 %>%
  mutate(Variant.ordered = factor(Variant,levels = ordering.pw90_ratio.df$Variant)) %>%
  mutate(Category = factor(ifelse(Clinical.significance == "Pathogenic",
                                  "Pathogenic/Likely pathogenic",
                                  Clinical.significance),
                           levels = c("Pathogenic/Likely pathogenic","Wildtype","Synonymous","Benign/Likely benign")))

#Figure 2c
pw90_screen_figure <- ggplot(
  normalized.metrics.TNNI3.DCM.extra,
       aes(x = Variant.ordered,y = pw90_ratio,color = Category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  geom_hline(yintercept = 1,color = "red",lty = 2)+
  stat_pvalue_manual(p.df,label = "p.adj.signif2",remove.bracket = TRUE,y.position = 1.6,
                     hjust = 0.5)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  scale_color_discrete(palette = c("#aa0000","#00ff00","#008a41","#00dd55"))+
  xlab("Variant")+ylab("Contractility Duration/Calcium Duration\n(Normalized)")


pw90_screen_figure

#Shape the panel so that there is room for panel 2d
plot.80.percent <- ggarrange(pw90_screen_figure,NA,widths = c(4,1))
plot.80.percent

ggsave("Figures/Figure2c.pdf",plot.80.percent,width = 13,height = 6.5)

#save data
pw90_screen_figure_data <- normalized.metrics.TNNI3.DCM.extra %>%
  select(plate.ID,plate.letter,well,Variant.ordered,pw90_ratio)

write.csv(pw90_screen_figure_data,"Files_For_Source_Data/Figure2c_data.csv",row.names = FALSE)

#############
#Panels 2e-h
#############

#Create per-variant data frame
violin.df <- normalized.metrics.TNNI3.DCM.extra %>%
  group_by(Variant) %>%
  summarise(across(where(is.numeric),function(x) {median(x,na.rm = TRUE)})) %>%
  mutate(Category = factor(ifelse(pw90_ratio > 1.02,"HCM-like",
                           ifelse(pw90_ratio<0.92,"DCM-like","WT/Benign")),
                           levels = c("WT/Benign","DCM-like","HCM-like")))

max.wt.pw90_ratio <- max(violin.df %>% filter(Category == "WT/Benign") %>% pull(pw90_ratio))
min.wt.pw90_ratio <- min(violin.df %>% filter(Category == "WT/Benign") %>% pull(pw90_ratio))

max.wt.amplitude <- max(violin.df %>% filter(Category == "WT/Benign") %>% pull(mean_amplitude.contractility))
min.wt.amplitude <- min(violin.df %>% filter(Category == "WT/Benign") %>% pull(mean_amplitude.contractility))

max.wt.valley <- max(violin.df %>% filter(Category == "WT/Benign") %>% pull(mean_valley.contractility))
min.wt.valley <- min(violin.df %>% filter(Category == "WT/Benign") %>% pull(mean_valley.contractility))

max.wt.peak <- max(violin.df %>% filter(Category == "WT/Benign") %>% pull(mean_peak_value.contractility))
min.wt.peak <- min(violin.df %>% filter(Category == "WT/Benign") %>% pull(mean_peak_value.contractility))

violin.pw90_ratio <-
ggplot(violin.df,
       aes(x = Category,y = pw90_ratio,color = Category)) + 
  geom_hline(yintercept = max.wt.pw90_ratio,col = "red",lty = 2)+
  geom_hline(yintercept = min.wt.pw90_ratio,col = "red",lty = 2)+
  geom_violin() + 
  geom_beeswarm()+
  theme_pubr()+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  xlab("")+ylab("Contractility Duration/\nCalcium Duration\n(Normalized)")

violin.amplitude <-
ggplot(violin.df,
       aes(x = Category,y = mean_amplitude.contractility,color = Category)) + 
  geom_hline(yintercept = max.wt.amplitude,col = "red",lty = 2)+
  geom_hline(yintercept = min.wt.amplitude,col = "red",lty = 2)+
  geom_violin() + 
  geom_beeswarm()+
  theme_pubr()+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  xlab("")+ylab("Contractile Amplitude\n(Normalized)")

violin.peak <-
ggplot(violin.df,
       aes(x = Category,y = mean_peak_value.contractility,color = Category)) + 
  geom_hline(yintercept = max.wt.peak,col = "red",lty = 2)+
  geom_hline(yintercept = min.wt.peak,col = "red",lty = 2)+
  geom_violin() + 
  geom_beeswarm()+
  theme_pubr()+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  xlab("")+ylab("Maximum Tension\n(Normalized)")

violin.valley <-
ggplot(violin.df,
       aes(x = Category,y = mean_valley.contractility,color = Category)) + 
  geom_hline(yintercept = max.wt.valley,col = "red",lty = 2)+
  geom_hline(yintercept = min.wt.valley,col = "red",lty = 2)+
  geom_violin() + 
  geom_beeswarm()+
  theme_pubr()+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  xlab("")+ylab("Diastolic Tension\n(Normalized)")

violin.plots <- ggarrange(violin.pw90_ratio,violin.peak,violin.amplitude,violin.valley,ncol = 4)

ggsave(plot = violin.plots,filename = "Figures/Figure2efgh.pdf",width = 13,height = 4)

violin.df.forfigures <- violin.df %>%
  select(Category,Variant,
         pw90_ratio,
         mean_peak_value.contractility,
         mean_amplitude.contractility,
         mean_valley.contractility)

write.csv(violin.df.forfigures,"Files_For_Source_Data/Figure2efgh_data.csv",row.names = FALSE)

