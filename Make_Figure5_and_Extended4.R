source("Utility/contractility_kinetical_loader.R")
library(ggpubr)
library(rstatix)
library(scales)
library(ggrepel)
library(ggpmisc)
library(DescTools) # for spearman
library(glmnet)


#load clinical data
clinical.data <- read.csv("Data/PatientData/proband_data.csv",stringsAsFactors = FALSE)

###
#Statistics
###

#linear models
lm.proportion.valley.contractility <- lm(Age.of.onset ~ proportion.valley.contractility,clinical.data)

lm.fall.rate <- lm(Age.of.onset ~ fall_rate.contractility,clinical.data)

lm.mean_rise_time.calcium <- lm(Age.of.onset ~ mean_rise_time.calcium,clinical.data)

lm.pw90_ratio <- lm(Age.of.onset ~ pw90_ratio,clinical.data)

lm.fall_ratio <- lm(Age.of.onset ~ fall_ratio,clinical.data)

#extract significance
significance.cluster.variables <- 
  data.frame(variable = c("proportion.valley.contractility",
                          "fall_rate.contractility",
                          "pw90_ratio",
                          "fall_ratio",
                          "mean_rise_time.calcium"),
             p = c(summary(lm.proportion.valley.contractility)$coefficients[2,4],
                   summary(lm.fall.rate)$coefficients[2,4],
                   summary(lm.pw90_ratio)$coefficients[2,4],
                   summary(lm.fall_ratio)$coefficients[2,4],
                   summary(lm.mean_rise_time.calcium)$coefficients[2,4]),
             r2 = c(summary(lm.proportion.valley.contractility)$r.squared,
                    summary(lm.fall.rate)$r.squared,
                    summary(lm.pw90_ratio)$r.squared,
                    summary(lm.fall_ratio)$r.squared,
                    summary(lm.mean_rise_time.calcium)$r.squared)) %>%
  adjust_pvalue("holm",p.col = "p",output.col = "p.adj") %>%
  select(variable,r2,p,p.adj)

write.csv(significance.cluster.variables %>% filter(variable %in% c("proportion.valley.contractility",
                                                                    "fall_rate.contractility",
                                                                    "pw90_ratio")),
          "Files_For_Source_Data/Figure5_pvalues.csv",row.names = FALSE)
write.csv(significance.cluster.variables %>% filter(variable %in% c("fall_ratio",
                                                                    "mean_rise_time.calcium")),
          "Files_For_Source_Data/Extended_Figure4_pvalues.csv",row.names = FALSE)

####

plot.age.of.onset.two.clusters <- 
  ggplot(clinical.data %>% filter(cluster.label != "HCM_like_2"),
         aes(x = cluster.label,y = Age.of.onset)) + 
  geom_violin(aes(fill = cluster.label),alpha = 0.1,scale = "width",width = 0.5)+
  geom_point(aes(color = cluster.label))+
  scale_color_manual(values = c("#f564e3","#619cff"))+
  scale_fill_manual(values = c("#f564e3","#619cff"))+
  stat_compare_means(method = "wilcox.test",comparisons = list(c(1,2)))+
  theme_pubr(base_size = 18)+
  theme(legend.position = "none")+
  xlab("")+ylab("Age of Onset")

plot.age.of.onset.two.clusters

plot.proportion.valley.contractility.age <-
  ggplot(clinical.data,
         aes(x = proportion.valley.contractility,y = Age.of.onset)) + 
  stat_poly_line()+
  geom_point(aes(color = cluster.label)) + 
  scale_color_manual(values = c("#f564e3","#00bfc4","#619cff"))+
  
  #stat_poly_eq(label.x = "right",label.y = "top",use_label("rr","p"))+
  annotate("text",
           x = max(clinical.data$proportion.valley.contractility),
           y = max(clinical.data$Age.of.onset),
           hjust = 1,
           vjust = 1,
           label = paste0("R2 = ",signif(significance.cluster.variables %>% filter(variable == "proportion.valley.contractility") %>% pull("r2"),2),
                          "\np_adj = ",signif(significance.cluster.variables %>% filter(variable == "proportion.valley.contractility") %>% pull("p.adj"),2)))+
  theme_pubr(base_size = 18)+
  theme(legend.position = "none")+
  xlab("Diastolic Tension/Peak Tension")+ylab("Age of Onset")


plot.fall_rate.contractility.age <-
  ggplot(clinical.data,
         aes(x = fall_rate.contractility,y = Age.of.onset)) + 
  stat_poly_line()+
  geom_point(aes(color = cluster.label)) + 
  scale_color_manual(values = c("#f564e3","#00bfc4","#619cff"))+
  annotate("text",
           x = min(clinical.data$fall_rate.contractility),
           y = max(clinical.data$Age.of.onset),
           hjust = 0,
           vjust = 1,
           label = paste0("R2 = ",signif(significance.cluster.variables %>% filter(variable == "fall_rate.contractility") %>% pull("r2"),2),
                          "\np_adj = ",signif(significance.cluster.variables %>% filter(variable == "fall_rate.contractility") %>% pull("p.adj"),2)))+
  theme_pubr(base_size = 18)+
  theme(legend.position = "none")+
  xlab("Relaxation Rate")+ylab("Age of Onset")


plot.pw90_ratio.age <-
  ggplot(clinical.data,
         aes(x = pw90_ratio,y = Age.of.onset)) + 
  stat_poly_line()+
  geom_point(aes(color = cluster.label)) + 
  scale_color_manual(values = c("#f564e3","#00bfc4","#619cff"))+
  annotate("text",
           x = max(clinical.data$pw90_ratio),
           y = max(clinical.data$Age.of.onset),
           hjust = 1,
           vjust = 1,
           label = paste0("R2 = ",signif(significance.cluster.variables %>% filter(variable == "pw90_ratio") %>% pull("r2"),2),
                          "\np_adj = ",signif(significance.cluster.variables %>% filter(variable == "pw90_ratio") %>% pull("p.adj"),2)))+
  theme_pubr(base_size = 18)+
  theme(legend.position = "none")+
  xlab("Contraction Duration/Calcium Duration")+ylab("Age of Onset")


plot.fall_ratio.age <-
  ggplot(clinical.data,
         aes(x = fall_ratio,y = Age.of.onset)) + 
  stat_poly_line()+
  geom_point(aes(color = cluster.label)) + 
  scale_color_manual(values = c("#f564e3","#00bfc4","#619cff"))+
  
  #stat_poly_eq(label.x = "right",label.y = "top",use_label("rr","p"))+
  annotate("text",
           x = max(clinical.data$fall_ratio),
           y = max(clinical.data$Age.of.onset),
           hjust = 1,
           vjust = 1,
           label = paste0("R2 = ",signif(significance.cluster.variables %>% filter(variable == "fall_ratio") %>% pull("r2"),2),
                          "\np_adj = ",signif(significance.cluster.variables %>% filter(variable == "fall_ratio") %>% pull("p.adj"),2)))+
  theme_pubr(base_size = 18)+
  theme(legend.position = "none")+
  xlab("Relaxation Time/Calcium Fall Time")+ylab("Age of Onset")


plot.mean_rise_time.calcium.age <-
  ggplot(clinical.data,
         aes(x = mean_rise_time.calcium,y = Age.of.onset)) + 
  stat_poly_line()+
  geom_point(aes(color = cluster.label)) + 
  scale_color_manual(values = c("#f564e3","#00bfc4","#619cff"))+
  
  #stat_poly_eq(label.x = "right",label.y = "top",use_label("rr","p"))+
  annotate("text",
           x = max(clinical.data$mean_rise_time.calcium),
           y = max(clinical.data$Age.of.onset),
           hjust = 1,
           vjust = 1,
           label = paste0("R2 = ",signif(significance.cluster.variables %>% filter(variable == "mean_rise_time.calcium") %>% pull("r2"),2),
                          "\np_adj = ",signif(significance.cluster.variables %>% filter(variable == "mean_rise_time.calcium") %>% pull("p.adj"),2)))+
  theme_pubr(base_size = 18)+
  theme(legend.position = "none")+
  xlab("Calcium Upstroke Time")+ylab("Age of Onset")


blank_plot <- ggplot()

figure.5 <- ggarrange(blank_plot,plot.age.of.onset.two.clusters,blank_plot,
                      plot.proportion.valley.contractility.age,plot.fall_rate.contractility.age,plot.pw90_ratio.age)

ggsave(plot = figure.5,filename = "Figures/Figure5.pdf",width = 15,height = 10)

figure.extended.4 <- ggarrange(plot.fall_ratio.age,plot.mean_rise_time.calcium.age)

ggsave(plot = figure.extended.4,filename = "Figures/Extended_Data_Figure4.pdf",width = 10,height = 5)



###output data
write.csv(clinical.data %>% select(-fall_ratio,-mean_rise_time.calcium),"Files_For_Source_Data/Figure5_data.csv",row.names = FALSE)
write.csv(clinical.data %>% select(-pw90_ratio,-proportion.valley.contractility,-fall_rate.contractility),"Files_For_Source_Data/Extended_Figure4_data.csv",row.names = FALSE)
