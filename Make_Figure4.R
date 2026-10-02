source("Utility/contractility_kinetical_loader.R")
library(ggpubr)
library(rstatix)

library(randomForest)
library(mclust)
library(patchwork)
library(cluster)

library(ggrepel)
library(ggbeeswarm)
library(ggpmisc)
library(scales)

#load data
well.metrics <- readRDS("Data/PathogenicScreen/CompiledData/per.well.metrics.rds")
normalized.metrics <- readRDS("Data/PathogenicScreen/CompiledData/per.well.metrics.WT.normalized.rds")

#load lookup table
virus.ID.lookup <- read.csv("Data/PathogenicScreen/Lookup/TNNI3screen_pathvariants_lookuptable.csv")

#join the lookup table to the normalized data
normalized.metrics.TNNI3 <- normalized.metrics %>%
  left_join(virus.ID.lookup %>% select(-X)) %>%
  mutate(Clinical.significance = ifelse(Clinical.significance == "Other","Pathogenic",Clinical.significance)) %>%# rename one variant
  filter(Gene == "TNNI3",!(Clinical.significance %in% c(NA,"Phosphomimetic")))

normalized.metrics.per.variant <- normalized.metrics.TNNI3 %>%
  group_by(Variant,Clinical.significance) %>%
  summarise(across(where(is.numeric),~median(.x,na.rm = TRUE))) %>%
  mutate(Clinical.significance = ifelse(Clinical.significance == "Other","Pathogenic",Clinical.significance)) %>% # rename one variant
  mutate(wt.or.benign = Clinical.significance %in% c("Synonymous","Benign/Likely benign","Wildtype")) %>%
  ungroup()

working.df <- normalized.metrics.per.variant %>%
  select(-contains("sd_"),#get rid of sds
         -mean_E,-mean_G,-mean_nu,-mean_K,#get rid of broken variables
         -contains("pw50"),#get rid of arbitrary variables
         -contains("pw30"), #get rid of arbitrary variables
         -contains("baseline"), #get rid of extra variables (encoded in valley)
         -keep_trace,-keep_trace_ppqtran,
         -Clinical.significance,-AA.pos
  )

working.df.clean <- na.omit(working.df)

#############

# Run random forest a number of times to ensure pattern is stable
set.seed(42)

variable.ranks <- rep(0,length(names(working.df.clean))-2) # -2 for Variant and wt.or.benign
variable.importance <- rep(0,length(names(working.df.clean))-2) # -2 for Variant and wt.or.benign
names(variable.ranks) <- names(working.df.clean %>% select(-Variant,-wt.or.benign))
names(variable.importance) <- names(working.df.clean %>% select(-Variant,-wt.or.benign))

nrf <- 100

for (i in 1:nrf) {
  rf.wt <- randomForest(wt.or.benign ~ .,data = working.df.clean %>% select(-Variant))
  
  #plot to look at importance
  sorted.importance <- sort(rf.wt$importance[,1],decreasing = TRUE)
  
  #normalize to sum of all importance so that you can compare between random forests
  total.importance <- sum(sorted.importance)
  sorted.importance <- sorted.importance / total.importance 
  
  sorted.names <- names(sorted.importance)
  
  for (j in 1:length(sorted.names)) {
    v <- sorted.names[j]
    variable.ranks[v] <- variable.ranks[v] + j
    variable.importance[v] <- variable.importance + sorted.importance[v]
  }
  #normalized importance
  
}

sorted.variable.ranks <- sort(variable.ranks)
sorted.variable.importance <- sort(variable.importance,decreasing = TRUE)

importance.df <- data.frame(Average_Importance = sorted.variable.importance/nrf,
                            variable = names(sorted.variable.importance),
                            importance.rank = 1:length(sorted.variable.importance))

named.variables.for.plot <- 
  c(proportion.valley.contractility = "Diastolic Tension/Peak Tension",
    mean_rise_time.calcium = "Upstroke Time Calcium",
    fall_rate.contractility = "Relaxation Rate",
    pw90_ratio = "Contraction Duration/Calcium Duration",
    fall_ratio = "Relaxation Time/Calcium Fall Time")


importance.legend <- paste(named.variables.for.plot[names(sorted.variable.importance)[1]],
                           named.variables.for.plot[names(sorted.variable.importance)[2]],
                           named.variables.for.plot[names(sorted.variable.importance)[3]],
                           named.variables.for.plot[names(sorted.variable.importance)[4]],
                           named.variables.for.plot[names(sorted.variable.importance)[5]],
                           sep = "\n")

point.conversion = 2.845

importance.plot <-
  ggplot(importance.df,
         aes(x = importance.rank,y = Average_Importance,
             label = ifelse(importance.rank<=5,named.variables.for.plot[variable],NA))) + 
  geom_point(point.size = 3)+
  geom_label_repel(size = 18/point.conversion,
                   min.segment.length = 0,
                   point.size = 3)+
  theme_pubr(base_size = 14)+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+
  xlab("Ranked Importance")+ylab("Average Importance")

importance.plot

important_variables <- names(sorted.variable.importance)[1:5]

#############

#####
# shuffle through variables
#####
working.df.with.groups <- working.df.clean %>%
  mutate(known.ids = ifelse(wt.or.benign,1,2)) %>%
  select(-wt.or.benign)

extra.in.wt <- c()
missing.from.wt <- c()
k.vec <- c()
plot.index <- c()
sil.vec <- c()
n.unique.dcm.clusters.vec <- c()

counter <- 1
for (k in 3:10) {
  X_selected <- working.df.with.groups[,important_variables]
  
  X_scaled <- scale(X_selected)
  
  d <- dist(X_scaled)
  hc <- hclust(d, method="ward.D2")
  working.df.with.groups$cluster <- cutree(hc, k=k)
  sil <- silhouette(working.df.with.groups$cluster,d)
  
  pca <- prcomp(X_scaled)
  working.df.with.groups$pca1 <- pca$x[,1]
  working.df.with.groups$pca2 <- pca$x[,2]
  
  #check fidelity of clustering
  #WT group should only contain WT and benign variants
  wt.cluster.id <- working.df.with.groups %>% filter(Variant == "WT") %>% pull(cluster)
  wt.cluster.members <- working.df.with.groups %>% filter(cluster == wt.cluster.id)
  number.extra.in.wt <- sum(wt.cluster.members$known.ids != 1)
  missing.in.wt <- sum((working.df.with.groups %>% filter(known.ids==1) %>% pull(cluster))!=wt.cluster.id)
  
  k.vec[counter] <- k
  extra.in.wt[counter] <- number.extra.in.wt
  missing.from.wt[counter] <- missing.in.wt
  plot.index[counter] <- counter
  sil.vec[counter] <- mean(sil[, "sil_width"])
  n.unique.dcm.clusters.vec[counter] <- length(unique(working.df.with.groups %>% 
                                                        filter(Variant %in% c("K36Q","E182K","N185K")) %>%
                                                        pull(cluster)))  
  
  knowns <- working.df.with.groups[working.df.with.groups$known.ids %in% c(1),]
  ari <- adjustedRandIndex(knowns$known.ids, knowns$cluster)
  counter <- counter + 1
}


cluster.info <- data.frame(nvar = 5,
                           ncluster = k.vec,
                           extra.in.wt = extra.in.wt,
                           missing.from.wt = missing.from.wt,
                           plot.index = plot.index,
                           silhouette.width = sil.vec,
                           dcm.clusters = n.unique.dcm.clusters.vec)

good.clusters <- cluster.info %>%
  filter(extra.in.wt == 0,missing.from.wt == 0)

good.clusters %>% arrange(desc(silhouette.width))

for (i in 1:length(good.clusters$plot.index)) {
  nvar <- good.clusters$nvar[i]
  ncluster <- good.clusters$ncluster[i]
}

#with this metric, nvar = 5, nclusters = 6 has the best silhouette
###########
#check best silhouette 
working.df.with.groups <- working.df.clean %>%
  mutate(known.ids = ifelse(wt.or.benign,1,2)) %>%
  select(-wt.or.benign)

k <- 6

X_selected <- working.df.with.groups[,important_variables]

X_scaled <- scale(X_selected)

d <- dist(X_scaled)
hc <- hclust(d, method="ward.D2")
working.df.with.groups$cluster <- cutree(hc, k=k)
sil <- silhouette(working.df.with.groups$cluster,d)

pca <- prcomp(X_scaled)
working.df.with.groups$pca1 <- pca$x[,1]
working.df.with.groups$pca2 <- pca$x[,2]
working.df.with.groups$pca3 <- pca$x[,3]
working.df.with.groups$pca4 <- pca$x[,4]
working.df.with.groups$pca5 <- pca$x[,5]

working.df.labeled <- working.df.with.groups %>%
  mutate(cluster.label = factor(
    case_when(
      cluster == 3 ~ "DCM_like_2",#b79f00
      cluster == 6 ~ "DCM_like_1",#f8766d
      cluster == 2 ~ "WT/Benign",#00ba38
      cluster == 1 ~ "HCM_like_1",#f564e3
      cluster == 5 ~ "HCM_like_2",#00bfc4
      cluster == 4 ~ "HCM_like_3"#619cff
    ),levels = c("WT/Benign","DCM_like_1","DCM_like_2","HCM_like_1","HCM_like_2","HCM_like_3"))) %>%
  mutate(cluster.color = case_when(
    cluster == 3 ~ "#b79f00",
    cluster == 6 ~ "#f8766d",
    cluster == 2 ~ "#00ba38",
    cluster == 1 ~ "#f564e3",
    cluster == 5 ~ "#00bfc4",
    cluster == 4 ~ "#619cff"
  ))

pca.plot12 <- ggplot(working.df.labeled, 
                     aes(x = pca1, y = pca2, color = cluster.color, label = Variant)) +
  #geom_point(size = 3) +
  geom_point()+
  geom_label_repel(aes(fill = col_mix(cluster.color,"#FFFFFF",0.9)),
                   box.padding = unit(0.1,"lines"),
                   point.size = 3,
                   min.segment.length = 0)+
  #labs(title = paste0("nvar = ",i," k = ",k))+
  theme_pubr(base_size = 14)+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+
  scale_color_identity()+scale_fill_identity()+
  xlab("PC1")+ylab("PC2")

pca.plot12

#this rearranges the important variables
important_variables_save <- important_variables

important_variables_ordered <- c("pw90_ratio",
                                 "fall_ratio",
                                 "proportion.valley.contractility",
                                 "fall_rate.contractility",
                                 "mean_rise_time.calcium")
important_variables_ordered.names <- c("Contraction Duration/\nCalcium Duration",
                                       "Relaxation Time/\nCalcium Fall Time",
                                       "Diastolic Tension/\nPeak Tension",
                                       "Relaxation Rate",
                                       "Upstroke Time Calcium"
)


combined_plot <-ggplot()
var.plots <- list()
for (j in 1:length(important_variables_ordered)) {
  v <- important_variables_ordered[j]
  vn <- important_variables_ordered.names[j]
  var.violin <- ggplot(working.df.labeled,
                        aes(x = cluster.label,y = .data[[v]],color = cluster.color)) + 
    geom_vline(xintercept = c(1.6,3.5),lty = "dotted",col = "black")+
    geom_violin(aes(fill = col_mix(cluster.color,"#FFFFFF",0.9)))+
    geom_beeswarm()+
    theme_pubr(base_size = 14)+
    theme(#axis.text.x = element_text(angle = 45,vjust = 0.5,hjust=0.5,size = 18),,
      axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1,size = 18),
      axis.title.y = element_text(size = 18),
      legend.position = "none")+
    scale_color_identity()+scale_fill_identity()+
    xlab("")+ylab(vn)
  
  var.plots[[j]] <- var.violin
}  


#New for revision
# cutoffs for groups in the proportion.valley.contractility/fall_ratio plot
diastole.cutoff <- 
  (max(working.df.labeled %>% filter(cluster.label == "HCM_like_1") %>% pull(proportion.valley.contractility))+
     min(working.df.labeled %>% filter(cluster.label == "HCM_like_2") %>% pull(proportion.valley.contractility)))/2

hcm.cutoff <-
  (max(working.df.labeled %>% filter(cluster.label == "WT/Benign") %>% pull(fall_ratio))+
     min(working.df.labeled %>% filter(cluster.label == "HCM_like_1") %>% pull(fall_ratio)))/2

dcm.cutoff <-
  (min(working.df.labeled %>% filter(cluster.label == "WT/Benign") %>% pull(fall_ratio))+
     max(working.df.labeled %>% filter(cluster.label == "DCM_like_1") %>% pull(fall_ratio)))/2

hcm.2.3.cutoff <-
  (min(working.df.labeled %>% filter(cluster.label == "HCM_like_3") %>% pull(fall_ratio))+
     max(working.df.labeled %>% filter(cluster.label == "HCM_like_2") %>% pull(fall_ratio)))/2

diastole.concal.df <-
  working.df.labeled %>% select(Variant,proportion.valley.contractility,fall_ratio,cluster.color)


#dummy points create artificial barriers so the labels don't cross lines
#they are not shown given Variant = ""
dummy_points.hcm <- data.frame(
  Variant = "",
  fall_ratio = hcm.cutoff,
  proportion.valley.contractility = seq(min(diastole.concal.df$proportion.valley.contractility),
                                        max(diastole.concal.df$proportion.valley.contractility),
                                        length.out = 100),
  cluster.color = "#FFFFFF"
)

dummy_points.dcm <- data.frame(
  Variant = "",
  fall_ratio = dcm.cutoff,
  proportion.valley.contractility = seq(min(diastole.concal.df$proportion.valley.contractility),
                                        max(diastole.concal.df$proportion.valley.contractility),
                                        length.out = 100),
  cluster.color = "#FFFFFF"
)

dummy_points.hcm.2.3 <- data.frame(
  Variant = "",
  fall_ratio = hcm.2.3.cutoff,
  proportion.valley.contractility = seq(min(diastole.concal.df$proportion.valley.contractility),
                                        max(diastole.concal.df$proportion.valley.contractility),
                                        length.out = 100),
  cluster.color = "#FFFFFF"
)

dummy_points.diastole <- data.frame(
  Variant = "",
  proportion.valley.contractility = diastole.cutoff,
  fall_ratio = seq(min(diastole.concal.df$fall_ratio)-0.01,
                   max(diastole.concal.df$fall_ratio),
                   length.out = 100),
  cluster.color = "#FFFFFF"
)

diastole.concal.df.withdummy <- diastole.concal.df %>%
  bind_rows(dummy_points.hcm,dummy_points.dcm,dummy_points.diastole,dummy_points.hcm.2.3)

line.alpha <- 0.5

diastole.concal.panel <-
  ggplot(diastole.concal.df.withdummy, 
         aes(x = proportion.valley.contractility, y = fall_ratio, color = cluster.color, label = Variant)) +
  #geom_point(size = 3) +
  
  geom_point()+
  #geom_label_repel(aes(fill = alpha(cluster.color,0.2)))+
  geom_label_repel(aes(fill = col_mix(cluster.color,"#FFFFFF",0.9)),
                   box.padding = unit(0.1,"lines"),
                   point.size = 3,
                   min.segment.length = 0)+
  geom_hline(yintercept = hcm.cutoff,lty = 2,col = "red",alpha = line.alpha)+
  geom_hline(yintercept = dcm.cutoff,lty = 2,col = "red",alpha = line.alpha)+
  geom_vline(xintercept = diastole.cutoff,lty = 2,col = "red",alpha = line.alpha)+
  geom_hline(yintercept = hcm.2.3.cutoff,lty = 2,col = "red",alpha = line.alpha)+
  #labs(title = paste0("nvar = ",i," k = ",k))+
  theme_pubr(base_size = 14)+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+
  scale_color_identity()+scale_fill_identity()+
  xlab("Diastolic Tension/Peak Tension")+ylab("Relaxation Time/Calcium Fall Time")

diastole.concal.panel


arranged.plots<-
  ggarrange(ggarrange(importance.plot,pca.plot12,diastole.concal.panel,ncol = 3,widths = c(2.5,3,3)),
            ggarrange(var.plots[[1]],var.plots[[2]],var.plots[[3]],var.plots[[4]],var.plots[[5]],ncol = 5,nrow = 1),
            nrow = 2,heights = c(2,1))


scale.factor.arrange <- 1.1
ggsave(plot = arranged.plots,filename = "Figures/Figure4.pdf",width = 19*scale.factor.arrange,height = 14.25*scale.factor.arrange)

##########


#Save data

figure.4.data <- working.df.labeled %>% select(Variant,
                                               cluster.label,
                                               pw90_ratio,
                                               fall_ratio,
                                               proportion.valley.contractility,
                                               fall_rate.contractility,
                                               mean_rise_time.calcium,
                                               pca1,
                                               pca2)

write.csv(figure.4.data,"Files_For_Source_Data/Figure4_data.csv",row.names = FALSE)
