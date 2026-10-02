library(tidyverse)
library(ggprism)
library(ggbeeswarm)
library(cowplot)

ClinVar.file <- "Data/ClinVarAnnotation/clinvar_result_20251004.txt"
ClinVar.df <- read.delim(ClinVar.file)
ClinVar.df <- filter(ClinVar.df,Gene.s. == "TNNI3")

TNNI3_domains <- read_csv("Data/ClinVarAnnotation/TNNI3_domains.csv")
AA_list <- read_csv("Data/ClinVarAnnotation/AA_list.csv")

Order_Clinical.significance <- c("Pathogenic/Likely pathogenic", "Benign/Likely benign", "Variant of Uncertain Signifiance")
Order_DomainName <- c("Cardiac-specific N-terminal extension", "Unspecified", "N-terminal conserved region", "TnT-binding region", "Switch region", "C-terminal mobile domain")
Order_DomainResidue <- c("[1-30]", "[31-41]", "[42-65]", "[66-136]", "[137-163]", "[164-210]")


ClinVar.df <- ClinVar.df %>% 
  mutate(AA.wt = ClinVar.df$Protein.change %>%  substr(1,1),
         AA.change = ClinVar.df$Protein.change %>%  substr(nchar(ClinVar.df$Protein.change), nchar(ClinVar.df$Protein.change)),
         AA.pos = ClinVar.df$Protein.change %>%  substr(2, nchar(ClinVar.df$Protein.change)-1) %>% as.numeric(),
         DCM = grepl("D|dilated",ClinVar.df$Condition.s.),
         HCM = grepl("H|hypertrophic",ClinVar.df$Condition.s.),
         RCM = grepl("R|restrictive",ClinVar.df$Condition.s.),
         Clinical.significance = sub("\\(.*","",ClinVar.df$Germline.classification),
         Domain = case_when(
           AA.pos >= 1 & AA.pos <= 30 ~ "Cardiac-specific N-terminal extension [1-30]",
           AA.pos >= 42 & AA.pos <= 65 ~ "N-terminal conserved region [42-65]",
           AA.pos >= 66 & AA.pos <= 136 ~ "TnT-binding region [66-136]",
           AA.pos >= 137 & AA.pos <= 163 ~ "Switch region [137-163]",
           AA.pos >= 164 & AA.pos <= 210 ~ "C-terminal mobile domain [164-210]",
           TRUE ~ "Unspecified [31-41]"
         )) %>% 
  mutate(Clinical.significance_simple = case_when(
    Clinical.significance == "Benign/Likely benign" ~ "Benign/Likely benign",
    Clinical.significance == "Likely benign" ~ "Benign/Likely benign",
    Clinical.significance == "Pathogenic" ~ "Pathogenic/Likely pathogenic",
    Clinical.significance == "Likely pathogenic" ~ "Pathogenic/Likely pathogenic",
    Clinical.significance == "Pathogenic/Likely pathogenic" ~ "Pathogenic/Likely pathogenic",
    Clinical.significance == "Conflicting classifications of pathogenicity" ~ "Variant of Uncertain Signifiance",
    Clinical.significance == "Uncertain significance" ~ "Variant of Uncertain Signifiance",
    TRUE ~ "Other"))

ClinVar.tidy_df <- ClinVar.df %>% select(Name, Gene.s., AA.wt, AA.change, AA.pos, DCM, HCM, RCM, Clinical.significance, Clinical.significance_simple, Domain)

####

ann1 <- data.frame(AA.pos = c(1,31,42,66,137,164), AAend = c(30,41,65,136,163,210), ypos = -3, lab1 = Order_DomainName, lab2 = Order_DomainResidue)

# Histogram plot (columns)
pathogenic.plot.barplot <- 
  ClinVar.tidy_df %>% 
  ggplot(aes(x = AA.pos, fill = Clinical.significance_simple)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous("Amino acid position", breaks = seq(0, 210, by = 10), guide = "prism_minor", minor_breaks = seq(0, 210, 1)) + 
  scale_y_continuous(breaks = seq(0, 5, by = 1)) +
  scale_fill_manual("Clinical significance", breaks = Order_Clinical.significance, 
                    values = c("red", "green3", "grey75")) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        text = element_text(size = 24),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.length.x = unit(0.25, "cm"),
        strip.background = element_rect(fill = "white"),
        legend.position = "top"
  )

p_ann <- ggplot(ann1) +
  geom_rect(
    aes(xmin = AA.pos, xmax = AAend, ymin = 0, ymax = 1),
    fill = "grey80", color = "black"
  ) +
  geom_text(
    aes(x = (AA.pos + AAend) / 2, y = 0.5, label = paste0(lab1,"\n",lab2)),
    size = 4
  ) +
  xlim(min(ClinVar.tidy_df$AA.pos), max(ClinVar.tidy_df$AA.pos)) +
  theme_void() +  # no axes etc.
  theme(plot.margin = margin(0, 1, 0, 1))  # 0=top,right,bot,left

labeled.plot <- plot_grid(pathogenic.plot.barplot, p_ann, ncol = 1, align = 'v', rel_heights = c(1, 0.13))

ggsave("Figures/Figure2a.pdf", plot = labeled.plot, width = 9000, height = 1600, units = "px")

