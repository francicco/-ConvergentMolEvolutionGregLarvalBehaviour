library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
library(dplyr)
library(tidyverse)
library(awtools)
library(UpSetR)
library(ggExtra)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(smatr)

########################## BUSTED-PH #################################
BUSTED<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/BUSTED-PH.Results.tsv", header=TRUE)
head(BUSTED)
length(BUSTED$Gene)

ordBUSTED <- BUSTED %>%
  arrange(PvalDiff)

# Perform Bonferroni correction on PvalDiff
ordBUSTED <- ordBUSTED %>%
  mutate(Ratio = OmegaTest/OmegaBackground)

ordBUSTED <- ordBUSTED %>%
  mutate(PvalDiff_Bonferroni = pmin(1, p.adjust(PvalDiff, method = "bonferroni")))

ordBUSTED <- ordBUSTED %>%
  mutate(PvalDiff_BH = pmin(1, p.adjust(PvalDiff, method = "BH")))

# Substitue 0 values with the smallest other value
ordBUSTED <- ordBUSTED %>%
  mutate(PvalDiff_BH = if_else(PvalDiff_BH == 0, min(PvalDiff_BH[PvalDiff_BH != 0]), PvalDiff_BH))

ordBUSTED <- ordBUSTED %>%
  mutate(LogPvalDiff_BH = log10(PvalDiff_BH))


length(ordBUSTED$Gene)

#write.table(ordBUSTED, file = "~/Documents/HeliconiiniProg/GregariusConvergence/BUSTED-PH.Results.Corrected.tsv", sep = "\t", row.names = FALSE)
ordBUSTED<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/BUSTED-PH.Results.Corrected.tsv", header=TRUE)
# Make 2 cathegories
# Different and episodic selection in Gregarious Putative Convergence
gene_SelecInGreg <- ordBUSTED %>%
  filter(PvalBackground > 0.01 & PvalTest < 0.01 & PvalDiff_BH < 0.05) %>%
  pull(Gene)

subst_SelecInGreg <- ordBUSTED %>%
  filter(Gene %in% gene_SelecInGreg) %>%
  mutate(Type = "Convergence") # mutate(Type = "DiffSelecInGreg")
# Different but no episodic selection in Gregarious
gene_NoSelecInBoth <- ordBUSTED %>%
  filter(PvalBackground > 0.01 & PvalTest > 0.01 & PvalDiff_BH < 0.05) %>%
  pull(Gene)

subst_NoSelecInBoth <- ordBUSTED %>%
  filter(Gene %in% gene_NoSelecInBoth) %>%
  mutate(Type = "Convergence") # mutate(Type = "DiffNoSelecInBoth")
# Different but episodic selection in Non-Gregarious  #### NO Convergence
gene_SelecInNoGreg <- ordBUSTED %>%
  filter(PvalBackground < 0.01 & PvalTest > 0.01 & PvalDiff_BH < 0.05) %>%
  pull(Gene)

subst_SelecInNoGreg <- ordBUSTED %>%
  filter(Gene %in% gene_SelecInNoGreg) %>%
  mutate(Type = "NoConvergence") # mutate(Type = "DiffSelecInNoGreg")
# Different and episodic selection in both #### NO Convergence
gene_BothSelec <- ordBUSTED %>%
  filter(PvalBackground < 0.01 & PvalTest < 0.01 & PvalDiff_BH < 0.05) %>%
  pull(Gene)

subst_BothSelec <- ordBUSTED %>%
  filter(Gene %in% gene_BothSelec) %>%
  mutate(Type = "NoConvergence") # mutate(Type = "BothSelec")
# No Difference  #### NO Convergence
gene_NoDiff <- ordBUSTED %>%
  filter(PvalDiff_BH > 0.05) %>%
  pull(Gene)

subst_NoDiff <- ordBUSTED %>%
  filter(Gene %in% gene_NoDiff) %>%
  mutate(Type = "NoConvergence") # mutate(Type = "NoDiff")

# Concatenate dataframes
length(subst_SelecInGreg$Gene)
length(subst_SelecInNoGreg$Gene)
length(subst_BothSelec$Gene)
length(subst_NoSelecInBoth$Gene)
length(subst_NoDiff$Gene)

BUSTED_extended <- bind_rows(subst_SelecInGreg,subst_NoSelecInBoth,subst_SelecInNoGreg,subst_BothSelec,subst_NoDiff)
#write.table(BUSTED_extended, file = "~/Documents/HeliconiiniProg/GregariusConvergence/BUSTED-PH.Results.Corrected.Extended.tsv", sep = "\t", row.names = FALSE)
length(BUSTED_extended$Gene)
# View the dataframe with the adjusted p-values
unique(BUSTED_extended$Type)
head(BUSTED_extended)

slope.FusedShort<-sma(OmegaBackground~OmegaTest*Type, data=BUSTED_extended, log='xy', multcomp=FALSE)
summary(slope.FusedShort)
NoConvEle<-round(slope.FusedShort$coef$NoConvergence$`coef(SMA)`[1],2)
NoConvSlop<-round(slope.FusedShort$coef$NoConvergence$`coef(SMA)`[2],2)
NoConvR2<-round(slope.FusedShort$r2$NoConvergence[1],2)
NoConveq <- paste("y=", NoConvSlop, "*x", NoConvEle, " | R2: ",NoConvR2, sep='')
NoConveq
ConvEle<-round(slope.FusedShort$coef$Convergence$`coef(SMA)`[1],2)
ConvSlop<-round(slope.FusedShort$coef$Convergence$`coef(SMA)`[2],2)
ConvR2<-round(slope.FusedShort$r2$Convergence[1],2)
Conveq <- paste("y=", ConvSlop, "*x", ConvEle, " | R2: ",ConvR2, sep='')
Conveq
plot(slope.FusedShort)

# Fit linear models
subst_Convergence <- BUSTED_extended %>%
  filter(Type == "Convergence")

subst_NoConvergence <- BUSTED_extended %>%
  filter(Type == "NoConvergence")

lm_modelsConv <- lm(OmegaBackground ~ OmegaTest, data = subst_NoConvergence)
lm_modelsConv
coef(lm_modelsConv$coefficients[0])
# Extract coefficients and format equation text
eq_texts <- paste("Equation:", round(coef(lm_models)[2], 2), "* x +", round(coef(lm_models)[1], 2))

?wilcox.test
wilcox.test(subst_NoConvergence$OmegaBackground, subst_Convergence$OmegaBackground, alternative='greater')
median(subst_NoConvergence$OmegaBackground)
median(subst_Convergence$OmegaBackground)
wilcox.test(subst_NoConvergence$OmegaTest, subst_Convergence$OmegaTest, alternative='less')
median(subst_NoConvergence$OmegaTest)
median(subst_Convergence$OmegaTest)

p_values <- c(0.02376, 5.688e-10)
p.adjust(p_values, method = "bonferroni")

#desired_order <- c("DiffSelecInGreg","DiffNoSelecInBoth","DiffSelecInNoGreg","BothSelec","NoDiff")  # Replace with your desired order
#desired_order <- c("Convergent","DiffSelecInNoGreg","BothSelec","NoDiff")  # Replace with your desired order
desired_order <- c("NoConvergence","Convergence")  # Replace with your desired order

# Convert the 'Type' variable to a factor with the desired order
BUSTED_extended$Type <- factor(BUSTED_extended$Type, levels = desired_order)
#my_colors <- c("#ff1f1f", "#7b1727", "#514d72", "#170057","black")
#my_colors <- c("#ff1f1f", "#514d72", "#170057","black")
my_colors <- c("#a5c4cf","#c82b27")

pmain <- ggplot(BUSTED_extended, aes(x = OmegaTest, y = OmegaBackground, color = Type))+
  geom_point(alpha=0.4)+
  scale_x_log10() +
  scale_y_log10() +
  stat_smooth(data = BUSTED_extended, aes(x = OmegaTest, color = Type),
              method = "lm", se = TRUE) +
  geom_text(data = subset(BUSTED_extended, OmegaBackground < 0.01 & Type == "Convergence"), 
            aes(label = Gene), 
            hjust = -0.1, vjust = 0.1) + # Add text labels next to the linear model
  #geom_text(x = 0.01, y = 0.5, label = NoConveq, color = c("#a5c4cf")) +
  #geom_text(x = 1, y = 0.5, label = Conveq, color = c("#c82b27")) +
  scale_color_manual(values = my_colors) + a_plex_theme(base_size=12,axis_title_size=12)
pmain

subst_ConOBsum
subst_ConOBsum["Median"]
subst_NoConOBsum["Median"]
xdens <- axis_canvas(pmain, axis = "x") +
  geom_vline(xintercept = median(subst_Convergence$OmegaTest), color = "#c82b27") +
  geom_vline(xintercept = median(subst_NoConvergence$OmegaTest), color = "#a5c4cf") +
  geom_density(data = BUSTED_extended, aes(x = OmegaTest, fill = Type),
               alpha = 0.3, size = 0.2)+
  scale_x_log10() +
  scale_fill_manual(values = my_colors) + a_plex_theme(base_size=12,axis_title_size=12)
xdens
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
  geom_vline(xintercept = median(subst_Convergence$OmegaBackground), color = "#c82b27") +
  geom_vline(xintercept = median(subst_NoConvergence$OmegaBackground), color = "#a5c4cf") +
  geom_density(data = BUSTED_extended, aes(x = OmegaBackground, fill = Type),
               alpha = 0.3, size = 0.2)+
  coord_flip()+
  scale_x_log10() +
  scale_fill_manual(values = my_colors) + a_plex_theme(base_size=12,axis_title_size=12)
ydens
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

########################## BUSTED-PH + CSUBST #################################
subst<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/All.csubst_cb_2_OCNa2s1.0omegaCa2s1.0_filtered.tsv", header=TRUE)
head(subst)

length((subst %>% distinct(OG))$OG)

ggplot(subst, aes(x = omegaCspe2any)) +
  geom_density() +
  xlim(0,5) +
  theme_minimal()

#data = subset(ordBUSTED, PvalBackground > 0.01 & PvalTest < 0.01 & PvalDiff_BH < 0.05

head(ordBUSTED)

# Filter BUSTED-PH for Pvalue_BH < 0.05
DiffSelPress <- BUSTED_extended %>%
  filter(Type == "Convergence") %>%
  pull(Gene)

length(DiffSelPress)

NoDiffSelPress <- BUSTED_extended %>%
  filter(Type == "NoConvergence") %>%
  pull(Gene)

length(NoDiffSelPress)

# Select data in df_A matching the filtered genes
subst_DiffSelPress <- subst %>%
  filter(OG %in% DiffSelPress) %>%
  mutate(Type = "DiffSelPress")

head(subst_DiffSelPress)

subst_NoDiffSelPress <- subst %>%
  filter(OG %in% NoDiffSelPress) %>%
  mutate(Type = "NoDiffSelPress")

head(subst_NoDiffSelPress)

# Concatenate dataframes
BUSTED_csubst <- bind_rows(subst_DiffSelPress, subst_NoDiffSelPress)

(BUSTED_csubst %>% distinct(Type))$Type
head(BUSTED_csubst)

# Find the largest finite value
BUSTED_csubst_filtered <- BUSTED_csubst %>%
  filter(!is.na(omegaCany2spe) & is.finite(omegaCany2spe))

BUSTED_csubst %>%
  filter(!is.na(omegaCany2spe) & is.finite(omegaCany2spe))

max_finite_value <- max(BUSTED_csubst_filtered$omegaCany2spe, na.rm = TRUE, finite = TRUE)
max_finite_value

# Replace infinite values with the largest finite value
BUSTED_csubst_filtered <- BUSTED_csubst %>%
  mutate(omegaCany2spe = ifelse(omegaCany2spe > max_finite_value, max_finite_value, omegaCany2spe))
max(BUSTED_csubst_filtered$omegaCany2spe)

desired_order <- c("NoDiffSelPress", "DiffSelPress")  # Replace with your desired order

# Convert the 'Type' variable to a factor with the desired order
BUSTED_csubst_filtered$Type <- factor(BUSTED_csubst_filtered$Type, levels = desired_order)
my_colors <- c("#a5c4cf","#c82b27")

ggplot(BUSTED_csubst_filtered, aes(x=Type, y=omegaCany2spe, fill=Type)) +
  geom_violin(alpha = 0.3, width = 0.6, color = NA) +
  scale_fill_manual(values = my_colors) +
  #geom_jitter(aes(colour = Type),alpha = 0.5, size = 0.5, width = 0.3) +
  scale_color_manual(values = my_colors) +
  geom_boxplot(width = 0.005, outlier.shape = NA, size = 0.1) +
  a_plex_theme(base_size=12,axis_title_size=12) +
  #scale_y_continuous(trans='log10') +
  #geom_text_repel(data=subst_DiffSelPress, aes(label = paste(OG,'|',GeneSymbol, sep = '')),
  #               point.padding = 0.1, box.padding = 1, max.overlaps = 100,
  #              segment.color = 'grey50', size = 2, color='black') +
  a_plex_theme(base_size=12,axis_title_size=12) + ylim(1,500) 

wilcox.test(subst_NoDiffSelPress$omegaCany2spe, subst_DiffSelPress$omegaCany2spe, alternative='less')
median(subst_NoDiffSelPress$omegaCany2spe)
median(subst_DiffSelPress$omegaCany2spe)
head(BUSTED_csubst_filtered)

ggplot(BUSTED_csubst_filtered, aes(x=Type, y=dNCany2spe, fill=Type)) +
  geom_violin(alpha = 0.3, width = 0.6, color = NA) +
  scale_fill_manual(values = my_colors) +
  #geom_jitter(aes(colour = Type),alpha = 0.5, size = 0.5, width = 0.3) +
  scale_color_manual(values = my_colors) +
  geom_boxplot(width = 0.005, outlier.shape = NA, size = 0.1) +
  a_plex_theme(base_size=12,axis_title_size=12) +
  #scale_y_continuous(trans='log10') +
  #geom_text_repel(data=subst_DiffSelPress, aes(label = paste(OG,'|',GeneSymbol, sep = '')),
  #               point.padding = 0.1, box.padding = 1, max.overlaps = 100,
  #              segment.color = 'grey50', size = 2, color='black') +
  a_plex_theme(base_size=12,axis_title_size=12) + ylim(1,250) 

median(subst_NoDiffSelPress$dNCany2spe)
median(subst_DiffSelPress$dNCany2spe)

wilcox.test(subst_NoDiffSelPress$dNCany2spe, subst_DiffSelPress$dNCany2spe, alternative='less')


ggplot(BUSTED_csubst_filtered, aes(x=Type, y=dSCany2spe, fill=Type)) +
  geom_violin(alpha = 0.3, width = 0.6, color = NA) +
  scale_fill_manual(values = my_colors) +
  #geom_jitter(aes(colour = Type),alpha = 0.5, size = 0.5, width = 0.3) +
  scale_color_manual(values = my_colors) +
  geom_boxplot(width = 0.005, outlier.shape = NA, size = 0.1) +
  a_plex_theme(base_size=12,axis_title_size=12) +
  #scale_y_continuous(trans='log10') +
  #geom_text_repel(data=subst_DiffSelPress, aes(label = paste(OG,'|',GeneSymbol, sep = '')),
  #               point.padding = 0.1, box.padding = 1, max.overlaps = 100,
  #              segment.color = 'grey50', size = 2, color='black') +
  ylim(0,10) 

median(subst_NoDiffSelPress$dSCany2spe)
median(subst_DiffSelPress$dSCany2spe)

wilcox.test(subst_NoDiffSelPress$dSCany2spe, subst_DiffSelPress$dSCany2spe, alternative='less')

# synonymous rates of 9 types of combinatorial substitutions
ggplot(BUSTED_extended_filtered, aes(x=Type, y=dSCany2spe)) + 
  geom_violin(aes(color = Type)) +
  geom_jitter(aes(color = Type)) +
  a_plex_theme(base_size=12,axis_title_size=12) +
  #scale_y_continuous(trans='log10') +
  #geom_text_repel(data=BUSTED_extended, aes(label = paste(OG,'|',GeneSymbol, sep = '')),
  #                point.padding = 0.1, box.padding = 1, max.overlaps = 100,
  #                segment.color = 'grey50', size = 2, color='black') +
  a_plex_theme(base_size=12,axis_title_size=12) + ylim(1,50) 


########################## RELAX #################################
RELAX<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/RELAX.Results.tsv", header=TRUE)
head(RELAX)

ordRELAX <- RELAX %>%
  arrange(Pvalue)

head(ordRELAX)

# Perform Bonferroni correction on Pvalue
ordRELAX <- ordRELAX %>%
  mutate(Pvalue_Bonferroni = pmin(1, p.adjust(Pvalue, method = "bonferroni")))

ordRELAX <- ordRELAX %>%
  mutate(Pvalue_BH = pmin(1, p.adjust(Pvalue, method = "BH")))

# Substitue 0 values with the smallest other value
ordRELAX <- ordRELAX %>%
  mutate(Pvalue_BH = if_else(Pvalue_BH == 0, min(Pvalue_BH[Pvalue_BH != 0]), Pvalue_BH))

ordRELAX <- ordRELAX %>%
  mutate(LogPvalue_BH = log10(Pvalue_BH))


# View the dataframe with the adjusted p-values
head(ordRELAX)
write.table(ordRELAX, file = "~/Documents/HeliconiiniProg/GregariusConvergence/RELAX.Results.Corrected.tsv", sep = "\t", row.names = FALSE)
length(ordRELAX$Pvalue)
RelaxedGenes <- ordRELAX %>% filter(Pvalue < 0.05 & K < 1) %>% pull(Gene)
length(RelaxedGenes)
IntensifiedGenes <- ordRELAX %>% filter(Pvalue < 0.05 & K > 1) %>% pull(Gene)
length(IntensifiedGenes)

ggplot(ordRELAX, aes(log(K,2), LogPvalue_BH)) +
  geom_point(data=subset(ordRELAX, Pvalue_BH < 0.05 & K > 1), color = c("#c82b27"),alpha=0.4) +  # Add points for each combination
  geom_point(data=subset(ordRELAX, Pvalue_BH < 0.05 & K < 1), color = c("#575b7b"),alpha=0.4) +  # Add points for each combination
  geom_point(data=subset(ordRELAX, Pvalue_BH > 0.05), color = c("#a5c4cf"),alpha=0.2) +  # Add points for each combination
  scale_y_reverse() +
  labs(x = "K in Gregarious branches", y = "Adjusted P-values (log10)") +
  geom_vline(xintercept = log(1,2), color = "black") +
  geom_hline(yintercept = -1.30103, linetype = "dashed", color = "red") +
  geom_text_repel(data=subset(ordRELAX, Pvalue_BH < 0.05 & K > 1), aes(label = paste(Gene,'|',GeneSymbol, sep = '')),
                  point.padding = 0.1, box.padding = 0.1, max.overlaps = 30,
                  segment.color = 'grey50', size = 2) +
  xlim(-20,20) +
  a_plex_theme(base_size=12,axis_title_size=12) 

ggplot(ordRELAX, aes(log(K,2), LogPvalue_BH)) +
  geom_point(data=subset(ordRELAX, Pvalue_BH < 0.05 & K > 1), color = c("#c82b27"),alpha=0.4) +  # Add points for each combination
  geom_point(data=subset(ordRELAX, Pvalue_BH < 0.05 & K < 1), color = c("#575b7b"),alpha=0.4) +  # Add points for each combination
  geom_point(data=subset(ordRELAX, Pvalue_BH > 0.05), color = c("#a5c4cf"),alpha=0.2) +  # Add points for each combination
  scale_y_reverse() +
  labs(x = "K in Gregarious branches", y = "Adjusted P-values (log10)") +
  geom_vline(xintercept = log(1,2), color = "black") +
  geom_hline(yintercept = -1.30103, linetype = "dashed", color = "red") +
#  geom_text_repel(data=subset(ordRELAX, Pvalue_BH < 0.05 & K > 1), aes(label = paste(Gene,'|',GeneSymbol, sep = '')),
 #                 point.padding = 0.1, box.padding = 0.1, max.overlaps = 30,
  #                segment.color = 'grey50', size = 2) +
  xlim(-20,20) +
  a_plex_theme(base_size=12,axis_title_size=12) 

########################## RELAX + CSUBST #################################
csubst<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/All.csubst_cb_2_OCNa2s1.0omegaCa2s1.0_filtered.tsv", header=TRUE)
head(csubst)

length((csubst %>% distinct(OG))$OG)

ggplot(csubst, aes(x = omegaCspe2any)) +
  geom_density() +
  xlim(0,5) +
  theme_minimal()
head(ordRELAX)
# Filter df_B for Pvalue_BH < 0.05
Intendified_genes <- ordRELAX %>%
  filter(Pvalue < 0.05 & K > 1) %>%
  pull(Gene)

length(Intendified_genes)

Relaxed_genes <- ordRELAX %>%
  filter(Pvalue < 0.05 & K < 1) %>%
  pull(Gene)

length(Relaxed_genes)

Neutral_genes <- ordRELAX %>%
  filter(Pvalue > 0.05) %>%
  pull(Gene)

length(Neutral_genes)

# Select data in df_A matching the filtered genes
subst_Intensified <- subst %>%
  filter(OG %in% Intendified_genes) %>%
  mutate(Type = "Intensified")

head(subst_Intensified)

subst_Relaxed <- subst %>%
  filter(OG %in% Relaxed_genes) %>%
  mutate(Type = "Relaxed")

head(subst_Relaxed)

subst_Neutral <- subst %>%
  filter(OG %in% Neutral_genes) %>%
  mutate(Type = "NoShift")


head(subst_Neutral)

# Concatenate dataframes
subst_extended <- bind_rows(subst_Intensified, subst_Relaxed, subst_Neutral)

(subst_extended %>% distinct(Type))$Type
head(subst_extended)

# Find the largest finite value
subst_extended_filtered <- subst_extended %>%
  filter(!is.na(omegaCany2spe) & is.finite(omegaCany2spe))

subst_extended %>%
  filter(!is.na(omegaCany2spe) & is.finite(omegaCany2spe))

max_finite_value <- max(subst_extended_filtered$omegaCany2spe, na.rm = TRUE, finite = TRUE)
max_finite_value

# Replace infinite values with the largest finite value
subst_extended_filtered <- subst_extended %>%
  mutate(omegaCany2spe = ifelse(omegaCany2spe > 1000, 1000, omegaCany2spe))
max(subst_extended_filtered$omegaCany2spe)

write.table(subst_extended_filtered, file = "~/Documents/HeliconiiniProg/GregariusConvergence/RELAX.Results.Corrected.Extended.tsv", sep = "\t", row.names = FALSE)

desired_order <- c("NoShift", "Relaxed", "Intensified")  # Replace with your desired order

# Convert the 'Type' variable to a factor with the desired order
subst_extended_filtered$Type <- factor(subst_extended_filtered$Type, levels = desired_order)


my_colors <- c("#a5c4cf","#575b7b","#c82b27")
?geom_boxplot
ggplot(subst_extended_filtered, aes(x=Type, y=omegaCany2spe, fill = Type)) +
  geom_violin(alpha = 0.3, width = 0.5, color = NA) +
  scale_fill_manual(values = my_colors) +
  geom_boxplot(width = 0.01, outlier.shape = NA, size = 0.1) +
  #geom_jitter(aes(colour = Type),alpha = 0.5, size = 0.5, width = 0.3) +
  scale_color_manual(values = my_colors) +
  a_plex_theme(base_size=12,axis_title_size=12) +
#  geom_text_repel(data=subst_Intensified, aes(label = paste(OG,'|',GeneSymbol, sep = '')),
 #                 point.padding = 0.1, box.padding = 1, max.overlaps = 100,
  #                segment.color = 'grey50', size = 2, color='black') +
  a_plex_theme(base_size=12,axis_title_size=12) + ylim(1,150) 


subst_Neutral <- subst_extended_filtered %>%
  filter(Type == "NoShift")

subst_Relaxed <- subst_extended_filtered %>%
  filter(Type == "Relaxed")

subst_Intensified <- subst_extended_filtered %>%
  filter(Type == "Intensified")

median(subst_Neutral$omegaCany2spe)
mean(subst_Neutral$omegaCany2spe)

median(subst_Relaxed$omegaCany2spe)
mean(subst_Relaxed$omegaCany2spe)

median(subst_Intensified$omegaCany2spe)
mean(subst_Intensified$omegaCany2spe)

# first sample is stochastically greater than the distribution of values in the second sample
wilcox.test(subst_Neutral$omegaCany2spe, subst_Relaxed$omegaCany2spe, alternative='greater')
# first sample is stochastically lesser than the distribution of values in the second sample
wilcox.test(subst_Neutral$omegaCany2spe, subst_Intensified$omegaCany2spe, alternative='less')
# first sample is stochastically greater than the distribution of values in the second sample
wilcox.test(subst_Intensified$omegaCany2spe,subst_Relaxed$omegaCany2spe, alternative='greater')

ggplot(subst_extended_filtered, aes(x=Type, y=dNCany2spe, fill = Type)) +
  geom_violin(alpha = 0.3, width = 0.5, color = NA) +
  scale_fill_manual(values = my_colors) +
  geom_boxplot(width = 0.01, outlier.shape = NA, size = 0.1) +
  #geom_jitter(aes(colour = Type),alpha = 0.5, size = 0.5, width = 0.3) +
  scale_color_manual(values = my_colors) +
  a_plex_theme(base_size=12,axis_title_size=12) +
  #  geom_text_repel(data=subst_Intensified, aes(label = paste(OG,'|',GeneSymbol, sep = '')),
  #                 point.padding = 0.1, box.padding = 1, max.overlaps = 100,
  #                segment.color = 'grey50', size = 2, color='black') +
  a_plex_theme(base_size=12,axis_title_size=12) + ylim(0,100) 


median(subst_Neutral$dNCany2spe)
median(subst_Relaxed$dNCany2spe)
median(subst_Intensified$dNCany2spe)

# first sample is stochastically greater than the distribution of values in the second sample
wilcox.test(subst_Neutral$dNCany2spe, subst_Relaxed$dNCany2spe, alternative='greater')
# first sample is stochastically lesser than the distribution of values in the second sample
wilcox.test(subst_Neutral$dNCany2spe, subst_Intensified$dNCany2spe, alternative='less')
# first sample is stochastically greater than the distribution of values in the second sample
wilcox.test(subst_Intensified$dNCany2spe,subst_Relaxed$dNCany2spe, alternative='greater')

ggplot(subst_extended_filtered, aes(x=Type, y=dSCany2spe, fill = Type)) +
  geom_violin(alpha = 0.3, width = 0.5, color = NA) +
  scale_fill_manual(values = my_colors) +
  geom_boxplot(width = 0.01, outlier.shape = NA, size = 0.1) +
  #geom_jitter(aes(colour = Type),alpha = 0.5, size = 0.5, width = 0.3) +
  scale_color_manual(values = my_colors) +
  a_plex_theme(base_size=12,axis_title_size=12) +
  #  geom_text_repel(data=subst_Intensified, aes(label = paste(OG,'|',GeneSymbol, sep = '')),
  #                 point.padding = 0.1, box.padding = 1, max.overlaps = 100,
  #                segment.color = 'grey50', size = 2, color='black') +
  a_plex_theme(base_size=12,axis_title_size=12) + ylim(0,10) 


median(subst_Neutral$dSCany2spe)
median(subst_Relaxed$dSCany2spe)
median(subst_Intensified$dSCany2spe)

# first sample is stochastically greater than the distribution of values in the second sample
wilcox.test(subst_Neutral$dSCany2spe, subst_Relaxed$dSCany2spe, alternative='greater')
# first sample is stochastically lesser than the distribution of values in the second sample
wilcox.test(subst_Neutral$dSCany2spe, subst_Intensified$dSCany2spe, alternative='less')
# first sample is stochastically greater than the distribution of values in the second sample
wilcox.test(subst_Intensified$dSCany2spe,subst_Relaxed$dSCany2spe, alternative='greater')



set.seed(123)
x <- rnorm(20, mean = 1)
y <- rnorm(20, mean = 0)

# Perform Wilcoxon rank sum test
wilcox.test(x, y, alternative = "greater")

# Print the test result
print(result)



ggplot(subst_extended, aes(x=Type, y=dSCany2spe)) +
  geom_violin(aes(color = Type)) +
  geom_jitter(aes(color = Type)) +
  a_plex_theme(base_size=12,axis_title_size=12) +
  #scale_y_continuous(trans='log10') +
  geom_text_repel(data=subst_Intensified, aes(label = paste(OG,'|',GeneSymbol, sep = '')),
                  point.padding = 0.1, box.padding = 1, max.overlaps = 100,
                  segment.color = 'grey50', size = 2, color='black') +
  a_plex_theme(base_size=12,axis_title_size=12) + ylim(1,10) 


length(subst_Neutral$dSCany2spe)
length(subst_Relaxed$dSCany2spe)
length(subst_Intensified$dSCany2spe)


median(subst_Neutral$dSCany2spe)
median(subst_Relaxed$dSCany2spe)
median(subst_Intensified$dSCany2spe)

wilcox.test(subst_Neutral$dSCany2spe, subst_Relaxed$dSCany2spe, alternative='less')
wilcox.test(subst_Neutral$dSCany2spe, subst_Intensified$dSCany2spe, alternative='less')
wilcox.test(subst_Relaxed$dSCany2spe, subst_Intensified$dSCany2spe, alternative='less')


head(gene_SelecInGreg)
head(IntensifiedGenes)


Cross <- subst %>%
  filter(OG %in% gene_SelecInGreg) %>%
  filter(OG %in% IntensifiedGenes)

Cross

subst_SelecInGreg <- ordBUSTED %>%
  filter(Gene %in% gene_SelecInGreg) %>%
  mutate(Type = "SelecInGreg")

subst_SelecInGreg <- ordBUSTED %>%
  filter(Gene %in% gene_SelecInGreg) %>%
  mutate(Type = "SelecInGreg")


csubstMtx<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/Matrix.Zscore.tsv", header=TRUE, row.names = 1)
csubstMtx[is.na(csubstMtx) | is.nan(csubstMtx) | is.infinite(csubstMtx)] <- 0
print(class(csubstMtx))
str(csubstMtx)
# Convert to matrix
mat <- as.matrix(csubstMtx)
?heatmap
heatmap(mat, Rowv = NA, Colv = NA, scale = "none")

########################## CSUBST #####################################
CSUBST<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/All.csubstOnTree.Selected.tsv", header=TRUE)
head(CSUBST)

# Calculate the median values for each branch
median_values <- tapply(CSUBST$omegaCany2spe, CSUBST$branches, median)
class(median_values)

# Order branches by median values
sort(median_values, decreasing = TRUE)
ordered_branches <- names(sort(median_values, decreasing = TRUE))
ordered_branches
# Create a factor with ordered levels
CSUBST$branches <- factor(CSUBST$branches, levels = ordered_branches)

summary_stats <- data.frame(
  branches = unique(CSUBST$branches),
  median = sapply(unique(CSUBST$branches), function(branch) median(CSUBST$omegaCany2spe[CSUBST$branches == branch])),
  lower = sapply(unique(CSUBST$branches), function(branch) quantile(CSUBST$omegaCany2spe[CSUBST$branches == branch], 0.25)),
  upper = sapply(unique(CSUBST$branches), function(branch) quantile(CSUBST$omegaCany2spe[CSUBST$branches == branch], 0.75))
)

summary_stats
write.table(summary_stats, file = "~/Documents/HeliconiiniProg/GregariusConvergence/All.csubstOnTree.Selected.Medians.tsv", sep = "\t", row.names = FALSE)

ggplot(CSUBST, aes(x = branches, y = omegaCany2spe)) +
  geom_jitter(size = 5, alpha = 0.3, color = "#e3782f", width = 0.2) +
  a_plex_theme(base_size = 12, axis_title_size = 12) +
  geom_crossbar(data = summary_stats, aes(x = branches, y = median, ymin = median, ymax = median),
                width = 1, color = "#e3782f", alpha = 0.5, position = position_dodge(width = 0.75)) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels

mediansCSUBST<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/All.csubstOnTree.Selected.Medians.Edit.tsv", header=TRUE)
mediansCSUBST

mediansCSUBST$branches <- factor(mediansCSUBST$branches, levels = ordered_branches)
mediansCSUBST

ggplot(mediansCSUBST, aes(x=branches, y=median, color=median)) + 
  scale_color_gradient(low = "#f1bc98", high = "#984a14") +
  geom_point(size = 5) + scale_y_log10() +
  a_plex_theme(base_size = 12, axis_title_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  

GenePerLineage<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/GenesPerLineagesPlot.tsv", header=TRUE)
head(GenePerLineage)

LineageOrder <- c('HburHwal','Hwal','Hbur','HdorHhebHxan','HdorHxan','Hdor',
                  'Hxan','Hheb','Haoe','HantHconHdemHeleHertHhewHleuHricHsapHsar',
                  'HantHconHeleHhewHleuHsapHsar','HconHeleHhewHleuHsapHsar',
                  'HantHconHeleHhewHsap','HhewHsap','Hsap','Hhew','HconHele',
                  'Hcon','Hele','Hant','HleuHsar','Hleu','Hric','Hsar','HdemHert',
                  'Hdem','Hert','Hhec','Evib','Djun')
GenePerLineage$Lineage <- factor(GenePerLineage$Lineage, levels = LineageOrder)

ggplot(GenePerLineage, aes(x=Lineage, y=Ngenes, fill=Ngenes)) + 
  geom_col() +
  scale_fill_gradient(low = "#a5c4cf", high = "#575b7b") +
  scale_y_log10() +
  a_plex_theme(base_size=12,axis_title_size=12) +
  a_plex_theme(base_size = 12, axis_title_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  


########################## UpSet plot #####################################
CSUBST<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/CSUBST.scOGs.txt", header=TRUE)

gene_Convergence <- subst_Convergence %>%
  filter(Type == "Convergence") %>%
  pull(Gene)

gene_Convergence<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/BUSTED-PH.scOGs.txt", header=FALSE)
RelaxedGenes<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/RELAXed.scOGs.txt", header=FALSE)
IntensifiedGenes<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/RELAXint.scOGs.txt", header=FALSE)
CSUBST<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/CSUBST.scOGs.txt", header=FALSE)


# If necessary, rename the column to "gene"
colnames(gene_Convergence)[1] <- "gene"  # or the correct index
colnames(RelaxedGenes)[1] <- "gene"  # or the correct index
colnames(IntensifiedGenes)[1] <- "gene"  # or the correct index
colnames(CSUBST)[1] <- "gene"  # or the correct index
length(IntensifiedGenes$gene)
length(CSUBST$gene)
length(gene_Convergence$gene)
length(RelaxedGenes$gene)


combined <- purrr::reduce(list(data.frame(gene=gene_Convergence, BUSTEDPH=1),
                        data.frame(gene=RelaxedGenes, Relaxed=1),
                        data.frame(gene=IntensifiedGenes, Intensified=1),
                        data.frame(gene=CSUBST$gene, CSUBST=1)
), full_join)

tail(combined)
combined[is.na(combined)] <- 0
tail(combined)
write.table(combined, file = "~/Documents/HeliconiiniProg/GregariusConvergence/UpSetdata.tsv", sep = "\t", row.names = FALSE)
head(combined)

IntSize <- combined %>% 
  filter(Intensified == 0 & CSUBST == 1 & BUSTEDPH == 1 & Relaxed == 1) %>%
  pull(gene)
length(IntSize)

IntSize <- combined %>% 
  filter(Intensified == 1 & CSUBST == 1 & BUSTEDPH == 0 & Relaxed == 0) %>%
  pull(gene)
length(IntSize)

IntSize <- combined %>% 
  filter(Intensified == 1 & CSUBST == 0 & BUSTEDPH == 1 & Relaxed == 0) %>%
  pull(gene)
length(IntSize)

IntSize <- combined %>% 
  filter(Intensified == 0 & CSUBST == 1 & BUSTEDPH == 1 & Relaxed == 0) %>%
  pull(gene)
length(IntSize)

IntSize <- combined %>% 
  filter(Intensified == 0 & CSUBST == 1 & BUSTEDPH == 1 & Relaxed == 1) %>%
  pull(gene)
length(IntSize)

IntSize <- combined %>% 
  filter(Intensified == 0 & CSUBST == 0 & BUSTEDPH == 1 & Relaxed == 1) %>%
  pull(gene)
length(IntSize)

IntSize <- combined %>% 
  filter(Intensified == 1 & CSUBST == 0 & BUSTEDPH == 0 & Relaxed == 0) %>%
  pull(gene)
length(IntSize)

IntSize <- combined %>% 
  filter(Intensified == 0 & CSUBST == 1 & BUSTEDPH == 0 & Relaxed == 0) %>%
  pull(gene)
length(IntSize)


head(combined)
?upset()
# Create UpSet plot
upset(combined, sets = c("BUSTEDPH", "Relaxed", "Intensified", "CSUBST"), order.by = "degree",
      sets.bar.color = c("#575b7b", "#c82b27", "#178b76","#e3782f"),
      main.bar.color = "#1c183d", mb.ratio = c(0.8, 0.2))


# Define set sizes
A <- 122  # Size of set A
B <- 192   # Size of set B
intersection <- 14  # Size of the intersection
total <- length(combined %>% 
                  filter(CSUBST == 1 | Intensified == 1) %>%
                  pull(gene))
# Calculate values for the contingency table
only_A <- A - intersection  # Elements in A but not in B
only_B <- B - intersection  # Elements in B but not in A
neither <- total - A - B + intersection  # Assuming a total population of 10,000 (adjust this as needed)

# Create the contingency table
contingency_table <- matrix(c(intersection, only_A, only_B, neither), 
                            nrow = 2, 
                            dimnames = list("A" = c("In B", "Not in B"),
                                            "B" = c("In A", "Not in A")))

contingency_table
# Perform Fisher's exact test
fisher_test_result <- fisher.test(contingency_table)

# Show the result
print(fisher_test_result)


# Define set sizes
A <- 122  # Size of set A
B <- 203   # Size of set B
intersection <- 25  # Size of the intersection
total <- length(combined %>% 
                  filter(BUSTEDPH == 1 | Intensified == 1) %>%
                  pull(gene))

# Calculate values for the contingency table
only_A <- A - intersection  # Elements in A but not in B
only_B <- B - intersection  # Elements in B but not in A
neither <- total - A - B + intersection  # Assuming a total population of 10,000 (adjust this as needed)

# Create the contingency table
contingency_table <- matrix(c(intersection, only_A, only_B, neither), 
                            nrow = 2, 
                            dimnames = list("A" = c("In B", "Not in B"),
                                            "B" = c("In A", "Not in A")))

contingency_table
# Perform Fisher's exact test
fisher_test_result <- fisher.test(contingency_table)

# Show the result
print(fisher_test_result)

# Define set sizes
A <- 192  # Size of set A
B <- 203   # Size of set B
intersection <- 2  # Size of the intersection
total <- length(combined %>% 
                  filter(BUSTEDPH == 1 | CSUBST == 1) %>%
                  pull(gene))
total
# Calculate values for the contingency table
only_A <- A - intersection  # Elements in A but not in B
only_B <- B - intersection  # Elements in B but not in A
neither <- 0 # <- total - A - B + intersection  # Assuming a total population of 10,000 (adjust this as needed)

# Create the contingency table
contingency_table <- matrix(c(intersection, only_A, only_B, neither), 
                            nrow = 2, 
                            dimnames = list("A" = c("In B", "Not in B"),
                                            "B" = c("In A", "Not in A")))

contingency_table
# Perform Fisher's exact test
fisher_test_result <- fisher.test(contingency_table)

# Show the result
print(fisher_test_result)



########################## PhyloP on Regulatory domains #################################
PhyloPRD<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/RegDomainPhyloPconaccOut.tsv", header=TRUE)

head(PhyloPRD)

lm_eqn <- function(OG){
  m <- lm(-NonGregarious_AccAverage ~ -Gregarious_AccAverage, PhyloPRD);
  eq <- substitute(italic(x) == a + b %.% italic(y)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

lm_eqn(OG)

fit=lm(-NonGregarious_AccAverage~-Gregarious_AccAverage,data=PhyloPRD)
summary(fit)



AccAverage <- ggplot(PhyloPRD, aes(x = -Gregarious_AccAverage, y = -NonGregarious_AccAverage, size=Nclades))+
  geom_point(alpha=0.2, stroke = 0) + geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") #+
  xlim(0.5,0.8) + ylim(0.5,0.8)
cor(-PhyloPRD$Gregarious_AccAverage, -PhyloPRD$NonGregarious_AccAverage, method = "pearson")
AccAverage

ggplot(PhyloPRD, aes(x = Gregarious_AccAverage, y = Gregarious_ConsAverage, size=Nclades))+
  geom_point(alpha=0.2, stroke = 0) + geom_smooth(method = "lm", se = FALSE, color = "red") #+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") #+
  xlim(-0.85,-0.55) + ylim(-0.85,-0.55)
cor(PhyloPRD$Gregarious_AccAverage, PhyloPRD$NonGregarious_AccAverage, method = "pearson")
AccAverage


AccMedian <- ggplot(PhyloPRD, aes(x = Gregarious_AccMedian, y = NonGregarious_AccMedian, size=Gregarious_NumSites))+
  geom_point(alpha=0.2, stroke = 0) + geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  xlim(-0.65,-0.3) + ylim(-0.65,-0.3) 
cor(PhyloPRD$Gregarious_AccMedian, PhyloPRD$NonGregarious_AccMedian, method = "pearson")
AccMedian

ConsAverage <- ggplot(PhyloPRD, aes(x = Gregarious_ConsAverage, y = NonGregarious_ConsAverage, size=Gregarious_NumSites))+
  geom_point(alpha=0.2, stroke = 0) + geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  xlim(0.3,0.55) + ylim(0.3,0.55)
cor(PhyloPRD$Gregarious_ConsAverage, PhyloPRD$NonGregarious_ConsAverage, method = "pearson")
ConsAverage

ConsMedian <- ggplot(PhyloPRD, aes(x = Gregarious_ConsMedian, y = NonGregarious_ConsMedian, size=Gregarious_NumSites))+
  geom_point(alpha=0.2, stroke = 0) + geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  xlim(0.25,0.55) + ylim(0.25,0.55)
cor(PhyloPRD$Gregarious_ConsMedian, PhyloPRD$NonGregarious_ConsMedian, method = "pearson")
ConsMedian

OverallAverage <- ggplot(PhyloPRD, aes(x = Gregarious_OverallAverage, y = NonGregarious_OverallAverage, size=Gregarious_NumSites))+
  geom_point(alpha=0.2, stroke = 0) + geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  xlim(-0.2,0.15) + ylim(-0.2,0.15)
cor(PhyloPRD$Gregarious_OverallAverage, PhyloPRD$NonGregarious_OverallAverage, method = "pearson")
OverallAverage

OverallMedian <- ggplot(PhyloPRD, aes(x = Gregarious_OverallMedian, y = NonGregarious_OverallMedian, size=Gregarious_NumSites))+
  geom_point(alpha=0.2, stroke = 0) + geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  xlim(-0.2,0.3) + ylim(-0.2,0.3)
cor(PhyloPRD$Gregarious_OverallMedian, PhyloPRD$NonGregarious_OverallMedian, method = "pearson")
OverallMedian

grid.arrange(AccAverage, AccMedian, ConsAverage, ConsMedian, OverallAverage, OverallMedian, nrow = 3, ncol = 2)

head(PhyloPRD)
length(PhyloPRD$Gene)

row_names <- PhyloPRD$Gene

# Assign the extracted column as row names to the dataframe
rownames(df) <- row_names

AccAverage_df <- subset(PhyloPRD, select = c('Gregarious_AccAverage', 'NonGregarious_AccAverage'))
# Assign the extracted column as row names to the dataframe
rownames(AccAverage_df) <- row_names
head(AccAverage_df)
length(AccAverage_df$Gregarious_AccAverage)

mean_vec <- colMeans(AccAverage_df)
mean_vec
cov_mat <- cov(AccAverage_df)
cov_mat

mahalanobis_dist <- mahalanobis(AccAverage_df, mean_vec, cov_mat)
head(mahalanobis_dist)

threshold <- qchisq(0.05, df = 2)
threshold

outliers <- AccAverage_df[mahalanobis_dist > threshold, ]
length(outliers$Gregarious_AccAverage)


mad <- function(x) {
  med <- median(x)
  mad_val <- median(abs(x- med))
  return(mad_val)
}

mad_X <- AccAverage_df$Gregarious_AccAverage
mad_Y <- AccAverage_df$NonGregarious_AccAverage

threshold_X <- 3 * mad_X
threshold_Y <- 3 * mad_Y

outliers <- AccAverage_df[abs(AccAverage_df$Gregarious_AccAverage - median(AccAverage_df$Gregarious_AccAverage)) > threshold_X
                          | (abs(AccAverage_df$NonGregarious_AccAverage - median(AccAverage_df$NonGregarious_AccAverage)) > threshold_Y), ]
length(outliers$Gregarious_AccAverage)



# Fit a linear model y = b0 + b1*x
model <- lm(NonGregarious_AccAverage ~ Gregarious_AccAverage, data = AccAverage_df)
#model <- lm(NonGregarious_AccAverage ~ poly(Gregarious_AccAverage, 2), data = AccAverage_df)
model
# Get the residuals from the linear model
residuals <- residuals(model)

# Set a threshold for identifying outliers
threshold <- 0.1  # Adjust this threshold as needed

# Identify indices of points where absolute residuals are above the threshold
outlier_indices <- which(abs(residuals) > threshold)

# Display identified outliers
cat("Outlier indices:", outlier_indices, "\n")

# Plot the data points with the linear model
plot(AccAverage_df$Gregarious_AccAverage, AccAverage_df$NonGregarious_AccAverage, pch = 16, col = "blue", main = "Identifying Outliers from Linear Model")
abline(model, col = "red")  # Add the fitted linear model line
# Highlight outliers on the plot
points(AccAverage_df$Gregarious_AccAverage[outlier_indices], AccAverage_df$NonGregarious_AccAverage[outlier_indices], pch = 16, col = "red")

# Generate new data for prediction
new_data <- data.frame(Gregarious_AccAverage = seq(min(AccAverage_df$Gregarious_AccAverage), max(AccAverage_df$Gregarious_AccAverage), length.out = 100))
new_data$Gregarious_AccAverage
# Predict using the model
predicted_values <- predict(model, newdata = new_data)


# Add fitted line
lines(new_data$Gregarious_AccAverage, predicted_values, col = "red", lwd = 2)

# Highlight outliers on the plot
points(AccAverage_df$Gregarious_AccAverage[outlier_indices], AccAverage_df$NonGregarious_AccAverage[outlier_indices], pch = 16, col = "red")
outlier_indices


#+
#  scale_x_log10() +
 # scale_y_log10() +
#  stat_smooth(data = PhyloPRD,
 #             method = "lm", se = TRUE) +
#  scale_color_manual(values = my_colors) + a_plex_theme(base_size=12,axis_title_size=12)
pmain


slope.FusedShort<-sma(AccAverage~ConsAverage*Group, data=PhyloPRD, log='xy', multcomp=FALSE)
summary(slope.FusedShort)
NoConvEle<-round(slope.FusedShort$coef$NoConvergence$`coef(SMA)`[1],2)
NoConvSlop<-round(slope.FusedShort$coef$NoConvergence$`coef(SMA)`[2],2)
NoConvR2<-round(slope.FusedShort$r2$NoConvergence[1],2)
NoConveq <- paste("y=", NoConvSlop, "*x", NoConvEle, " | R2: ",NoConvR2, sep='')
NoConveq
ConvEle<-round(slope.FusedShort$coef$Convergence$`coef(SMA)`[1],2)
ConvSlop<-round(slope.FusedShort$coef$Convergence$`coef(SMA)`[2],2)
ConvR2<-round(slope.FusedShort$r2$Convergence[1],2)
Conveq <- paste("y=", ConvSlop, "*x", ConvEle, " | R2: ",ConvR2, sep='')
Conveq
plot(slope.FusedShort)

# Fit linear models
subst_Convergence <- BUSTED_extended %>%
  filter(Type == "Convergence")

subst_NoConvergence <- BUSTED_extended %>%
  filter(Type == "NoConvergence")

lm_modelsConv <- lm(OmegaBackground ~ OmegaTest, data = subst_NoConvergence)
lm_modelsConv
coef(lm_modelsConv$coefficients[0])
# Extract coefficients and format equation text
eq_texts <- paste("Equation:", round(coef(lm_models)[2], 2), "* x +", round(coef(lm_models)[1], 2))

?wilcox.test
wilcox.test(subst_NoConvergence$OmegaBackground, subst_Convergence$OmegaBackground, alternative='greater')
median(subst_NoConvergence$OmegaBackground)
median(subst_Convergence$OmegaBackground)
wilcox.test(subst_NoConvergence$OmegaTest, subst_Convergence$OmegaTest, alternative='less')
median(subst_NoConvergence$OmegaTest)
median(subst_Convergence$OmegaTest)

p_values <- c(0.02376, 5.688e-10)
p.adjust(p_values, method = "bonferroni")

#desired_order <- c("DiffSelecInGreg","DiffNoSelecInBoth","DiffSelecInNoGreg","BothSelec","NoDiff")  # Replace with your desired order
#desired_order <- c("Convergent","DiffSelecInNoGreg","BothSelec","NoDiff")  # Replace with your desired order
desired_order <- c("NoConvergence","Convergence")  # Replace with your desired order

# Convert the 'Type' variable to a factor with the desired order
BUSTED_extended$Type <- factor(BUSTED_extended$Type, levels = desired_order)
#my_colors <- c("#ff1f1f", "#7b1727", "#514d72", "#170057","black")
#my_colors <- c("#ff1f1f", "#514d72", "#170057","black")
my_colors <- c("#a5c4cf","#c82b27")

pmain <- ggplot(BUSTED_extended, aes(x = OmegaTest, y = OmegaBackground, color = Type))+
  geom_point(alpha=0.4)+
  scale_x_log10() +
  scale_y_log10() +
  stat_smooth(data = BUSTED_extended, aes(x = OmegaTest, color = Type),
              method = "lm", se = TRUE) +
  geom_text(data = subset(BUSTED_extended, OmegaBackground < 0.01 & Type == "Convergence"), 
            aes(label = Gene), 
            hjust = -0.1, vjust = 0.1) + # Add text labels next to the linear model
  #geom_text(x = 0.01, y = 0.5, label = NoConveq, color = c("#a5c4cf")) +
  #geom_text(x = 1, y = 0.5, label = Conveq, color = c("#c82b27")) +
  scale_color_manual(values = my_colors) + a_plex_theme(base_size=12,axis_title_size=12)
pmain

subst_ConOBsum
subst_ConOBsum["Median"]
subst_NoConOBsum["Median"]
xdens <- axis_canvas(pmain, axis = "x") +
  geom_vline(xintercept = median(subst_Convergence$OmegaTest), color = "#c82b27") +
  geom_vline(xintercept = median(subst_NoConvergence$OmegaTest), color = "#a5c4cf") +
  geom_density(data = BUSTED_extended, aes(x = OmegaTest, fill = Type),
               alpha = 0.3, size = 0.2)+
  scale_x_log10() +
  scale_fill_manual(values = my_colors) + a_plex_theme(base_size=12,axis_title_size=12)
xdens
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
  geom_vline(xintercept = median(subst_Convergence$OmegaBackground), color = "#c82b27") +
  geom_vline(xintercept = median(subst_NoConvergence$OmegaBackground), color = "#a5c4cf") +
  geom_density(data = BUSTED_extended, aes(x = OmegaBackground, fill = Type),
               alpha = 0.3, size = 0.2)+
  coord_flip()+
  scale_x_log10() +
  scale_fill_manual(values = my_colors) + a_plex_theme(base_size=12,axis_title_size=12)
ydens
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)


########################## Volcano Plots #################################
AvDj_DEGs<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/CallumStuff/Data/Av_Dj_nonsig.csv", header=TRUE)

head(AvDj_DEGs)

AvDj_DEGs$Significance <- ifelse(AvDj_DEGs$adj.P.Val < 0.05 & AvDj_DEGs$logFC > 1, "Up-regulated",
                          ifelse(AvDj_DEGs$adj.P.Val < 0.05 & AvDj_DEGs$logFC < -1, "Down-regulated", "Not significant"))

head(AvDj_DEGs$Significance)

length(AvDj_DEGs %>% filter(Significance == "Up-regulated") %>% pull(OG))
length(AvDj_DEGs %>% filter(Significance == "Down-regulated") %>% pull(OG))


#15x10 landscape
ggplot(AvDj_DEGs, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(size = 2, stroke = 0, alpha = 0.7) +
  scale_color_manual(values = c("Up-regulated" = "#c82b27", "Down-regulated" = "#006999", "Not significant" = "#bfd4e6")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#6a7d8d") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#6a7d8d") +
  labs(title = "A. vanillae vs D. juno",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  xlim(-16, 16) +
  a_plex_theme(base_size=12,axis_title_size=12)

HhHd_DEGs<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/CallumStuff/Data/Hh_Hd_nonsig.csv", header=TRUE)

head(HhHd_DEGs)

HhHd_DEGs$Significance <- ifelse(HhHd_DEGs$adj.P.Val < 0.05 & HhHd_DEGs$logFC > 1, "Up-regulated",
                                 ifelse(HhHd_DEGs$adj.P.Val < 0.05 & HhHd_DEGs$logFC < -1, "Down-regulated", "Not significant"))

length(HhHd_DEGs %>% filter(Significance == "Up-regulated") %>% pull(OG))
length(HhHd_DEGs %>% filter(Significance == "Down-regulated") %>% pull(OG))

ggplot(HhHd_DEGs, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(size = 2, stroke = 0, alpha = 0.7) +
  scale_color_manual(values = c("Up-regulated" = "#c82b27", "Down-regulated" = "#006999", "Not significant" = "#bfd4e6")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#6a7d8d") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#6a7d8d") +
  labs(title = "H. hecale vs H. doris",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  xlim(-16, 16) +
  a_plex_theme(base_size=12,axis_title_size=12)

HeHs_DEGs<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/CallumStuff/Data/He_Hs_nonsig.csv", header=TRUE)

head(HeHs_DEGs)

HeHs_DEGs$Significance <- ifelse(HeHs_DEGs$adj.P.Val < 0.05 & HeHs_DEGs$logFC > 1, "Up-regulated",
                                 ifelse(HeHs_DEGs$adj.P.Val < 0.05 & HeHs_DEGs$logFC < -1, "Down-regulated", "Not significant"))
length(HeHs_DEGs %>% filter(Significance == "Up-regulated") %>% pull(OG))
length(HeHs_DEGs %>% filter(Significance == "Down-regulated") %>% pull(OG))

ggplot(HeHs_DEGs, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(size = 2, stroke = 0, alpha = 0.7) +
  scale_color_manual(values = c("Up-regulated" = "#c82b27", "Down-regulated" = "#006999", "Not significant" = "#bfd4e6")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#6a7d8d") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#6a7d8d") +
  labs(title = "H. erato vs H. sara",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  xlim(-16, 16) +
  a_plex_theme(base_size=12,axis_title_size=12)
