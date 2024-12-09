
#BiocManager::install("ggHex")
#install.packages("ggHex")

library(seqinr)
library(karyoploteR)
library(ggplot2)
library(awtools)
library(dplyr)
library(cowplot)
library(ggHex)
library(viridis)
library(ggbeeswarm)
library(ggrepel)
library(qqman)

########################## BUSTED-PH + PhyloAcc #################################
BustPhyloACC<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/Elem_lik.AllscOGs.M1.StrictGregarious.tsv", header=TRUE)
head(BustPhyloACC)
BustPhyloACC$Condition

BustPhyloACCcneeDiff <- BustPhyloACC %>%
  filter(Condition == "Diff") %>%
  pull(originalId)

BustPhyloACCDiff <- BustPhyloACC %>%
  filter(originalId %in% BustPhyloACCcneeDiff)
#
BustPhyloACCcneeNoDiff <- BustPhyloACC %>%
  filter(Condition == "NoDiff") %>%
  pull(originalId)

BustPhyloACCNoDiff <- BustPhyloACC %>%
  filter(originalId %in% BustPhyloACCcneeNoDiff)

head(BustPhyloACCcneeDiff)

wilcox.test(BustPhyloACCNoDiff$AccelRate, BustPhyloACCDiff$AccelRate, alternative='less')
median(BustPhyloACCNoDiff$AccelRate)
median(BustPhyloACCDiff$AccelRate)

wilcox.test(BustPhyloACCNoDiff$ConservedRate, BustPhyloACCDiff$ConservedRate, alternative='greater')
median(BustPhyloACCNoDiff$ConservedRate)
median(BustPhyloACCDiff$ConservedRate)



BustPhyloACC <- bind_rows(BustPhyloACCNoDiff, BustPhyloACCDiff)

BustPhyloACC$Condition <- factor(BustPhyloACC$Condition, levels = c("NoDiff","Diff"))

my_colors <- c("#a5c4cf","#c82b27")

ggplot(BustPhyloACC, aes(x=Condition, y=AccelRate, fill=Condition)) +
  geom_violin(alpha = 0.5, width = 0.5, color = NA) +
  scale_fill_manual(values = my_colors) +
  #geom_jitter(aes(colour = Condition),alpha = 0.5, size = 0.5, width = 0.3) +
  scale_color_manual(values = my_colors) +
  geom_boxplot(width = 0.03, outlier.shape = NA, size = 0.1) +
  a_plex_theme(base_size=12,axis_title_size=12) #+ ylim(1,550) 

ggplot(BustPhyloACC, aes(x=Condition, y=ConservedRate, fill=Condition)) +
  geom_violin(alpha = 0.5, width = 0.5, color = NA) +
  scale_fill_manual(values = my_colors) +
  #geom_jitter(aes(colour = Condition),alpha = 0.5, size = 0.5, width = 0.3) +
  scale_color_manual(values = my_colors) +
  geom_boxplot(width = 0.03, outlier.shape = NA, size = 0.1) +
  a_plex_theme(base_size=12,axis_title_size=12) #+ ylim(0.02,0.25) 


pmain <- ggplot(BustPhyloACC, aes(x = AccelRate, y = ConservedRate, color = Condition))+
  geom_point(alpha=0.4)+
  scale_x_log10() +
  scale_y_log10() +
  stat_smooth(data = BustPhyloACC, aes(x = AccelRate, color = Condition),
              method = "lm", se = TRUE) +
  #geom_text(data = subset(BustPhyloACC, ConservedRate < 0.01 & Condition == "Diff"), 
   #         aes(label = Gene), 
    #        hjust = -0.1, vjust = 0.1) + # Add text labels next to the linear model
  #geom_text(x = 0.01, y = 0.5, label = NoConveq, color = c("#a5c4cf")) +
  #geom_text(x = 1, y = 0.5, label = Conveq, color = c("#c82b27")) +
  scale_color_manual(values = my_colors) + a_plex_theme(base_size=12,axis_title_size=12)
pmain

slope.FusedShort<-sma(ConservedRate~AccelRate*Condition, data=BustPhyloACC, log='xy', multcomp=FALSE)
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
par(mar=c(1, 1, 1, 1) + 0.1)
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


########################## Gene Enrichment + PhyloAcc #################################
GenePhyloACC<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/GeneEnrichment.aCNEE.tsv", header=TRUE)
head(GenePhyloACC)

GenePhyloACC <- GenePhyloACC %>%
  filter(nCNEEs > 3)
  
head(GenePhyloACC)

ordGenePhyloACC <- GenePhyloACC %>%
  arrange(Pval)

head(ordGenePhyloACC)

# Perform Bonferroni correction on Pvalue
ordGenePhyloACC <- ordGenePhyloACC %>%
  mutate(Pvalue_Bonferroni = pmin(1, p.adjust(Pval, method = "bonferroni")))

ordGenePhyloACC <- ordGenePhyloACC %>%
  mutate(Pvalue_BH = pmin(1, p.adjust(Pval, method = "BH")))

# Substitue 0 values with the smallest other value
ordGenePhyloACC <- ordGenePhyloACC %>%
  mutate(Pvalue_BH = if_else(Pvalue_BH == 0, min(Pvalue_BH[Pvalue_BH != 0]), Pvalue_BH))

ordGenePhyloACC <- ordGenePhyloACC %>%
  mutate(LogPvalue_BH = log10(Pvalue_BH))


# View the dataframe with the adjusted p-values
head(ordGenePhyloACC)
write.table(ordGenePhyloACC, file = "~/Documents/HeliconiiniProg/GregariusConvergence/GeneEnrichment.aCNEE.3More.Corrected.tsv", sep = "\t", row.names = FALSE)

BestGenePhyloACC <- ordGenePhyloACC %>%
  filter(Pvalue_BH < 0.05)

BestGenePhyloACC <- as.data.frame(BestGenePhyloACC)
BestGenePhyloACC

# Convert OG to factor and reorder levels by nCNEEs and Pvalue_BH
BestGenePhyloACC$OG <- factor(BestGenePhyloACC$OG, levels = unique(BestGenePhyloACC$OG[order(BestGenePhyloACC$nCNEEs, -BestGenePhyloACC$LogPvalue_BH)]))
  
# Create the dot plot
ggplot(BestGenePhyloACC, aes(x = OG, y = nCNEEs, size = -LogPvalue_BH, color = ObsFreq/ExpFreq)) +
  geom_point() +  # Add points
  labs(title = "Gene Enrichment", x = "Gene", y = "nCNEEs") +  # Add labels
  scale_colour_gradientn( colours = c("#4888a2","#edb81d","#d19c2c","#d32f27")) +
  a_plex_theme(base_size=12,axis_title_size=12) +
  coord_flip()  # Apply a minimal theme and flip the coordinates
  
  
  
################ Print the alignment of a selected CNEE #######################################################
source("~/Documents/HeliconiiniProg/GregariusConvergence/drawAlign_function.R")
treeData <- prepare_data(tree_path = "~/Documents/HeliconiiniProg/GregariusConvergence/neut_Eisa_corrected.mod", 
                         species_name = "~/Documents/HeliconiiniProg/GregariusConvergence/species_names.txt", 
                         common_name = "~/Documents/HeliconiiniProg/GregariusConvergence/ButterfliesNames.txt")

targets = c("Hwal","Hbur","Hxan","Hdor","Hheb","Haoe","Hhec","Hsap","Hhew","Hcon","Hele","Hant","Hsar","Hleu",
            "Hric","Hdem","Hert","Evib","Djun")
conserved = c("Dple","Bany","Jcoe","Mcin","Smor","Pdid","Dpha","Diul","Ptel","Avfl","Avpe","Avcr",
              "Eisa","Elam","Eali","Elyb","Etal",
              "Htel","Hcly","Hhor","Hher","Hpet","Herd","Heet","Hlat","Hhim","Hper","Hcha",
              "Hhie",
              "Hege","Hnat","Hbes","Hism","Hnum","Heth","Hhel","Hatt","Hpar","Helv",
              "Hmel","Hcyd","Hpac","Htim","Hheu")
#bed <- read.delim("~/Dropbox/Heliconius_Project/CNEEanalysis/Data/Eisa.phyloacc.CNEE179234.M1.Heliconius.bed", header=F)
fasta <- read.alignment(file = "~/Documents/HeliconiiniProg/GregariusConvergence/CNEE109386471076.aln.fasta", format = "fasta")
fasta <- read.alignment(file = "~/Documents/HeliconiiniProg/GregariusConvergence/CNEE613611015333.aln.fasta", format = "fasta")
fasta <- read.alignment(file = "~/Documents/HeliconiiniProg/GregariusConvergence/CNEE613691015833.aln.fasta", format = "fasta")
fasta <- read.alignment(file = "~/Documents/HeliconiiniProg/GregariusConvergence/CNEE614331019565.aln.fasta", format = "fasta")


species_name = "~/Documents/HeliconiiniProg/GregariusConvergence/species_names.txt"
species <- read.table(species_name, header=F, sep="\t")
align <- as.matrix(fasta)
plotAlign <- function(align, treeData, target_species = NULL, legend = "top"){
  species <- treeData$tree$tip.label
  cols_sp <- rep(1, length(species))
  names(cols_sp) <- species
  if (!is.null(target_species)) cols_sp[target_species] <- 4
  #cols <- c("azure2", "#1e56b9", "#FFFFFF", "#ff3500")
  cols <- c("azure2", "#14a9dd", "#e2e2e2", "#f93838")
  
  element1 <- align[species, ]
  
  
  ele_cons <- apply(element1, 2, function(x) { 
    y <- xtabs(~ x); 
    cb <- names(y)[order(y, decreasing = TRUE)]
    cb <- cb[which(cb %in% c("a", "c", "g", "t", "A", "C", "G", "T"))[1]]
    z <- rep("substitution", length(x)) 
    z[x == cb] <- "consensus";
    z[x == "-"] <- "indel";
    z[x == "n" | x == "N" | x == "*"] <- "N"; #
    return(z);
  })
  
  rownames(ele_cons) <- treeData$tip
  dat_m <- melt(ele_cons) 
  dat_m$Var1 <- factor(dat_m$Var1, levels = rev(unique(dat_m$Var1)))
  
  p <- ggplot(dat_m, aes(as.factor(Var2), Var1)) + 
    geom_tile(aes(fill = factor(value, levels = c("N", "consensus", "indel", "substitution")))) + 
    scale_fill_manual(values = cols, name = "", drop = FALSE) +
    ylab(NULL) + xlab(paste(ncol(element1), "bp")) + 
    theme(panel.background = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 7, colour = rev(cols_sp), face = "italic"), 
          axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.title.x = element_text(size = rel(2.5)),
          axis.title.y = element_text(size = rel(2)),
          legend.text = element_text(size = 12), legend.position = legend,
          plot.background = element_rect(fill = "transparent")) 
  
  plot(p)
}
plotAlign(align, treeData, target_species=targets)
##############################################################################################################


########################### Annotations #########################################
bed_file <- "~/Documents/HeliconiiniProg/GregariusConvergence/Near.Sema1a.CNEEs.bed"
Chr <- c('Hmel210001o')
Start <- c(11249000)
End <- c(11339000)

#bed_file <- "~/Documents/HeliconiiniProg/GregariusConvergence/Near.Sema2a.CNEEs.bed"
#Chr <- c('Hmel217001o')
#Start <- c(6118000)
#End <- c(6139000)

bed_data <- read.delim(bed_file, header=FALSE, sep="\t")
bed_data
bed=data.frame(r=1:length(bed_data$V1), 
               x1=bed_data$V2, 
               x2=bed_data$V3, 
               y1=rep(1, each = length(bed_data$V1)),
               y2=rep(2, each = length(bed_data$V1)), 
               t=bed_data$V4)
length(bed$t)
bed

bed <- bed %>%
  filter(x1 > Start)

bed <- bed %>%
  filter(x2 < End)

par(mfrow = c(3, 2))

length(bed$t)

CNEEs <- ggplot() + 
  scale_x_continuous(name="x") + 
  scale_y_continuous(name="y") +
  geom_rect(bed, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="black", color="black", alpha=0.5) +
  a_plex_theme(base_size=6,axis_title_size=6) +
  xlim(Start,End) + theme(legend.position = "none")
CNEEs

bed_file <- "~/Documents/HeliconiiniProg/GregariusConvergence/Hmel.AccCNEEs.M1.sema-1a.bed"
#bed_file <- "~/Documents/HeliconiiniProg/GregariusConvergence/Hmel.AccCNEEs.M1.sema-2a.bed"
bed_data <- read.delim(bed_file, header=FALSE, sep="\t")
bed_data
bed=data.frame(r=1:length(bed_data$V1), 
               x1=bed_data$V2, 
               x2=bed_data$V3, 
               y1=rep(1, each = length(bed_data$V1)),
               y2=rep(2, each = length(bed_data$V1)), 
               t=bed_data$V4)
bed

aCNEEs <- ggplot() + 
  scale_x_continuous(name="x") + 
  scale_y_continuous(name="y") +
  geom_rect(bed, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color="black", alpha=0.5) +
  a_plex_theme(base_size=6,axis_title_size=6) +
  xlim(Start,End) + theme(legend.position = "none")
aCNEEs

gff_file <- "~/Documents/HeliconiiniProg/GregariusConvergence/Hmel.v3.2.annotation.CAT.gff3"
gff_data <- read.delim(gff_file, header=FALSE, sep="\t")
head(gff_data, n=5)

CDS <- subset(gff_data, (gff_data$V3 == "exon"))
head(CDS)
gff3=data.frame(r=1:length(CDS$V1),
                chr=CDS$V1,
                x1=CDS$V4, 
                x2=CDS$V5, 
                y1=rep(1, each = length(CDS$V1)),
                y2=rep(2, each = length(CDS$V1)), 
                t=CDS$V9)

head(gff3)
gff3 <- subset(gff3, (gff3$chr == Chr & gff3$x1 > Start & x2 < End))
gff3
genestr <- ggplot() + 
  scale_x_continuous(name="x") + 
  scale_y_continuous(name="y") +
  geom_rect(gff3, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="black", alpha=0.5) +
  a_plex_theme(base_size=6,axis_title_size=6) +
  xlim(Start,End)
genestr

library(rtracklayer)
library(Gviz)

# Set range you want from the bw file example)
range <- paste(Chr,':',Start,'-',End,sep='')
range

PhastCons <- import.bw(paste("~/Documents/HeliconiiniProg/GregariusConvergence/",Chr,".scores.bw", sep=""), which = GRanges(range))
PhyloP <- import.bw(paste("~/Documents/HeliconiiniProg/GregariusConvergence/Hmel_",Chr,".phyloP.bw", sep=""), which = GRanges(range))

head(PhyloP$score)
head(col_scale)

par(mfrow = c(1, 1))
PhyloP$score

#plot((start(PhastCons) + end(PhastCons)), PhastCons$score, type='h', col='#0099ce',lwd = 0.1, ylim = c(0,+1))
#plot((start(PhyloP) + end(PhyloP)), PhyloP$score, type='h', col = col_scale, ylim = c(-10,+10), lwd = 0.1)

# Generate color gradient based on the score values
color_gradient <- colorRampPalette(c("#aacfe9", "#70b4dc", "#008bbf", "#006293"))(length(PhastCons$score))

# Plot with a color gradient
plot((start(PhastCons) + end(PhastCons)),
     PhastCons$score, 
     type = 'h', 
     col = color_gradient[rank(PhastCons$score)], # Apply gradient based on rank of the score
     lwd = 0.1, 
     ylim = c(0, 1))

# Plot with ggplot2
PhastCons_df <- as.data.frame(PhastCons)
head(PhastCons_df)
PconsPlot <- ggplot(PhastCons_df, aes(x = start, y = score, color = score)) +
  geom_segment(aes(xend = start, yend = 0), size = 0.1) +  # Draw vertical lines from y = 0 to y = score
  scale_color_gradientn(colors = c("#aacfe9","#aacfe9","#aacfe9", "#aacfe9","#aacfe9","#70b4dc", "#008bbf", "#006293")) +  # Apply the color gradient
  ylim(0, 1) +  # Set y-axis limits
  labs(x = "Coordinates (start + end)", y = "PhyloP Score") +  # Axis labels
  theme_minimal()
PconsPlot

# Generate color gradient based on the score values
color_gradient <- colorRampPalette(c("#a70000","#d80221","#e8f2f9","#d9e9f5","#aacfe9", "#70b4dc", "#008bbf", "#006293"))(length(PhastCons$score))

# Rank the scores to map them to the gradient
color_mapped <- color_gradient[rank(PhyloP$score)]

# Plot with color gradient applied to each score
plot((start(PhyloP) + end(PhyloP)),
     PhyloP$score, 
     type = 'h', 
     col = color_mapped,   # Apply the color gradient based on the score
     ylim = c(-10, 10), 
     lwd = 0.1)

# Plot with ggplot2
PhyloP_df <- as.data.frame(PhyloP)
head(PhyloP_df)
PhyloPlot <- ggplot(PhyloP_df, aes(x = start, y = score, color = score)) +
  geom_segment(aes(xend = start, yend = 0), size = 0.1) +  # Draw vertical lines from y = 0 to y = score
  scale_color_gradientn(colors = c("#a70000","#d80221","#e8f2f9", "#008bbf", "#006293")) +  # Apply the color gradient
  ylim(-10, 10) +  # Set y-axis limits
  labs(x = "Coordinates (start + end)", y = "PhyloP Score") +  # Axis labels
  theme_minimal()  + theme(legend.position = "none")
PhyloPlot

ggarrange(CNEEs, aCNEEs, PhyloPlot, genestr,
          ncol = 1, nrow = 4)

############# Manhattan plot ################
ManPlotM1<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/Manhattan.M1.Strict.tsv", header=TRUE)
head(ManPlotM1)

par(mfrow = c(1, 1))

# Calculate average values for each chromosome
average_values <- aggregate(ManPlotM1$BF1, by=list(ManPlotM1$Chr), FUN=median)
average_values
# Plot the Manhattan plot
manhattan(ManPlotM1, chr="Chr", bp="BP", snp="NearestTSS", p="BF1",logp = FALSE,
          ylab = "Bayes factor between M1 and M0", col = c("#4888a2", "#edb81d"), suggestiveline = F, 
          genomewideline = FALSE, annotatePval = 50, cex = c(0.5), annotateTop = FALSE)#, ylim = c(0.1, 7), annotatePval = 5, annotateTop = FALSE)# cex = 0.4)

average_values <- aggregate(ManPlotM1$BF2, by=list(ManPlotM1$Chr), FUN=mean)
average_values
manhattan(ManPlotM1, chr="Chr", bp="BP", snp="NearestTSS", p="BF2",logp = FALSE,
          ylab = "Bayes factor between M1 and M2", col = c("#4888a2", "#edb81d"), suggestiveline = F, 
          genomewideline = FALSE, annotatePval = 20, cex = c(0.5), annotateTop = FALSE)#,  annotatePval = 5, annotateTop = FALSE)# cex = 0.4)

average_values <- aggregate(ManPlotM1$accelrate, by=list(ManPlotM1$Chr), FUN=median)
average_values
manhattan(ManPlotM1, chr="Chr", bp="BP", snp="NearestTSS", p="accelrate",logp = FALSE,
          ylab = "Acceleration Rate (M1)", col = c("#4888a2", "#edb81d"), suggestiveline = F, 
          genomewideline = FALSE, annotatePval = 3, cex = c(0.5), annotateTop = FALSE)#, ylim = c(0.1, 7), annotatePval = 5, annotateTop = FALSE)# cex = 0.4)


scatPlotM1<-read.table("~/Documents/HeliconiiniProg/GregariusConvergence/ChrACNEEstats.tsv", header=TRUE)
scatPlotM1$AccRate <- as.numeric(as.character(scatPlotM1$AccRate))
scatPlotM1

mean(scatPlotM1$MedianAccRate)

ggplot(scatPlotM1, aes(x = ChrLen, y = naCNEE, color = MedianAccRate, size = naCNEEperChr)) +
  geom_point(alpha = 1) +
  stat_smooth(aes(group = 1), method = "lm", se = FALSE) + 
  geom_text(aes(label = Chr), size = 3, color = "white") +  # You need to specify size for geom_text
  scale_colour_gradientn(colours = c("#4888a2","#edb81d","#d19c2c","#d32f27")) +
  a_plex_theme(base_size = 12, axis_title_size = 12)


model <- lm(naCNEE ~ ChrLen, data = scatPlotM1)
summary(model)$r.squared
cor(y=scatPlotM1$naCNEE, x=scatPlotM1$ChrLen, method = "pearson")

residuals <- residuals(model)
residuals

# Calculate mean and standard deviation of residuals
residual_mean <- mean(residuals)
residual_sd <- sd(residuals)

# Define threshold for outliers (e.g., 1.5 standard deviations)
threshold <- 1.5

# Identify outliers based on threshold
outliers <- which(abs(residuals - residual_mean) > threshold * residual_sd)

# View indices of outlier observations
outliers


############# Copy number plot ################


BiocManager::install("regioneR")


library(karyoploteR)
regions <- createRandomRegions(nregions=10000, length.mean = 1e6, mask=NA, non.overlapping=FALSE)
kp <- plotKaryotype()
kpPlotCoverage(kp, data=regions)


kfile<-'~/Dropbox/Heliconius_Project/ChrLocalHistory/EisaKaryotype.tsv'
tfile<-'~/Dropbox/Heliconius_Project/ChrLocalHistory/EisaWGA_TopologyDistributions.tsv'

kfile<-'~/Documents/HeliconiiniProg/ChrLocalHistory/HmelKaryotype.tsv'
tfile<-'~/Documents/HeliconiiniProg/ChrLocalHistory/Hmel221001o.TopologyDistributions.tsv'
CnRatiofile<-'~/Documents/HeliconiiniProg/GregariusConvergence/AllCNVratio.BedGraph.10k.5ks.mean.bed'
#mfile<-'~/Dropbox/Heliconius_Project/ChrLocalHistory/Eisa1Z.MAST.topologyDist.tsv'

custom.genome <- toGRanges(data.frame(read.table(kfile, header=TRUE, sep="\t")))

TopDist<-read.table(tfile, header=TRUE, sep="\t")
head(TopDist)
regs <- toGRanges(TopDist$Coord)
colors <- TopDist$color

kp <- plotKaryotype(genome = custom.genome, plot.type = 1)
kpPlotRegions(kp, data=regs, r0=-0.35, r1=0.1, col = colors, border=colors, lwd=1)


CNRatioDist<-read.table(CnRatiofile, header=TRUE, sep="\t")

kpPlotDensity(kp, data=CNRatioDist)

# Plot heatmap on chromosomes
?kpPlotHeatmap
kpPlotHeatmap(kp, data = CNRatioDist, chr = CNRatioDist$Chr,
              start = CNRatioDist$Start, end = CNRatioDist$End, values = CNRatioDist$CNratio,
              col = colorRampPalette(c("blue", "white", "red"))(50),
              border = NA)

kpAddBaseNumbers(kp, tick.dist = 1000000, tick.len = 7, tick.col="black", cex=1,
                 minor.tick.dist = 100000, minor.tick.len = 3, minor.tick.col = "gray")



blockfile<-'~/Dropbox/Heliconius_Project/ChrLocalHistory/BlockTopologyDist_new.tsv'
BlockTopDist<-read.table(blockfile, header=TRUE, sep="\t")
head(BlockTopDist)

Top1 <- subset(BlockTopDist, BlockTopDist$Topology == "Topology1_AoedeBasalMelpomene")
Top2 <- subset(BlockTopDist, BlockTopDist$Topology == "Topology2_Aoede+HeliconiusSS")
Top3 <- subset(BlockTopDist, BlockTopDist$Topology == "Topology3_AoedeNestedMelpomene")
Top4 <- subset(BlockTopDist, BlockTopDist$Topology == "Topology4_AoedeNestedEratoSaraSapho")
Top5 <- subset(BlockTopDist, BlockTopDist$Topology == "Topology5_AoedeNestedHeliconius")


ggplot(data=Top1, aes(x=BlockLength, y=Frequency)) +
  geom_bar(stat="identity")
ggplot(data=Top1, aes(x=BlockLength, y=Frequency)) +
  geom_bar(stat="identity",fill = "#2d60ad") + ylim(c(0,3000)) + xlim(c(1,350000)) +
  a_plex_theme(base_size=12,axis_title_size=12)
ggplot(data=Top2, aes(x=BlockLength, y=Frequency)) +
  geom_bar(stat="identity",fill = "#64293b") + ylim(c(0,3000)) + xlim(c(1,350000)) +
  a_plex_theme(base_size=12,axis_title_size=12)
ggplot(data=Top3, aes(x=BlockLength, y=Frequency)) +
  geom_bar(stat="identity",fill = "#448563") + ylim(c(0,3000)) + xlim(c(1,350000)) +
  a_plex_theme(base_size=12,axis_title_size=12)
ggplot(data=Top4, aes(x=BlockLength, y=Frequency)) +
  geom_bar(stat="identity",fill = "#7f802f") + ylim(c(0,3000)) + xlim(c(1,350000)) +
  a_plex_theme(base_size=12,axis_title_size=12)
ggplot(data=Top5, aes(x=BlockLength, y=Frequency)) +
  geom_bar(stat="identity",fill = "#595a5c") + ylim(c(0,3000)) + xlim(c(1,350000)) +
  a_plex_theme(base_size=12,axis_title_size=12)



############ N(x) distribution ############
file<-'~/Dropbox/Heliconius_Project/ChrLocalHistory/ChromosomeGammaDistributions.tsv'
GammaDIst <- read.table(file, header=TRUE, sep="\t")
head(GammaDIst)

ChrOrder=c('Eisa1Z00',
           'Eisa2100',
           'Eisa1800',
           'Eisa2000',
           'Eisa1700',
           'Eisa1900',
           'Eisa1600',
           'Eisa0700',
           'Eisa0500',
           'Eisa0800',
           'Eisa0300',
           'Eisa0200Eisa2700',
           'Eisa0200',
           'Eisa2700',
           'Eisa1000Eisa3100',
           'Eisa1000',
           'Eisa3100',
           'Eisa1200Eisa2800',
           'Eisa1200',
           'Eisa2800',
           'Eisa0600Eisa2200',
           'Eisa0600',
           'Eisa2200',
           'Eisa1300Eisa2500',
           'Eisa1300',
           'Eisa2500',
           'Eisa1100Eisa2600',
           'Eisa1100',
           'Eisa2600',
           'Eisa1500Eisa2900',
           'Eisa1500',
           'Eisa2900',
           'Eisa1400Eisa2400',
           'Eisa1400',
           'Eisa2400',
           'Eisa0400Eisa2300',
           'Eisa0400',
           'Eisa2300',
           'Eisa0900Eisa3000',
           'Eisa0900',
           'Eisa3000')

GammaDIst$Chr <- factor(GammaDIst$Chr,levels = ChrOrder)

df_median <- GammaDIst %>% group_by(Chr) %>% summarise(median = median(Gamma))
df_median

df2_agg <- df_median %>% 
  group_by(Chr) %>% 
  summarize(g_median = mean(median))

df3 <- merge(GammaDIst, df2_agg, by = "Chr", all.x = TRUE)
head(df3)
ggplot(df3, aes(x=Chr, y=Gamma, fill=g_median)) +
  geom_violin() +
  scale_fill_gradient2(low = "darkgreen", high = "darkred") +
  a_plex_theme(base_size=12,axis_title_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


file<-'~/Dropbox/Heliconius_Project/ChrLocalHistory/Astrals/AstralCFdistributionsPerChr.tsv'

CFdist <- read.table(file, header=TRUE, sep="\t")
head(CFdist)

df_median <- CFdist %>% group_by(Chr) %>% summarise(median = median(CF))
df_median

Med_agg <- df_median %>% 
  group_by(Chr) %>% 
  summarize(g_median = mean(median))

CFdist_agg <- merge(CFdist, Med_agg, by = "Chr", all.x = TRUE)
head(CFdist_agg)


ggplot(CFdist_agg, aes(x=Chr, y=log(CF), fill=g_median)) +
  geom_violin() +
  a_plex_theme(base_size=12,axis_title_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


NodeOrder=c('Heliconiinae','Heliconiini','PhDrPoDy','DioneAgraulisEueidesHeliconius',
            'DioneAgraulis','Agraulis','EueidesHeliconius','Eueides','Heliconius',
            'Erato','SaraSapho','EratoSaraSapho','AoeDorWalMelSil','DorWalMelSil','Doris',
            'WalMelSil','Wallacei','MelpomeneSilvaniformis')

NodeOrder=c('MelpomeneSilvaniformis','Wallacei','WalMelSil','Doris','DorWalMelSil',
            'AoeDorWalMelSil','SaraSapho','Erato','EratoSaraSapho','Heliconius',
            'Eueides','EueidesHeliconius','Agraulis','DioneAgraulis',
            'DioneAgraulisEueidesHeliconius','PhDrPoDy','Heliconiini','Heliconiinae')

CFdist$Node <- factor(CFdist$Node, levels = NodeOrder)

ggplot(CFdist, aes(x=Node, y=CF)) +
  geom_violin() +
  geom_beeswarm(aes(color=CF),size = 1, alpha = 0.5,dodge.width=.8,cex=0.7) +
  geom_label_repel(aes(label = ifelse(Chr=="Eisa1Z00",as.character(paste('Z ->',round(CF,2))),'')),
                   point.padding = 0.1, box.padding = 1, label.padding = 0.1,
                   segment.color = 'grey50', size = 2, max.overlaps = 50) +
  geom_label_repel(aes(label = ifelse(Chr=="AllAutosomes",as.character(paste('Autosomes ->',round(CF,2))),'')),
                   point.padding = 0.1, box.padding = 1, label.padding = 0.1,
                   segment.color = 'grey50', size = 2, max.overlaps = 50) +
  a_plex_theme(base_size=12,axis_title_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip()
