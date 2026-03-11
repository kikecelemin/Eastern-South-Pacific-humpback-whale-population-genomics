library(ggpubr)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(car)
install.packages("car")
setwd("C:/Users/ecele/OneDrive/Desktop/KIKE/04_Other_projects/Humpback/Analysis/09_ROH/Downsample")

roh_plink_5x <- read.table(file = "26_Humpback_whale_annotated_filtered.hom", header = T)
# prepare ROH .hom data output into run-length categories #############################################
# small
ROHsmall<-roh_plink_5x[roh_plink_5x$KB >= 100 & roh_plink_5x$KB <= 300,]
ROHsmall_cat<-aggregate(ROHsmall$KB ~ ROHsmall$IID, FUN = "sum")
colnames(ROHsmall_cat)<-c("Sample_ID", "sumKB")
ROHsmall_cat$Category<-"Short"
ROHsmall_cat$Category<-"0.1-0.5"
# medium
ROHmedium<-roh_plink_5x[roh_plink_5x$KB > 300 & roh_plink_5x$KB <= 1000,]
ROHmedium_cat<-aggregate(ROHmedium$KB ~ ROHmedium$IID, FUN = "sum")
colnames(ROHmedium_cat)<-c("Sample_ID", "sumKB")
ROHmedium_cat$Category<-"Medium"
ROHmedium_cat$Category<-"0.3-1"
# long
ROHlong<-roh_plink_5x[roh_plink_5x$KB > 1000,]
ROHlong_cat<-aggregate(ROHlong$KB ~ ROHlong$IID, FUN = "sum")
colnames(ROHlong_cat)<-c("Sample_ID", "sumKB")
ROHlong_cat$Category<-"Long"
ROHlong_cat$Category<-">1"
# write 3 categories in .csv files
write.csv(ROHsmall_cat,"ROH_small_category.csv")
write.csv(ROHmedium_cat,"ROH_medium_category.csv")
write.csv(ROHlong_cat,"ROH_long_category.csv")
# this was done once and now the csv file is modified
# subsequently, this only needs to be imported

# Plot ROH tracts per sample
data <- read.csv("ROH_input2.csv")

# Define the desired sample order
desired_order <- c(
  "MnEcu01", "MnEcu06", "MnEcu09", "MnEcu14", "MnEcu19", "MnEcu30", "MnEcu32", "MnEcu33",
  "MnEcu37", "MnEcu39", "MnMg136", "MnMg137", "MnMg139", "MnMg143", "MnMg146", "MnMg147",
  "MnMg149", "MnMg151", "MnMg152", "MnMg153", "MnMg154", "MnAnt02", "MnAnt03", "MnAnt04", 
  "MnAnt07", "MnAnt12"
)

# Reshape data to long format and calculate proportions
long_data <- data %>%
  pivot_longer(cols = Short_ROH:Long_ROH, 
               names_to = "ROH_Type", 
               values_to = "Mb_value") %>%
  mutate(Proportion = Mb_value / Sum_ROH) %>%
  mutate(ROH_Type = factor(ROH_Type, levels = c("Short_ROH", "Medium_ROH", "Long_ROH"))) %>%
  mutate(Sample = factor(Sample, levels = desired_order)) # Reorder samples

# Plot the vertical barplot
ggbarplot(long_data, 
          x = "Sample", 
          y = "Mb_value", 
          fill = "ROH_Type", 
          color = "black", 
          palette = "Set2", 
          position = position_stack(reverse = TRUE)) +
  theme_classic() +
  labs(x = "Sample",
       y = "Sum ROH (Mb)",
       fill = "ROH Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ scale_y_continuous(limits = c(0,250), expand = c(0, 0))

# Plot per region
FROH <- read.csv("ROH_input3.csv")
FROH$Region<- fct_relevel(FROH$Region, "ALL", "ECU", "MAG", "ANT")

p<-ggplot(FROH,aes(x=Region,y=FROH,fill=Region))
p.box<-p+geom_boxplot(notch=FALSE, outlier.shape = NA)+geom_jitter(position=position_jitter(width=.1, height=0))+ scale_fill_manual(values = c("MAG"="#11823b", "ANT"="#a70000", "ECU"="#dab600", "ALL"="#82AEC0"))+ theme_classic()+theme(legend.position = "none")
p.box

p2<-ggplot(FROH,aes(x=Region,y=F,fill=Region))
p2.box<-p2+geom_boxplot(notch=FALSE, outlier.shape = NA)+geom_jitter(position=position_jitter(width=.1, height=0))+ scale_fill_manual(values = c("MAG"="#11823b", "ANT"="#a70000", "ECU"="#dab600", "ALL"="#82AEC0"))+ theme_classic()+theme(legend.position = "none")
p2.box

# Plot per genetic cluster
FROH <- read.csv("ROH_input4.csv")
FROH$Region<- fct_relevel(FROH$Region, "Turquoise", "Purple")

anova_model_het <- aov(Heterozygosity ~ Region, data = FROH)
summary(anova_model_het)

anova_model_FROH <- aov(FROH ~ Region, data = FROH)
summary(anova_model_FROH)

anova_model_F <- aov(F ~ Region, data = FROH)
summary(anova_model_F)

#FROH
p<-ggplot(FROH,aes(x=Region,y=FROH,fill=Region))
p.box<-p+geom_boxplot(notch=FALSE, outlier.shape = NA)+geom_jitter(position=position_jitter(width=.1, height=0))+ scale_fill_manual(values = c("Turquoise"="#008080", "Purple"="#800080"))+ theme_classic()+theme(legend.position = "none")
p.box

#F
p2<-ggplot(FROH,aes(x=Region,y=F,fill=Region))
p2.box<-p2+geom_boxplot(notch=FALSE, outlier.shape = NA)+geom_jitter(position=position_jitter(width=.1, height=0))+ scale_fill_manual(values = c("Turquoise"="#008080", "Purple"="#800080"))+ theme_classic()+theme(legend.position = "none")
p2.box

#Heterozygosity
p3<-ggplot(FROH,aes(x=Region,y=Heterozygosity,fill=Region))
p3.box<-p3+geom_boxplot(notch=FALSE, outlier.shape = NA)+geom_jitter(position=position_jitter(width=.1, height=0))+ scale_fill_manual(values = c("Turquoise"="#008080", "Purple"="#800080"))+ theme_classic()+theme(legend.position = "none")
p3.box

###################Correlations

ggplot(data= FROH, aes(x=Heterozygosity,y=FROH))+ geom_point(color="#82AEC0")+ theme_classic()+geom_smooth(method=lm,color="#a70000")
lm_Het_FROH=lm(Heterozygosity~FROH, data = FROH)
summary(lm_Het_FROH)

ggplot(data= FROH, aes(x=Heterozygosity,y=F))+ geom_point(color="#82AEC0")+ theme_classic()+geom_smooth(method=lm,color="#a70000")
lm_Het_F=lm(Heterozygosity~F, data = FROH)
summary(lm_Het_F)

ggplot(data= FROH, aes(x=F,y=FROH))+ geom_point(color="#82AEC0")+ theme_classic()+geom_smooth(method=lm,color="#a70000")
lm_F_FROH=lm(F~FROH, data = FROH)
summary(lm_F_FROH)

###################ANOVAs
FROH <- read.csv("ROH_input5.csv")
FROH$Region<- fct_relevel(FROH$Region, "ECU", "MAG", "ANT")

anova_model_het <- aov(Heterozygosity ~ Region, data = FROH)
summary(anova_model_het)
tukey_het<-TukeyHSD(anova_model_het)
tukey_het

anova_model_FROH <- aov(FROH ~ Region, data = FROH)
summary(anova_model_FROH)
tukey_FROH<-TukeyHSD(anova_model_FROH)
tukey_FROH

anova_model_F <- aov(F ~ Region, data = FROH)
summary(anova_model_F)
tukey_F<-TukeyHSD(anova_model_F)
tukey_F

