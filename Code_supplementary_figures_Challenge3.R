###############################################
######### SUPPLEMENTARY FIGURE CHALLENGE 3 ####
###############################################
## author: Francisco J. Pardo-Palacios, f.pardo.palacios@gmail.com
## author:Ana Conesa, ana.conesa@csic.es
## Last modified: March 9th 2023
###############################################

# Install and load packages
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(scales)
library(patchwork)
library(gridExtra)
library(grid)
library(RColorConesa)

#### set theme for plots
pub_theme <- theme_pubclean(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=13),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=13) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=10), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15.5)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  theme(legend.position = "bottom")

old.libplat.palette = c( "cDNA+ONT"="#74CDF0", "cDNA+PacBio"="#EE446F", "cDNA+Illumina"="#FFCF71", 
                         "CapTrap+ONT"="#7482F0", "R2C2+ONT"="#74F0D9", "dRNA+ONT"="#13BF5E", "CapTrap+PacBio"="#d14141")

libplat.palette = c( "cDNA-PacBio"="#c06636", "CapTrap-PacBio"="#802417", "cDNA-Illumina"="#e8b960",  "Freestyle-Freestyle"="#ce9344",
                     "cDNA-ONT"="#646e3b", "CapTrap-ONT"="#17486f", "R2C2-ONT"="#508ea2", "dRNA-ONT"="#2b5851"
)
cat.palette = c( "FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "GenicGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                 "Intergenic" = "darksalmon", "GenicIntron"="#41B6C4")


#Supplementary SX1 BUSCO Analysis Genome

BUSCO.data <- data.frame(BUSCO = c("Complete", "Complete_single", "Complete_Duplicated", "Framented", "Missing"), 
                    Value = c(9685, 9637, 48, 473, 1208), stringsAsFactors = FALSE)
BUSCO.data$Percentage <- round(BUSCO.data$Value/11366 * 100,1)

BU <- ggplot(BUSCO.data, aes(x=BUSCO, y=Percentage, fill = BUSCO)) +
               geom_bar(stat="identity") +
         pub_theme + 
         ylab("% BUSCOs") +
         ggtitle("BUSCO Analysis Manatee draft genome") +
         theme(plot.title = element_text(hjust = 0.5)) +
         theme(axis.text.x = element_text(angle=90)) 
BU
pdf("Supplementary_figure_SX31b.pdf")
annotate_figure(BU)
dev.off()

  #### Panel B: Overview of manatee and mouse submissions
# Manatee plots
# manatee_code <- read.csv("manatee/manatee_code.txt", sep=",", header = T )
# manatee_code$Lib_Plat <- apply(manatee_code, 1, function(x){
#   paste(x["Library_Preps"], x["Platform"], sep = "-")
# })
# manatee_code$Lib_DC=apply(cbind(manatee_code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
# manatee_code$Label <-apply(cbind(manatee_code[,c("Platform","Library_Preps", "Data_Category")]), 1, paste, collapse="-")

# A <- ggplot(manatee_code, aes(x=Tool, fill=Lib_Plat)) +
#   geom_bar(stat = "count", position = "stack") +
#   geom_text(aes(label = ..count.. , group=Lib_Plat), stat = "count", position = position_stack(vjust = 0.5) ,size = 4) +
#   pub_theme + 
#   scale_fill_manual(values = libplat.palette, limits=c("cDNA-ONT", "cDNA-PacBio","cDNA-Illumina")) +
#   ylab("# of submissions") +
#   ggtitle("Challenge 3 submissions for Manatee data") +
#   theme(plot.title = element_text(hjust = 0.5))


# ES_code <- read.csv("ES/ES_code.txt", sep=",", header = T )
# ES_code$Lib_Plat <- apply(ES_code, 1, function(x){
#   paste(x["Library_Preps"], x["Platform"], sep = "-")
# })
# ES_code$Lib_DC=apply(cbind(ES_code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
# ES_code$Label <-apply(cbind(ES_code[,c("Platform","Library_Preps", "Data_Category")]), 1, paste, collapse="-")
# 
# B <- ggplot(ES_code, aes(x=Tool, fill=Lib_Plat)) +
#   geom_bar(stat = "count", position = "stack") +
#   geom_text(aes(label = ..count.. , group=Lib_Plat), stat = "count", position = position_stack(vjust = 0.5) ,size = 4) +
#   pub_theme + 
#   scale_fill_manual(values = libplat.palette, limits=force ) +
#   ylab("# of submissions") +
#   theme(legend.title=element_blank())+
#   ggtitle("Challenge 3 submissions for Mouse ES data") +
#   theme(plot.title = element_text(hjust = 0.5))

# suppl = "SX2"
# figureSX2 <- ggarrange(A,B,
#                       labels = c( "a)", "b)"),
#                      ncol = 1, nrow = 2, common.legend = TRUE, legend="bottom") +
#   theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
# mylegend <- paste0("     Supplementary figure ", suppl, ". Submissions Challenge 3. a) Manatee, b) Mouse ES.")
# pdf("Supplementary_figure_SX2.pdf")
# annotate_figure(figureSX2,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
# dev.off()

## Supplementary SX3

mouse_ch3 <- read.csv("ES_challenge1/ES_challenge1_metrics.summary_table_SC.csv", sep=",", header = T )
mouse_ch1 <- read.csv("ES_challenge1/ES.summary_table_SC.csv", sep=",", header = T )
ES_code$Sample <- paste(ES_code$Alias, ES_code$Library_Preps, ES_code$Platform, ES_code$Data_Category, sep = "-")
code <- read.csv("ES_challenge1/code.csv", sep=",", header = T)
code$Sample <- paste(code$Alias, code$Library_Preps, code$Platform, code$Data_Category, sep = "-")
codes <- merge (ES_code, code, by.x = "Sample" , by.y = "Sample")
head(codes)
merged.data <- merge(codes, mouse_ch3, by.x = "pipelineCode.x", by.y = "ID")
head(mouse_ch3)

merged.data2 <- merge(merged.data, mouse_ch1, by.x = "pipelineCode.y", by.y = "ID")
A <- merged.data2[,c(3,20:28)]; B<- merged.data2[,c(3,30:38)]
colnames(A) = colnames(B) <-  colnames(mouse_ch3)[-2]
data.f <- as.data.frame(rbind (B,A))
data.f$Challenge <- c(rep("Challenge 1", nrow(A)), rep("Challenge 3", nrow(B)))
data.ff <- melt(data.f)
head(data.ff)
  
C <- ggplot(data.ff, aes(x = Challenge, y = value, fill = variable)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = cat.palette) +
  facet_wrap(~ID, ncol = 4, as.table = FALSE) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.position="bottom")

suppl = "SX3"
mylegend <- paste0("     Supplementary figure ", suppl, ". Detection of SQANTI categories for the same tools in Challenge 1 and 3. \n     Challenge 1 used the reference annotation and Challenge 3 did not. Ba = Bambu, IQ = IsoQuant.")
pdf("Supplementary_figure_SX3.pdf")
annotate_figure(C,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

