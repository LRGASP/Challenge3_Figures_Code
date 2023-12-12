###############################################
######### SUPPLEMENTARY FIGURE CHALLENGE 3 ####
###############################################
## author: Francisco J. Pardo-Palacios, f.pardo.palacios@gmail.com
## author:Ana Conesa, ana.conesa@csic.es
## Last modified: March 9th 2023
###############################################

# Install and load packages
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(scales)
library(patchwork)
library(gridExtra)
library(grid)
library(RColorConesa)
library(ggthemes)
library(dplyr)
library(data.table)
library(UpSetR)
library(stringr) # Load

outdir = "output/extended"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

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

#Extended Data Fig. 63b BUSCO Analysis Genome

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
pdf(paste0(outdir, "/Extended_Fig._63b.pdf"))
annotate_figure(BU)
dev.off()

## Extended Data Fig. 64
###################

# note: this generates labels in wrong order, from top to bottom they should be
#       Bambu, RNA_Bloom, rnaSPAdes, StringTie2_isoquant
# correct in edited slide, but not here

ES_mapping <- data.frame(
  Label = c("illumina_1", "ONT_1", "ONT_10", "ONT_11", "ONT_2", "ONT_3", "ONT_4", "ONT_5", "ONT_6", "ONT_7", "ONT_8", "ONT_9", "PB_1", "PB_2", "PB_3", "PB_4", "PB_5"),
  mapping_transcript = c(87.87, 100, 98.68, 99.98, 99.86, 100, 99.99, 99.88, 99.83, 99.02, 95.36, 99.99, NA, 99.98, 99.99, 99.99, 95.76),
  Lib_Plat = c("cDNA-Illumina", "CapTrap-ONT", "dRNA-ONT", "R2C2-ONT", "cDNA-ONT", "dRNA-ONT", "R2C2-ONT", "cDNA-ONT", "cDNA-ONT", "cDNA-ONT", "cDNA-ONT", "dRNA-ONT", "CapTrap-PacBio", "cDNA-PacBio", "CapTrap-PacBio", "cDNA-PacBio", "cDNA-PacBio"),
  Data_Category = c("SO", "LO", "LS", "LO", "LO", "LO", "LO", "LO", "LS", "LO", "LS", "LO", "LO", "LO", "LO", "LO", "LS"),
  Tool = c("rnaSPAdes", "Bambu", "rnaSPAdes", "StringTie2\nIsoQuant", "Bambu", "Bambu", "Bambu", "RNA\nBloom", "RNA\nBloom", "StringTie2\nIsoQuant", "rnaSPAdes", "StringTie2\nIsoQuant", "Bambu", "Bambu", "StringTie2\nIsoQuant", "StringTie2\nIsoQuant", "rnaSPAdes")
)

manatee_mapping <- data.frame(
  Label = c("illumina1", "ONT2", "ONT3", "ONT4", "ONT5", "PB1", "PB2", "PB3"),
  mapping_transcript = c(75.82, 99.95, 100, 98.6, 99.86, 99.74, NA, 98.08),
  Lib_Plat = c("cDNA-Illumina", "cDNA-ONT", "cDNA-ONT", "cDNA-ONT", "cDNA-ONT", "cDNA-PacBio", "cDNA-PacBio", "cDNA-PacBio"),
  Data_Category = c("SO", "LO", "LO", "LS", "LS", "LO", "LO", "LS"),
  Tool = c("rnaSPAdes", "RNA\nBloom", "StringTie2\nIsoQuant", "rnaSPAdes", "RNA\nBloom", "Bambu", "StringTie2\nIsoQuant", "rnaSPAdes")
)

ES_mapping$mapping_transcript[is.na(ES_mapping$mapping_transcript)] <- 98
manatee_mapping$mapping_transcript[is.na(manatee_mapping$mapping_transcript)] <- 98

SX3.1 <- ggplot(ES_mapping, aes(x=Label, y=as.numeric(mapping_transcript), color=Lib_Plat, shape=Data_Category))  + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(mapping_transcript), color=Lib_Plat), size=2) +
  geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
  pub_theme+
  theme_pubclean(flip=TRUE)+
  theme( axis.text.y  = element_blank(),
         axis.text.x = element_text(size=14),
         axis.ticks.y = element_blank()) +
  scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
  scale_color_manual(values =  libplat.palette, name="Library-Platform") +
  scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                            "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS"),
                   labels=c("SO", rep("LO", 6), rep("LS", 3)))+
  facet_grid(Tool ~., drop = TRUE, scales="free_y" ) +
  xlab("") + 
  theme(axis.title=element_text(size=16),
        strip.text.y = element_blank(), 
        axis.line.y = element_blank()) +
  theme(legend.position = "none") +
  scale_y_reverse("", sec.axis = sec_axis(~ . , breaks = NULL, name = "Mouse ES"),
                  label = unit_format(unit = "%"), limits=c(100 ,-4), expand = expansion(mult = c(0.1,0)))+
  coord_flip()

SX3.2 <- ggplot(manatee_mapping, aes(x=Label, y=as.numeric(mapping_transcript), color=Lib_Plat, shape=Data_Category)) + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(mapping_transcript), color=Lib_Plat), size=2) +
  geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
  pub_theme+
  theme_pubclean(flip=TRUE)+
  theme( axis.text.y  = element_blank(),
         axis.text.x = element_text(size=14),
         axis.ticks.y = element_blank()) +
  scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
  scale_color_manual(values =  libplat.palette, name="Library-Platform") +
  scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                            "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS"),
                   labels=c("SO", rep("LO", 6), rep("LS", 3)))+
  facet_grid(Tool ~., drop = TRUE, scales="free_y" ) +
  xlab("") + 
  scale_y_continuous("", sec.axis = sec_axis(~ . , breaks = NULL, name = "Manatee"), 
                     label = unit_format(unit = "%"),  limits=c(-4, 100), expand = expansion(mult = c(0,0.1)))+
  theme(strip.text.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.title=element_text(size=16)) +
  theme(legend.position = "none") +
  coord_flip() 

SX.mid3<- ggplot(ES_mapping,aes(x=1,y=Tool))+
  ggtitle("")+
  scale_x_continuous(expand=c(0,0),limits=c(1,1))+
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        #axis.text.y = element_text(vjust=0.5, hjust = 0.5, size=14),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA)) 
        #theme_wsj()+ scale_colour_wsj("colors6")

gg3.1 <- ggplot_gtable(ggplot_build(SX3.1))
gg3.2 <- ggplot_gtable(ggplot_build(SX3.2))
g.mid <- ggplot_gtable(ggplot_build(SX.mid3))

Ch3S2 <- grid.arrange(gg3.1,g.mid, gg3.2,ncol=3,widths=c(4/9,1/9,4/9),
                    top = textGrob("Mapping rate (%)",gp=gpar(fontsize=18,font=1)))

ggsave(file=paste0(outdir, "/Ch3S2.svg"), plot=Ch3S2, width=8, height=5)

suppl = "64"
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Mapping rate of transcript detected by Challenge 3 submissions.")
pdf(paste0(outdir, "/Extended_Fig._64.pdf"))
annotate_figure(Ch3S2,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()


## Extended Data Fig. 65
###################
ES_code <- read.csv("Challenge3_Figures_Data/ES_challenge1/ES_code.txt", sep=",", header = T )
mouse_ch3 <- read.csv("Challenge3_Figures_Data/ES_challenge1/ES_challenge1_metrics.summary_table_SC.csv", sep=",", header = T )
mouse_ch1 <- read.csv("Challenge3_Figures_Data/ES_challenge1/ES.summary_table_SC.csv", sep=",", header = T )
ES_code$Sample <- paste(ES_code$Alias, ES_code$Library_Preps, ES_code$Platform, ES_code$Data_Category, sep = "-")
code <- read.csv("Challenge3_Figures_Data/ES_challenge1/code.csv", sep=",", header = T)
code$Sample <- paste(code$Alias, code$Library_Preps, code$Platform, code$Data_Category, sep = "-")
codes <- merge (ES_code, code, by.x = "Sample" , by.y = "Sample")
merged.data <- merge(codes, mouse_ch3, by.x = "pipelineCode.x", by.y = "ID")

merged.data2 <- merge(merged.data, mouse_ch1, by.x = "pipelineCode.y", by.y = "ID")
A <- merged.data2[,c(3,17:25)]; B<- merged.data2[,c(3,27:35)]
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

ggsave(file=paste0(outdir, "/Ch3S3.svg"), plot=Ch3S2, width=8, height=5)

suppl = "65"
mylegend <- paste0("     Extended Data Fig. ", suppl, ". SQANTI category classification of transcript models detected by the same tools in Challenge 1 and 3.\n     Challenge 1 predictions used the reference annotation and Challenge 3 predictions did not.\n     Ba = Bambu, IQ = StringTie2/IsoQuant.")
pdf(paste0(outdir, "/Extended_Fig._65.pdf"))
annotate_figure(C,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig. 66
###################

ES_coding <- data.frame(
  Label = c("illumina_1", "ONT_1", "ONT_10", "ONT_11", "ONT_2", "ONT_3", "ONT_4", "ONT_5", "ONT_6", "ONT_7", "ONT_8", "ONT_9", "PB_1", "PB_2", "PB_3", "PB_4", "PB_5"),
  coding_transcript = c(14.75, 75.65, 35.02, 92.24, 75.73, 91.76, 93.4, 80.33, 76.65, 58.85, 28.43, 88.88, 88.6, 94.1, 85.59, 92.47, 25.53),
  Lib_Plat = c("cDNA-Illumina", "CapTrap-ONT", "dRNA-ONT", "R2C2-ONT", "cDNA-ONT", "dRNA-ONT", "R2C2-ONT", "cDNA-ONT", "cDNA-ONT", "cDNA-ONT", "cDNA-ONT", "dRNA-ONT", "CapTrap-PacBio", "cDNA-PacBio", "CapTrap-PacBio", "cDNA-PacBio", "cDNA-PacBio"),
  Data_Category = c("SO", "LO", "LS", "LO", "LO", "LO", "LO", "LO", "LS", "LO", "LS", "LO", "LO", "LO", "LO", "LO", "LS"),
  Tool = c("rnaSPAdes", "Bambu", "rnaSPAdes", "StringTie2\nIsoQuant", "Bambu", "Bambu", "Bambu", "RNA\nBloom", "RNA\nBloom", "StringTie2\nIsoQuant", "rnaSPAdes", "StringTie2\nIsoQuant", "Bambu", "Bambu", "StringTie2\nIsoQuant", "StringTie2\nIsoQuant", "rnaSPAdes")
)

manatee_coding <- data.frame(
  Label = c("illumina1", "ONT2", "ONT3", "ONT4", "ONT5", "PB1", "PB2", "PB3"),
  coding_transcript = c(4.7, 49.64, 66.9, 7.73, 44.53, 64.35, 68.8, 10.06),
  Lib_Plat = c("cDNA-Illumina", "cDNA-ONT", "cDNA-ONT", "cDNA-ONT", "cDNA-ONT", "cDNA-PacBio", "cDNA-PacBio", "cDNA-PacBio"),
  Data_Category = c("SO", "LO", "LO", "LS", "LS", "LO", "LO", "LS"),
  Tool = c("rnaSPAdes", "RNA\nBloom", "StringTie2\nIsoQuant", "rnaSPAdes", "RNA\nBloom", "Bambu", "StringTie2\nIsoQuant", "rnaSPAdes")
)

SX4.1 <- ggplot(ES_coding, aes(x=Label, y=as.numeric(coding_transcript), color=Lib_Plat, shape=Data_Category))  + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(coding_transcript), color=Lib_Plat), size=2) +
  geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
  pub_theme+
  theme_pubclean(flip=TRUE)+
  theme( axis.text.y  = element_blank(),
         axis.text.x = element_text(size=14),
         axis.ticks.y = element_blank()) +
  scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
  scale_color_manual(values =  libplat.palette, name="Library-Platform") +
  scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                            "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS"),
                   labels=c("SO", rep("LO", 6), rep("LS", 3)))+
  facet_grid(Tool ~., drop = TRUE, scales="free_y" ) +
  xlab("") + 
  theme(axis.title=element_text(size=16),
        strip.text.y = element_blank(), 
        axis.line.y = element_blank()) +
  theme(legend.position = "none") +
  scale_y_reverse("", sec.axis = sec_axis(~ . , breaks = NULL, name = "Mouse ES"),
                  label = unit_format(unit = "%"), limits=c(100 ,-3), expand = expansion(mult = c(0.1,0)))+
  coord_flip()

SX4.2 <- ggplot(manatee_coding, aes(x=Label, y=as.numeric(coding_transcript), color=Lib_Plat, shape=Data_Category)) + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(coding_transcript), color=Lib_Plat), size=2) +
  geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
  pub_theme+
  theme_pubclean(flip=TRUE)+
  theme( axis.text.y  = element_blank(),
         axis.text.x = element_text(size=14),
         axis.ticks.y = element_blank()) +
  scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
  scale_color_manual(values =  libplat.palette, name="Library-Platform") +
  scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                            "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS"),
                   labels=c("SO", rep("LO", 6), rep("LS", 3)))+
  facet_grid(Tool ~., drop = TRUE, scales="free_y" ) +
  xlab("") + 
  scale_y_continuous("", sec.axis = sec_axis(~ . , breaks = NULL, name = "Manatee"), 
                     label = unit_format(unit = "%"),  limits=c(-1, 100), expand = expansion(mult = c(0,0.1)))+
  theme(strip.text.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.title=element_text(size=16)) +
  theme(legend.position = "none") +
  coord_flip() 

SX.mid3<- ggplot(ES_coding,aes(x=1,y=Tool))+
  ggtitle("")+
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(1,1))+
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA))

gg4.1 <- ggplot_gtable(ggplot_build(SX4.1))
gg4.2 <- ggplot_gtable(ggplot_build(SX4.2))
g.mid <- ggplot_gtable(ggplot_build(SX.mid3))

Ch3S4 <- grid.arrange(gg4.1,g.mid, gg4.2,ncol=3,widths=c(4/9,1/9,4/9),
                    top = textGrob("Transcripts with coding potential (%)",gp=gpar(fontsize=18,font=1)))

ggsave(file=paste0(outdir, "/Ch3S4.svg"), plot=Ch3S4, width=8, height=5)

suppl = "66"
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Coding potential of transcripts detected by Challenge 3 submissions.")
pdf(paste0(outdir, "/Extended_Fig._66_b.pdf"))
annotate_figure(Ch3S2,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Fig. 67
###################

library("ggpubr")

cat.palette = c( "FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "GenicGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                 "Intergenic" = "darksalmon", "GenicIntron"="#41B6C4")

# Structural category

# this was original from, however the RData file was huge, so the relevant data was saved
#    load("SIRVs_manatee_paper.RData")
#    SCmanateeSIRVs <- df %>% group_by(sample) %>% dplyr::count(structural_category) %>% mutate(prop=n/sum(n)*100) 
#    the above df is saved as the file 'nb_SIRV_reads_manatee_by_SQANTI_category.csv' 

#    how_many_SIRV_transcripts_with_RM <- df %>% ungroup() %>% dplyr::filter(structural_category=="full-splice_match" & abs(diff_to_TSS)<50 & abs(diff_to_TTS)<50) %>% group_by(sample) %>% summarize(count_distinct = n_distinct(associated_transcript)) 
#    the above df is saved as the file 'how_many_SIRV_transcripts_with_RM.csv'

SCmanateeSIRVs <- read.table("Challenge3_Figures_Data/SIRVs/nb_SIRV_reads_manatee_by_SQANTI_category.csv",sep=" ",header=TRUE)

# frame uses different names for categories than other graphs, map them
catmap <- c("antisense" = "Antisense",  
            "full-splice_match" = "FSM",        
            "fusion" = "Fusion",     
            "genic" = "GenicGenomic",
            "incomplete-splice_match" = "ISM",        
            "intergenic" = "Intergenic", 
            "novel_in_catalog" = "NIC",        
            "novel_not_in_catalog" = "NNC",        
            "genic_intron" = "GenicIntron")

SCmanateeSIRVs$structural_category_label <- catmap[SCmanateeSIRVs$structural_category]

fig67a <- ggplot(SCmanateeSIRVs, aes(x=sample,y=prop,fill=structural_category_label,
                                     group = structural_category_label, colour = structural_category_label)) + 
  geom_bar(position="stack", stat="identity",color="black",width = 0.5) + 
  theme_bw() + xlab("Manatee samples") +
  theme(aspect.ratio=0.8) +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.text = element_text(size= 8),legend.key.size = unit(0.5, 'cm'), legend.title=element_blank()) + 
  scale_fill_manual(name = "Structural_category", values = cat.palette , 
                         limits = names(cat.palette)) + 
  ylab("% SQANTI categories in SIRV reads") + theme(text=element_text(size=10))
  

# Plot for # SIRV transcripts detected with reference-match read (out of 69)
how_many_SIRV_transcripts_with_RM <- read.table("Challenge3_Figures_Data/SIRVs/how_many_SIRV_transcripts_with_RM.csv",header=TRUE)

fig67b <- ggplot(how_many_SIRV_transcripts_with_RM, aes(x=sample, y=count_distinct,fill=sample)) +
  geom_bar(stat="identity", position="dodge") + 
  theme_bw() + theme(aspect.ratio=0.8) +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),  legend.text = element_text(size=8), legend.key.size = unit(0.5, 'cm'), legend.title=element_blank()) + 
  ylab("# SIRV transcripts detected with RM reads") + xlab("Manatee sample") + 
  theme(text=element_text(size=10)) + 
  scale_fill_manual(values=c("#9DACBB","#9DACBB","#9DACBB","#104E8B","#104E8B","#104E8B")) + 
  geom_hline(yintercept=69, linetype="dashed",color = "darkred", linewidth=1) + 
  theme(legend.position = "none") + scale_y_continuous(breaks=seq(0,70,5)) + 
  geom_text(data = how_many_SIRV_transcripts_with_RM, aes(label = count_distinct), vjust=-0.5, hjust=0.5)

fig67 <- ggarrange(fig67a,fig67b,
                      labels = c( "a)", "b)"),
                      ncol = 2, nrow = 1, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 

suppl = "67"
mylegend <- paste0("     Extended Fig. ", suppl, ". SQANTI3 analysis of SIRV reads in manatee samples. a) SQANTI3 categories for reads mapping to SIRVs in cDNA-PacBio and cDNA-ONT replicates. \n     b) Number of SIRV transcripts with at least one Reference Match (RM) read in cDNA-PacBio and cDNA-ONT replicates")
pdf(paste0(outdir, "/Extended_Fig._67.pdf"), width = 10, height = 6)
annotate_figure(fig67,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()
