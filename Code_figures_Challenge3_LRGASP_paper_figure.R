#########################################
######### FIGURE 4: CHALLENGE 3 #########
#########################################
## author: Francisco J. Pardo-Palacios, f.pardo.palacios@gmail.com
## author:Ana Conesa, ana.conesa@csic.es
## Last modified: March 16th 2023
#######################################

# Install and load packages
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(scales)
library(patchwork)
library(gridExtra)
library(grid)
library(RColorConesa)
library(fmsb)
library(MetBrewer)

#### set theme for plots
pub_theme <- theme_pubclean(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y  = element_text(vjust=0.5, size=14) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=10), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=14)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  theme(legend.position = "none")

old.libplat.palette = c( "cDNA+ONT"="#74CDF0", "cDNA+PacBio"="#EE446F", "cDNA+Illumina"="#FFCF71", 
                         "CapTrap+ONT"="#7482F0", "R2C2+ONT"="#74F0D9", "dRNA+ONT"="#13BF5E", "CapTrap+PacBio"="#d14141")

libplat.palette = c( "cDNA-PacBio"="#c06636", "CapTrap-PacBio"="#802417", "cDNA-Illumina"="#e8b960",  "Freestyle-Freestyle"="#ce9344",
                     "cDNA-ONT"="#646e3b", "CapTrap-ONT"="#17486f", "R2C2-ONT"="#508ea2", "dRNA-ONT"="#2b5851"
)
cat.palette = c( "FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "GenicGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                 "Intergenic" = "darksalmon", "GenicIntron"="#41B6C4")


### Panel 4a
#############

# Process the code file
ES_code <- read.csv("ES/ES_code.txt", sep=",", header = T )
ES_code$Lib_Plat <- apply(ES_code, 1, function(x){
  paste(x["Library_Preps"], x["Platform"], sep = "-")
})
ES_code$Lib_DC=apply(cbind(ES_code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
ES_code$Label <-apply(cbind(ES_code[,c("Library_Preps", "Platform", "Data_Category")]), 1, paste, collapse="-")

ES_summary_table <- read.csv("ES_challenge1/ES_challenge1_metrics.summary_table_SC.csv", header = T)
ES_summary_table <- merge(ES_summary_table, ES_code, by.x="ID", by.y="pipelineCode")
ES_summary_table$Tool <- gsub("-", "\n", ES_summary_table$Tool)

melted_ES_summary <- ES_summary_table %>%
  pivot_longer(c("FSM","ISM","NIC","NNC","Antisense","Fusion","GenicGenomic","GenicIntron","Intergenic"))
melted_ES_summary$Tool <- gsub("-", "\n", melted_ES_summary$Tool)

melted_ES_summary$name <- melted_ES_summary$name %>%  factor(levels=c("FSM", "ISM", "NIC", "NNC",
                                                                      "Antisense", "Fusion", "GenicGenomic", "GenicIntron",
                                                                      "Intergenic"),
                                                             labels=c("FSM", "ISM", "NIC", "NNC",
                                                                      "Antisense", "Fusion", "GenicGenomic", "GenicIntron",
                                                                      "Intergenic"))

p.A1 <- ggplot(ES_summary_table, aes(x=Label, y=total, fill=Lib_Plat))+
  geom_bar(stat="identity", width=0.7)+
  facet_grid( .~Tool , scales="free_x", space="free_x") +
  scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
  pub_theme +
  ylab("Number of \ntranscripts") +
  xlab("")+
  theme(  axis.text.x = element_blank(),
         strip.text.x = element_text(size=18)) +
  scale_y_continuous(label = unit_format(unit = "K", scale = 0.001), expand = expansion(mult = c(0,0.1)))+
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  theme(plot.margin=unit(c(0,0,0,0), "cm"))


p.A2.Labels= c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
               "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS")
p.A2 <- ggplot(melted_ES_summary, aes(x=Label, y=value, fill=name)) + 
  geom_bar(position = "fill", stat="identity", width = 0.7) +
  facet_grid( .~Tool , scales="free_x", space="free_x") +
  scale_fill_manual(values =  cat.palette, name="" ) +
  scale_x_discrete(breaks=p.A2.Labels,
                   labels=c("SO", rep("LO", 6), rep("LS", 3)))+
  xlab("")+ 
  ylab("Structural Category \nDistribution")+
  pub_theme +
  theme(text = element_text(size = 9)) +
  theme( strip.text.x = element_blank(),
         legend.text = element_text(size=14)) +
  scale_y_continuous(label = unit_format(unit = "%", scale = 100), expand = expansion(mult = c(0,0.1)))+
  #labs(tag = "LO:Long Only\nLS:Long and Short", hjust = 0) +
  theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1 )) +
  theme(plot.margin=unit(c(-1,0,0,0), "cm"), plot.tag.position = c(0.95,0.95), ) 

pA <- p.A1 / p.A2 + 
  plot_layout(heights = c(1, 1), ncol = 1)

ggsave(file="panel4a.svg", plot=pA, width=12, height=6)


##### version to show mono-exons with alpha
df_SC_monoexons <- read.csv("ES_challenge1/ES_challenge1_metrics.monoexons_SC.csv", sep=",", header = T) %>% t() %>% as.data.frame()
df_SC_monoexons$total_mono <- apply(df_SC_monoexons,1,sum)
df_SC_monoexons$ID <- row.names(df_SC_monoexons)
total_mono <- df_SC_monoexons$total_mono

ES_summary_table <- merge(ES_summary_table, df_SC_monoexons[, c("ID", "total_mono")], by="ID")
ES_summary_table$total_multi <- apply(ES_summary_table,1, function(x){
  as.numeric(x["total"]) - as.numeric(x["total_mono"])
})

melted_ES_summary_exons <- ES_summary_table %>%
  pivot_longer(c("total_multi", "total_mono"))
melted_ES_summary_exons$Tool <- gsub("-", "\n", melted_ES_summary_exons$Tool)

p.A1_exons <- ggplot(melted_ES_summary_exons, aes(x=Label, y=value, alpha=name, fill=Lib_Plat))+
  geom_bar(stat="identity", width=0.7)+
  facet_grid( .~Tool , scales="free_x", space="free_x") +
  scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
  scale_alpha_manual(values=c(0.6,1), breaks=c("total_mono", "total_multi"), labels=c("Mono-exons", "Multi-exons"))+
  pub_theme +
  ylab("Num. total isoforms") +
  xlab("")+
  theme( axis.text.x = element_blank(),
         strip.text.x = element_text(size=18)) +
  scale_y_continuous(label = unit_format(unit = "K", scale = 0.001), expand = expansion(mult = c(0,0.1)))+
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

pA_alt <- p.A1_exons / p.A2 + 
  plot_layout(heights = c(1, 1))

ggsave(file="panel4a_alt.svg", plot=pA_alt, width=12, height=6)

### End Panel 4a
################

### Panel 4b
############

# Process the code file
manatee_code <- read.csv("manatee/manatee_code.txt", sep=",", header = T )
manatee_code$Lib_Plat <- apply(manatee_code, 1, function(x){
  paste(x["Library_Preps"], x["Platform"], sep = "-")
})
manatee_code$Lib_DC=apply(cbind(manatee_code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
manatee_code$Label <-apply(cbind(manatee_code[,c("Library_Preps", "Platform", "Data_Category")]), 1, paste, collapse="-")

# Obtain metrics for manatee
manatee_metrics <- read.csv("manatee/manatee_challenge3_metrics.challenge3_metrics.csv", sep=",",header = T) %>% t()
mono_exons_df<-data.frame(Row.names=c("ONT1","PB1","ONT2", "illumina1","ONT3","ONT4","PB2","PB3","ONT5"), 
                          Num.monoexons=c(0,7,79174,635028,3703,470801,52,246565,85291))
rownames(mono_exons_df) <- mono_exons_df$Row.names
manatee_metrics <- merge(manatee_metrics, mono_exons_df,by.x=0, by.y="Row.names")
manatee_metrics$Num.multiexons <- apply(manatee_metrics,1, function(x){
  as.numeric(x["Number of transcripts"]) - as.numeric(x["Num.monoexons"])
})
manatee_metrics <- merge(manatee_metrics, manatee_code, by.x="Row.names", by.y="pipelineCode")
manatee_metrics$Tool <- gsub("-", "\n", manatee_metrics$Tool)
colnames(manatee_metrics) <- make.names(colnames(manatee_metrics))

melted_manatee <- manatee_metrics %>%
  pivot_longer(c("Num.monoexons", "Num.multiexons"))
melted_manatee$Tool <- gsub("-", "\n", melted_manatee$Tool)

# Process info number isoforms per locus
trx_locus.list <- list()
for (i in c("ONT1", "ONT2","ONT3","ONT4","ONT5","PB1","PB2","PB3","illumina1")){
  trx_file=paste0("manatee/manatee_isoforms_per_gene/",i,"_cpm_vs_trans.tsv")
  df=read.csv(trx_file,sep="\t",header = T)
  num_trx_locus <- data.frame(Cat=c("1","2-3","4-5", ">6"), 
                              Num_locus=c(which(df$n_isoforms==1) %>% length(), 
                                          which(df$n_isoforms==2 |  df$n_isoforms==3) %>% length(),
                                          which(df$n_isoforms==4 |  df$n_isoforms==5) %>% length(),
                                          which(df$n_isoforms>5) %>% length()
                              )
  )
  trx_locus.list[[i]] <- num_trx_locus
}

trx_locus.df <- bind_rows(trx_locus.list, .id = "pipeline")
trx_locus.df <- merge(trx_locus.df, manatee_code, by.x="pipeline", by.y="pipelineCode")
trx_locus.df$Tool <- gsub("-", "\n", trx_locus.df$Tool)

melted_manatee <- melted_manatee %>% filter(Row.names!="ONT1")
pF4.1 <- ggplot(melted_manatee, aes(x=Row.names, y=value, alpha=name, fill=Lib_Plat))+
  geom_bar(stat="identity", position = "stack", width=0.7)+
  facet_grid( .~Tool , scales="free_x", space="free_x") +
  scale_fill_manual(values =  libplat.palette, name="Library-Platoform") +
  scale_alpha_manual(values=c(0.6,1), breaks=c("Num.monoexons", "Num.multiexons"), labels=c("Mono-exons", "Multi-exons"))+
  scale_x_discrete(breaks=c("PB1","ONT3","PB2","ONT2","ONT5","illumina1","ONT4", "PB3"),
                   labels=c("n=1.9K","n=63K","n=25K","n=179K","n=177K","n=916K","n=543K","n=294K")) +
  pub_theme +
  ylab("Num. total isoforms") +
  xlab("")+
  theme(strip.text.x = element_text(size=18)) +
  scale_y_continuous(label = unit_format(unit = "K", scale = 0.001), expand = expansion(mult = c(0,0.1)))+
  theme(legend.position = "none") +
  theme(plot.margin=unit(c(0,0,-0.5,0), "cm"))

trx_locus.df <- trx_locus.df %>% filter(pipeline!="ONT1")
pF4.2 <- ggplot(trx_locus.df, aes(x=pipeline, y=Num_locus, fill=Cat)) + 
  geom_bar(position = "fill", stat="identity", width = 0.7) +
  facet_grid( .~Tool , scales="free_x", space="free_x", drop=TRUE) +
  scale_fill_met_d(name = "Cassatt1", direction = -1) +
  scale_x_discrete(breaks=c("PB1","ONT3","PB2","ONT2","ONT5","illumina1","ONT4", "PB3"),
                   labels=c("LO", "LO", "LO", "LO", "LS","SO", "LS","LS")) +
  ylab("Number of transcripts \n per gene")+
  pub_theme +
  theme( strip.text.x = element_blank(),
         axis.title.x = element_blank(),
         legend.text = element_text(size=14),
         legend.title = element_blank()) +
  scale_y_continuous(label = unit_format(unit = "%", scale = 100), expand = expansion(mult = c(0,0.1)))+
  theme(legend.position = "bottom") +
  theme(plot.margin=unit(c(-1,0,0,0), "cm"))

pF <- pF4.1 / pF4.2 +
  plot_layout(heights = c(1, 1), ncol = 1)
ggsave(file="panel4b.svg", plot=pF, width=12, height=6)


### End Panel 4b
################

### Panel 4c
############

######### length distribution mouse
dist.list_ES <- list()
for (i in c("ONT_1", "ONT_2","ONT_3","ONT_4","ONT_5",
            "ONT_6", "ONT_7","ONT_8","ONT_9","ONT_10",
            "ONT_11","PB_1","PB_2","PB_3",
            "PB_4", "PB_5", "illumina_1")){
  length_file=paste0("ES/length_distributions/",i,"_length_dist.tsv")
  df=read.csv(length_file,sep="\t",header = T)
  dist.list_ES[[i]] <- as.numeric(df$length) %>% as.data.frame()
}

dist_df_ES <- bind_rows(dist.list_ES, .id = "pipeline")
colnames(dist_df_ES)<-c("pipeline","length")
dist_df_ES <- merge(dist_df_ES, ES_code, by.x="pipeline", by.y="pipelineCode")
dist_df_ES$Tool <- gsub("-", "\n", dist_df_ES$Tool)

pF3 <- ggplot(dist_df_ES, aes(x=Label, y=length, fill=Lib_Plat))+
  geom_violin()+
  geom_boxplot(width=0.15, color="white", alpha=0.8, outlier.shape = NA) +
  facet_grid( .~Tool , scales="free_x", space="free_x") +
  scale_fill_manual(values =  libplat.palette, name="Library-Platoform") +
  scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                            "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS"),
                   labels=c("SO", rep("LO", 6), rep("LS", 3)))+
  xlab("")+ 
  ylab("Length (bp), log10")+
  pub_theme +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(size=16)) +
  scale_y_continuous(trans='log10', breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)) )
ggsave(file="panel4c.svg", plot=pF3, width=12, height=6)

### End Panel 4c
################

### Panel 4d
############
dist.list_manatee <- list()
for (i in c("ONT1", "ONT2","ONT3","ONT4","ONT5","PB1","PB2","PB3","illumina1")){
  length_file=paste0("manatee/length_distributions/",i,"_length_dist.tsv")
  df=read.csv(length_file,sep="\t",header = T)
  dist.list_manatee[[i]] <- as.numeric(df$length) %>% as.data.frame()
}

dist_df_manatee <- bind_rows(dist.list_manatee, .id = "pipeline")
colnames(dist_df_manatee)<-c("pipeline","length")
dist_df_manatee <- merge(dist_df_manatee, manatee_code, by.x="pipeline", by.y="pipelineCode")
dist_df_manatee$Tool <- gsub("-", "\n", dist_df_manatee$Tool)

sample_size = dist_df_manatee %>% group_by(pipeline) %>% summarize(num=n())

pF2 <- ggplot(dist_df_manatee, aes(x=Label, y=length, fill=Lib_Plat))+
  geom_violin()+
  geom_boxplot(width=0.15, color="white", alpha=0.8, outlier.shape = NA) +
  facet_grid( .~Tool , scales="free_x", space="free_x") +
  scale_fill_manual(values =  libplat.palette, name="Library-Platoform") +
  scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                            "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS"),
                   labels=c("SO", rep("LO", 6), rep("LS", 3)))+
  xlab("")+ 
  ylab("Length (bp), log10")+
  pub_theme +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(size=16)) +
  scale_y_continuous(trans='log10', breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)) )

ggsave(file="panel4d.svg", plot=pF2, width=12, height=6)

### Panel 4e1
#############

ES_metrics_perc <- read.csv("ES/ES_challenge3_metrics.challenge3_metrics_perc.csv", sep=",",header = T) %>% t()
ES_metrics_perc <- merge(ES_metrics_perc, ES_code, by.x=0, by.y="pipelineCode")
ES_metrics_perc$Tool <- gsub("-", "\n", ES_metrics_perc$Tool)
colnames(ES_metrics_perc) <- make.names(colnames(ES_metrics_perc))
ES_metrics_perc[12, "Mapping.transcripts"] <- 100

pC_SJ <- ggplot(ES_metrics_perc, aes(x=Label, y=as.numeric(Splice.Junctions.with.short.read.coverage), color=Lib_Plat, shape=Data_Category)) + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(Splice.Junctions.with.short.read.coverage), color=Lib_Plat), size=2) +
  geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
  pub_theme+
  theme(axis.title.y = element_text(size=18),
        axis.text.y  = element_text(size=18) ) +
  scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
  scale_color_manual(values =  libplat.palette, name="Library-Platform") +
  scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                            "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS"),
                  labels=c("SO", rep("LO", 6), rep("LS", 3)))+
  facet_grid(.~ Tool, drop = TRUE, scales="free_x", space = "free_x") +
  xlab("") + ylab("% Splice Junctions without short read coverage")+
  scale_y_continuous(label = label_percent(scale = 1), limits=c(0,100), expand = expansion(mult = c(0,0.1)))+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 18))

ggsave(file="panel4e1.svg", plot=pC_SJ, width=8, height=8)

### End Panel 4e1
################

### Panel 4e2
#############

CAGE_QUANT_csv <- read.csv("ES_challenge1/ES_challenge1_metrics.CAGE_QuantSeq.csv", header = T)
CAGE_QUANT <- merge(CAGE_QUANT_csv, ES_code, by.x=0, by.y="pipelineCode")
CAGE_QUANT$total <- apply(CAGE_QUANT, 1,function(x){
  as.numeric(x["CAGE"])+as.numeric(x["noCAGE"])
})
CAGE_QUANT$CAGE_perc <- apply(CAGE_QUANT, 1,function(x){
  as.numeric(x["CAGE"])/as.numeric(x["total"]) *100
})
CAGE_QUANT$Quant_perc <- apply(CAGE_QUANT, 1,function(x){
  as.numeric(x["QuantSeq"])/as.numeric(x["total"]) *100
})

CAGE_QUANT$Tool <- gsub("-", "\n", CAGE_QUANT$Tool)

pCAGE <- ggplot(CAGE_QUANT, aes(x=Label, y=as.numeric(CAGE_perc), color=Lib_Plat, shape=Data_Category)) + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(CAGE_perc), color=Lib_Plat), size=2) +
  geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
  pub_theme+
  theme(axis.title.y = element_text(size=18),
        axis.text.y  = element_text(size=18) ) +
  scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
  scale_color_manual(values =  libplat.palette, name="Library-Platform") +
  scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                            "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS"),
                   labels=c("SO", rep("LO", 6), rep("LS", 3)))+
  facet_grid(.~ Tool, drop = TRUE, scales="free_x", space = "free_x") +
  xlab("") + ylab("% Transcripts supported by CAGE")+
  scale_y_continuous(label = label_percent(scale = 1), limits=c(0,100), expand = expansion(mult = c(0,0.1)))+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 18))

ggsave(file="panel4e2.svg", plot=pCAGE, width=8, height=8)


### End Panel 4e2
################

### Panel 4e3
#############

pC_nonCan <- ggplot(ES_metrics_perc, aes(x=Label, y=as.numeric(Non.canonical.Splice.Junctions), color=Lib_Plat, shape=Data_Category))  + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(Non.canonical.Splice.Junctions), color=Lib_Plat), size=2) +
  geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
  pub_theme+
  theme(axis.title.y = element_text(size=18),
        axis.text.y  = element_text(size=18) ) +
  scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
  scale_color_manual(values =  libplat.palette, name="Library-Platform") +
  scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                            "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS"),
                   labels=c("SO", rep("LO", 6), rep("LS", 3)))+
  facet_grid(.~ Tool, drop = TRUE, scales="free_x", space = "free_x") +
  xlab("") + ylab("% Non-canonical Splice Junctions")+
  scale_y_continuous(label = label_percent(scale = 1), limits=c(0,100), expand = expansion(mult = c(0,0.1)))+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 18))

ggsave(file="panel4e3.svg", plot=pC_nonCan, width=8, height=8)

### End Panel 4e3
#################

### Panel 4e4
#############

pQuant <- ggplot(CAGE_QUANT, aes(x=Label, y=as.numeric(Quant_perc), color=Lib_Plat, shape=Data_Category)) + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(Quant_perc), color=Lib_Plat), size=2) +
  geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
  pub_theme+
  theme(axis.title.y = element_text(size=18),
        axis.text.y  = element_text(size=18) ) +
  scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
  scale_color_manual(values =  libplat.palette, name="Library-Platform") +
  scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                            "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS"),
                   labels=c("SO", rep("LO", 6), rep("LS", 3)))+
  facet_grid(.~ Tool, drop = TRUE, scales="free_x", space = "free_x") +
  xlab("") + ylab("% Transcripts supported by Quant-seq")+
  scale_y_continuous(label = label_percent(scale = 1), limits=c(0,100), expand = expansion(mult = c(0,0.1)))+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 18))

ggsave(file="panel4e4.svg", plot=pQuant , width=8, height=8)

### End Panel 4e4
#################


### Panel 4f1
#############

manatee_BUSCO <- read.csv("manatee/manatee_challenge3_metrics.BUSCO_metrics.csv", sep=",",header = T) %>% t()
manatee_BUSCO <- merge(manatee_BUSCO, manatee_code, by.x=0, by.y="pipelineCode")
manatee_BUSCO$Tool <- gsub("-", "\n", manatee_BUSCO$Tool)
colnames(manatee_BUSCO) <- make.names(colnames(manatee_BUSCO))
manatee_BUSCO <- manatee_BUSCO %>% filter(Row.names!="ONT1")

ES_BUSCO <- read.csv("ES/ES_challenge3_metrics.BUSCO_metrics.csv", sep=",",header = T) %>% t()
ES_BUSCO <- merge(ES_BUSCO, ES_code, by.x=0, by.y="pipelineCode")
ES_BUSCO$Tool <- gsub("-", "\n", ES_BUSCO$Tool)
colnames(ES_BUSCO) <- make.names(colnames(ES_BUSCO))

total_busco <- 11366
get_perc_busco_found <- function(x){
  r=(total_busco - as.numeric(x["Missing.BUSCOs"]))*100/total_busco 
  round(r, digits = 2)
}

get_perc_busco_complete <- function(x){
  c <- as.numeric(x["Complete.and.single.copy.BUSCOs"]) + as.numeric(x["Complete.and.duplicated.BUSCOs"]) 
  r=c*100/total_busco 
  round(r, digits = 2)
}

get_perc_busco_fragmented <- function(x){
  c <- as.numeric(x["Fragmented.BUSCOs"]) 
  r=c*100/total_busco 
  round(r, digits = 2)
}

get_perc_busco_missing <- function(x){
  c <- as.numeric(x["Missing.BUSCOs"]) 
  r=c*100/total_busco 
  round(r, digits = 2)
}

ES_BUSCO$BUSCO_found <- apply(ES_BUSCO,1,get_perc_busco_found)
manatee_BUSCO$BUSCO_found <- apply(manatee_BUSCO,1,get_perc_busco_found)

ES_BUSCO$BUSCO_complete <- apply(ES_BUSCO,1,get_perc_busco_complete)
manatee_BUSCO$BUSCO_complete <- apply(manatee_BUSCO,1,get_perc_busco_complete)

ES_BUSCO$BUSCO_fragmented <- apply(ES_BUSCO,1,get_perc_busco_fragmented)
manatee_BUSCO$BUSCO_fragmented <- apply(manatee_BUSCO,1,get_perc_busco_fragmented)

ES_BUSCO$BUSCO_missing <- apply(ES_BUSCO,1,get_perc_busco_missing)
manatee_BUSCO$BUSCO_missing <- apply(manatee_BUSCO,1,get_perc_busco_missing)

pD7.1<- ggplot(ES_BUSCO, aes(x=Label, y=as.numeric(BUSCO_found), color=Lib_Plat, shape=Data_Category))  + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(BUSCO_found), color=Lib_Plat), size=2) +
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
                  label = unit_format(unit = "%"), limits=c(80,-2), expand = expansion(mult = c(0.1,0)))+
  coord_flip()

pD7.2 <- ggplot(manatee_BUSCO, aes(x=Label, y=as.numeric(BUSCO_found), color=Lib_Plat, shape=Data_Category)) + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(BUSCO_found), color=Lib_Plat), size=2) +
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
                     label = unit_format(unit = "%"),  limits=c(-1, 80), expand = expansion(mult = c(0,0.1)))+
  theme(strip.text.y = element_blank(), 
       axis.line.y = element_blank(),
       axis.title=element_text(size=16)) +
  theme(legend.position = "none") +
  coord_flip() 

p.mid7<- ggplot(ES_BUSCO,aes(x=1,y=Tool))+
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

gg1.7 <- ggplot_gtable(ggplot_build(pD7.1))
gg2.7 <- ggplot_gtable(ggplot_build(pD7.2))
g.mid <- ggplot_gtable(ggplot_build(p.mid7))

pD7 <- grid.arrange(gg1.7,g.mid, gg2.7,ncol=3,widths=c(4/9,1/9,4/9),
                    top = textGrob("BUSCO genes found (%)",gp=gpar(fontsize=18,font=1)))

ggsave(file="panel4f1.svg", plot=pD7, width=8, height=5)

### End Panel 4f1
################

### Panel 4f2
#############

pD7.1f <- ggplot(ES_BUSCO, aes(x=Label, y=as.numeric(BUSCO_fragmented), color=Lib_Plat, shape=Data_Category))  + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(BUSCO_fragmented), color=Lib_Plat), size=2) +
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
                  label = unit_format(unit = "%"), limits=c(20,-2), expand = expansion(mult = c(0.1,0)))+
  coord_flip()

pD7.2f <- ggplot(manatee_BUSCO, aes(x=Label, y=as.numeric(BUSCO_fragmented), color=Lib_Plat, shape=Data_Category)) + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(BUSCO_fragmented), color=Lib_Plat), size=2) +
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
                     label = unit_format(unit = "%"),  limits=c(-1, 20), expand = expansion(mult = c(0,0.1)))+
  theme(axis.title=element_text(size=16),
        strip.text.y = element_blank(), 
        axis.line.y = element_blank()) +
  theme(legend.position = "none") +
  coord_flip() 

gg1.7f <- ggplot_gtable(ggplot_build(pD7.1f))
gg2.7f <- ggplot_gtable(ggplot_build(pD7.2f))

pD7f <- grid.arrange(gg1.7f,g.mid, gg2.7f,ncol=3,widths=c(4/9,1/9,4/9),
                     top = textGrob("BUSCO genes fragmented (%)",gp=gpar(fontsize=18,font=1)))

ggsave(file="panel4f2.svg", plot=pD7f, width=8, height=5)

### End Panel 4f2
#################

### Panel 4f3
#############

pD7.1m <- ggplot(ES_BUSCO, aes(x=Label, y=as.numeric(BUSCO_missing), color=Lib_Plat, shape=Data_Category))  + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(BUSCO_missing), color=Lib_Plat), size=2) +
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
                  label = unit_format(unit = "%"), limits=c(100,-3), expand = expansion(mult = c(0.1,0)))+
  coord_flip()

pD7.2m <- ggplot(manatee_BUSCO, aes(x=Label, y=as.numeric(BUSCO_missing), color=Lib_Plat, shape=Data_Category)) + 
  geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(BUSCO_missing), color=Lib_Plat), size=2) +
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
                     label = unit_format(unit = "%"),  limits=c(-2, 100), expand = expansion(mult = c(0,0.1)))+
  theme(axis.title=element_text(size=16),
        strip.text.y = element_blank(), 
        axis.line.y = element_blank()) +
  theme(legend.position = "none") +
  coord_flip() 

gg1.7m<- ggplot_gtable(ggplot_build(pD7.1m))
gg2.7m <- ggplot_gtable(ggplot_build(pD7.2m))

pD7m <- grid.arrange(gg1.7m,g.mid, gg2.7m,ncol=3,widths=c(4/9,1/9,4/9),
                     top = textGrob("BUSCO genes missing (%)",gp=gpar(fontsize=20,font=1)))

ggsave(file="panel4f3.svg", plot=pD7m, width=8, height=5)


### End Panel 4f3
#################

#### Panel 4g
#############

p.leg <- ggplot(melted_ES_summary_exons, aes(x=Label, y=value, alpha=name, fill=Lib_Plat))+
  geom_bar(stat="identity", width=0.7)+
  facet_grid( .~Tool , scales="free_x", space="free_x") +
  scale_fill_manual(values = c(libplat.palette,c("LO","LS","SO")) , name="Library-Platform") +
  scale_alpha_manual(values = c(0.6,1), breaks=c("total_mono", "total_multi"), labels=c("Mono-exon", "Multi-exon"),name="") +
  pub_theme +
  ylab("Num. total isoforms") +
  xlab("")+
  theme( axis.text.y  = element_text( size=18),
         axis.title.y  = element_text( size=18),
         axis.text.x = element_blank(),
         strip.text.x = element_text(size=18)) +
  scale_y_continuous(label = unit_format(unit = "K", scale = 0.001), expand = expansion(mult = c(0,0.1)))+
  theme(legend.position = "right",
        axis.text.x = element_blank()) +
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

leg <- get_legend(p.leg)
p_leg <- as_ggplot(leg)

ggsave(file="Panel4g.svg", plot=p_leg, width=5, height=8)

### End Panel 4g
################

### Panel 4h1
#############

ES_SIRV_metrics <- read.csv("ES_challenge1/ES_challenge1_metrics.SIRVS_metrics.csv", sep=",",header = T) %>% t()
ES_SIRV_metrics <- merge(ES_SIRV_metrics, ES_code, by.x=0, by.y="pipelineCode")
ES_SIRV_metrics$"1/Red" <- round(1/ES_SIRV_metrics$Redundancy,2)
ES_SIRV_metrics$Label <-apply(cbind(ES_SIRV_metrics[,c("Library_Preps", "Platform","Data_Category")]), 1, paste, collapse="_")

ES_SIRV_metrics2 <- data.frame (ES_SIRV_metrics$Sensitivity,
                                ES_SIRV_metrics$`Positive Detection Rate`,
                                ES_SIRV_metrics$Precision,
                                ES_SIRV_metrics$`Non Redundant Precision`,
                                ES_SIRV_metrics$`False Discovery Rate`,
                             ES_SIRV_metrics$`1/Red`,
                             ES_SIRV_metrics$Tool,
                             ES_SIRV_metrics$Label
)

colnames(ES_SIRV_metrics2) <- c("Sen", "PDR", "Pre", "nrPre", "FDR", "1/Red", "Tool", "Label")
mycolors = c( "cDNA_PacBio_LO"="#d8527c", "cDNA_PacBio_LS"="#9a133d", 
              "cDNA_ONT_LO"="#6996e3", "cDNA_ONT_LS"="#1a318b" )

tool = sort(unique(ES_SIRV_metrics2$Tool))
plat_lib <- unique(ES_SIRV_metrics2$Label)

pdf("panel4h1.pdf")
par(mar=c(0,0.5,1.5,0) + 0.1)
#layout(matrix(c(1:12,13,13,13,13),  byrow = TRUE, ncol = 4,nrow = 4)) # for main figure
layout(matrix(c(1:4,5,5),  byrow = TRUE, ncol = 2,nrow = 3))

for ( i in 1:4) {
  sel <- ES_SIRV_metrics2$Tool == tool[i]
  subms <- ES_SIRV_metrics2[sel,c(1:6)]
  subms <- rbind(rep(1,ncol(subms)) , rep(0,ncol(subms)), subms)
  radarchart(subms, axistype=1 , title= tool[i], 
             #custom polygon
             pcol= mycolors[ES_SIRV_metrics2[sel,"Label"]], plwd=2 , plty=1.3,
             #custom the grid
             cglcol = "grey", cglty=1, axislabcol="grey", cglwd = 0.8, caxislabels=seq(0,0.25,1), 
             #custom labels
             #vlabels = rep("",6), palcex = 0, cex.main=1.5, paxislabels = 0
             vlabels = colnames(subms) , palcex = 1, cex.main=1.2, paxislabels = 0, col = c(2,2,2,3,3,1)
  )
}

# Add a legend
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("blue","black", "green", "orange", "pink")
l <- legend(x = "top",inset = 0, legend = names(mycolors), horiz = TRUE, bty = "n", pch=20 , col=mycolors , text.col = "black", cex=1.2, pt.cex=3)
dev.off()

### End Panel 4h1
#################

### Panel 4h2
#############
manatee_SIRV_metrics <- read.csv("manatee/manatee_challenge3_metrics.SIRVS_metrics.csv", sep=",",header = T) %>% t()
manatee_SIRV_metrics <- merge(manatee_SIRV_metrics, manatee_code, by.x=0, by.y="pipelineCode")
manatee_SIRV_metrics$"1/Red" <- round(1/manatee_SIRV_metrics$Redundancy,2)
manatee_SIRV_metrics$Label <-apply(cbind(manatee_SIRV_metrics[,c("Library_Preps", "Platform","Data_Category")]), 1, paste, collapse="_")

mana_metrics2 <- data.frame (manatee_SIRV_metrics$Sensitivity,
                             manatee_SIRV_metrics$`Positive Detection Rate`,
                             manatee_SIRV_metrics$Precision,
                             manatee_SIRV_metrics$`Non Redundant Precision`,
                             manatee_SIRV_metrics$`False Discovery Rate`,
                             manatee_SIRV_metrics$`1/Red`,
                             manatee_SIRV_metrics$Tool,
                             manatee_SIRV_metrics$Label
)

colnames(mana_metrics2) <- c("Sen", "PDR", "Pre", "nrPre", "FDR", "1/Red", "Tool", "Label")
mycolors = c( "cDNA_PacBio_LO"="#d8527c", "cDNA_PacBio_LS"="#9a133d", 
              "cDNA_ONT_LO"="#6996e3", "cDNA_ONT_LS"="#1a318b" )

tool = sort(unique(mana_metrics2$Tool))
plat_lib <- unique(mana_metrics2$Label)

pdf("panel4h2.pdf")
par(mar=c(0,0.5,1.5,0) + 0.1)
#layout(matrix(c(1:12,13,13,13,13),  byrow = TRUE, ncol = 4,nrow = 4)) # for main figure
layout(matrix(c(1:4,5,5),  byrow = TRUE, ncol = 2,nrow = 3))

for ( i in 1:4) {
  sel <- mana_metrics2$Tool == tool[i]
  subms <- mana_metrics2[sel,c(1:6)]
  subms <- rbind(rep(1,ncol(subms)) , rep(0,ncol(subms)), subms)
  radarchart(subms, axistype=1 , title= tool[i], 
             #custom polygon
             pcol= mycolors[mana_metrics2[sel,"Label"]], plwd=2 , plty=1.3,
             #custom the grid
             cglcol = "grey", cglty=1, axislabcol="grey", cglwd = 0.8, caxislabels=seq(0,0.25,1), 
             #custom labels
             #vlabels = rep("",6), palcex = 0, cex.main=1.5, paxislabels = 0
             vlabels = colnames(subms) , palcex = 1, cex.main=1.2, paxislabels = 0, col = c(2,2,2,3,3,1)
  )
}

# Add a legend
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("blue","black", "green", "orange", "pink")
l <- legend(x = "top",inset = 0, legend = names(mycolors), horiz = TRUE, bty = "n", pch=20 , col=mycolors , text.col = "black", cex=1.2, pt.cex=3)
dev.off()

### End Panel 4h2
#################

