# This R script performs comprehensive analysis of phosphorylation sites in Alzheimer's disease (AD) pathology, integrating phosphoproteomics data with AlphaFold2-predicted structural and charge distribution features. The code systematically investigates how structural context influences phosphorylation changes across different AD stages (HPC, MCI, AD) and specifically analyzes tau protein phosphorylation sites. Key analytical components include: processing AD stage-specific phosphoproteomics data to calculate log2 fold changes relative to control samples; integrating structural topology metrics (CF10QS) with basic charge distribution features (BasicQS) to create composite phosphorylation accessibility scores (phos_acc); implementing specialized functions for charged residue distance calculations and quintile-based interval plotting; generating comprehensive visualizations including interval bar plots showing relationships between AD stage-specific phosphorylation changes and structural features; comparative box plots analyzing structural feature distributions across phosphorylation change extremes; and conducting focused analysis of tau protein phosphorylation sites with specialized structural metrics. This multi-scale analysis provides insights into how structural and electrostatic environments influence phosphorylation dynamics in AD pathology, revealing potential structural determinants of phosphorylation changes in neurodegenerative disease progression and tau protein hyperphosphorylation.

library(tidyr)
library(dplyr)
library(ggpubr)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggdensity)
library(RColorBrewer)
library(ggrepel)
library(bio3d)
options(scipen=999)
options(stringsAsFactors = FALSE)
setwd("~/Other_project/AA_topology/reproducibility/f.PTM&AD")

source("../0.script/_plot_fig.R")
source("../0.script/_DN_function.R")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

load(paste0(data_folder,"AF2_AA_TopoTrend_dedup.rdata"))
topo_feature$AA<-bio3d::aa321(topo_feature$AA)

topo_feature$id<-paste0(topo_feature$AA,topo_feature$Pos)

load(paste0(data_folder,"AF2_Human_charged_Dist.rdata"))

phos_ab<-read.table(paste0(data_folder,"PTM_and_AD/AD_Phosphoproteome.csv"),quote="\"",header=T,sep=",",fill=TRUE)

phos_ab$adLFC<-round(log2((phos_ab$AD1+phos_ab$AD2)/(phos_ab$LPC1+phos_ab$LPC2)),4)
phos_ab$pspLFC<-round(log2((phos_ab$PSP1+phos_ab$PSP2)/(phos_ab$LPC1+phos_ab$LPC2)),4)
phos_ab$mciLFC<-round(log2((phos_ab$MCI1+phos_ab$MCI2)/(phos_ab$LPC1+phos_ab$LPC2)),4)
phos_ab$hpcLFC<-round(log2((phos_ab$HPC1+phos_ab$HPC2)/(phos_ab$LPC1+phos_ab$LPC2)),4)



Accession<-str_split_fixed(phos_ab$Protein.Accession,"\\|",n=3) 
phos_ab$Uniprot<-Accession[,2]

phos_ab_sep<-separate_longer_delim(phos_ab,c(Phosphosites), delim = ";")

phos_ab_sites<-unique(phos_ab_sep)


######
protein_in<-unique(phos_ab_sites$Uniprot)
####
charged_aa_dist<-charged_aa_dist[which(charged_aa_dist$Uniprot%in%protein_in),]

#######


phos_ab_topo<-inner_join(phos_ab_sites,topo_feature,by=join_by("Uniprot"=="Uniprot","Phosphosites"=="id"))

phos_ab_topo<-left_join(phos_ab_topo,charged_aa_dist[,c("Uniprot","id","logBasic","logAcidic","BasicQS","AcidicQS")],by=join_by("Uniprot"=="Uniprot","Phosphosites"=="id"))

phos_ab_topo<-data.frame(phos_ab_topo)
phos_ab_topo$phos_acc<-phos_ab_topo$CF10QS+phos_ab_topo$BasicQS

write.table(phos_ab_topo,
            paste0("AD_stage_phos_topo.csv"),
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)

######

df<-phos_ab_topo

p<-1
plot_list<-list()

tri_colsets<-matrix(rep(c("#377EB8","#B3CDE3","#4DAF4A","#CCEBC5","#E41A1C","#FBB4AE"),3),byrow = T,ncol=2)
for(xvar in c("hpcLFC","mciLFC","adLFC")){
  
  for(yvar in c("CF10QS","BasicQS","phos_acc")){
    
    plot_list[[p]]<-plot_bar_interval_10(df,xvar,yvar,yvar,tri_colsets[p,])
    p<-p+1
  }
  
}

pdf_file<-paste0(output_folder,"AD_stage_phos_topo_bar_intervalue.pdf")
p<-ggarrange(plotlist=plot_list,ncol=3,nrow=3)

ggsave(pdf_file,p,width=12,height=8,limitsize = FALSE)

######

for(yvar in c("CF10QS","phos_acc")){
  p<-1
  plot_list<-list()
  for(group_var in c("hpcLFC","mciLFC","adLFC")){
    
    plot_list[[p]]<-qu_boxplot(df,0.2,group_var,yvar,group_var)
    p<-p+1
  }
  
  pdf_file<-paste0(output_folder,"AD_stage_phos_",yvar,"_boxplot.pdf")
  p<-ggarrange(plotlist=plot_list,ncol=3,nrow=1,common.legend = T)
  
  ggsave(pdf_file,p,width=6,height=4,limitsize = FALSE)
  
}


#########


charged_dist_calc<-function(charged_dist){
  
  Basic_AA<-c("K","R","H")
  Acidic_AA<-c("D","E")
  Charged_AA<-c("K","R","H","D","E")
  
  aa_DN_value<-DN_sep2(charged_dist,TRUE,c(Basic_AA,Acidic_AA))
  
  charged_aa_dist<-charged_dist[,c("Uniprot","chain_len","Pos","AA")]
  charged_aa_dist$Basic<-matrixStats::rowMins(as.matrix(aa_DN_value[,Basic_AA]))
  charged_aa_dist$Acidic<-matrixStats::rowMins(as.matrix(aa_DN_value[,Acidic_AA]))
  
  charged_aa_dist$Charged<-matrixStats::rowMins(as.matrix(aa_DN_value[,Charged_AA]))
  
  charged_aa_dist$BasicQS<-round(rank(charged_aa_dist$Basic,na.last="keep")/length(charged_aa_dist$Basic),3)
  charged_aa_dist$AcidicQS<-round(rank(charged_aa_dist$Acidic,na.last="keep")/length(charged_aa_dist$Acidic),3)
  charged_aa_dist$ChargedQS<-round(rank(charged_aa_dist$Charged,na.last="keep")/length(charged_aa_dist$Charged),3)
  
  charged_aa_dist$id<-paste0(charged_aa_dist$AA,charged_aa_dist$Pos)
  charged_aa_dist
}


plot_bar_interval_5<-function(df,xvar,yvar,title_label,color_set){
  
  df$x<-df[,xvar]
  df$y<-df[,yvar]
  
  df<-df[which(!is.na(df$x)&!is.na(df$y)),]
  
  quantile_labels <- paste0("(", seq(0, 80, by = 20), "%,", seq(20, 100, by = 20), "%]")
  
  df_plot <- df %>%
    mutate(quantile_group = ntile(x, 5)) %>%
    group_by(quantile_group) %>%
    summarise(
      x_median = median(x),  
      y_mean = mean(y),       
      y_q1 = quantile(y, 0.25), # 
      y_q3 = quantile(y, 0.75)  # 
    ) %>%
    ungroup() %>%
    mutate(
      group_label = factor(quantile_group, 
                           levels = 1:5, 
                           labels = quantile_labels)
    )
  
  plot_title<-paste0(title_label, " ",nrow(df)," sites")
  
  p<-ggplot(df_plot, aes(x = x_median)) +
    labs(title=plot_title)+
    # ylim(c(y_min,y_max))+
    geom_errorbar(
      aes(ymin = y_q1, ymax = y_q3),
      width = diff(range(df_plot$x_median)) * 0.05,  
      color= color_set[2],
      size = 1
    ) +
    geom_point(aes(y = y_mean), size = 4,color= color_set[1]) +
    geom_line(aes(y = y_mean), linewidth = 1.2,color= color_set[1]) +
    labs(x = xvar, 
         y = yvar) +
    scale_y_continuous(
      expand = expansion(mult = c(0.1, 0.1))
    ) +
    # scale_fill_brewer(palette = "Dark2") +
    theme_minimal(base_size = 12) +
    theme(
      plot.margin = unit(c(20, 5, 20, 5), "points"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 10),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      plot.caption = element_text(color = "gray50", hjust = 0)
    )  
  p
}


topo_df<-read.table(paste0(data_folder,"PTM_and_AD/AA_TopoAttr_P10636-8_1.txt"),quote="\"",header=T,sep="\t",fill=TRUE)
topo_df$AA<-aa321(topo_df$AA)

phos<-read.table(paste0(data_folder,"PTM_and_AD/Tau_phosphosites.txt"),quote="\"",header=T,sep="\t",fill=TRUE)

dist<-read.table(paste0(data_folder,"PTM_and_AD/AA_DistNearest_P10636-8_1_A.csv"),quote="\"",header=T,sep=",",fill=TRUE)
dist$Uniprot<-"P10636-8"

charged_dist<-charged_dist_calc(dist)



phos$AA<-substr(phos$Tau_phosphosites,1,1)
phos$Pos<-substr(phos$Tau_phosphosites,2,nchar(phos$Tau_phosphosites))

phos$Pos<-as.numeric(phos$Pos)

phos_topo<-left_join(phos,topo_df,by=c("AA","Pos"))

phos_topo<-left_join(phos_topo,charged_dist[,c("AA","Pos","BasicQS","AcidicQS")],by=c("AA","Pos"))

phos_topo$phos_acc<-phos_topo$CF10QS+phos_topo$BasicQS

data_df<-phos_topo

tri_colsets<-matrix(c("#377EB8","#B3CDE3","#4DAF4A","#CCEBC5","#E41A1C","#FBB4AE"),byrow = T,ncol=2)

p<-1
plot_list<-list()
xvar<-"LFC"
for(yvar in c("CF10QS","BasicQS","phos_acc")){
  
  plot_list[[p]]<-plot_bar_interval_5(data_df,xvar,yvar,yvar,tri_colsets[p,])
  p<-p+1
}

pdf_file<-paste0(output_folder,"AD_tau_phos_topo_bar_intervalue.pdf")
p<-ggarrange(plotlist=plot_list,ncol=3,nrow=1)

ggsave(pdf_file,p,width=12,height=2.5,limitsize = FALSE)

write.table(phos_topo,
            "AD_tau_phos_topo.csv",
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)





