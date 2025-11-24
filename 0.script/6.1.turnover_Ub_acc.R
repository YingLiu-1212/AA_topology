# This R script performs comprehensive analysis of ubiquitination sites in relation to protein structural topology, charge distribution, and protein degradation kinetics. The code integrates experimental ubiquitination data with AlphaFold2-predicted structural features to investigate how structural context influences ubiquitin-mediated protein degradation. Key analytical components include: processing ubiquitination site data with associated protein half-life measurements and occupancy percentages; integrating structural topology metrics (CF10QS) with charge distribution features (AcidicQS) to create composite structural-charge accessibility scores (ub_acc); generating multiple visualization types including interval bar plots showing relationships between protein half-life and structural features across decile groups; comparative box plots analyzing structural feature distributions across proteasome activity clusters (DOWN, NC, UP) and half-life categories (VERY FAST, FAST, SLOW, VERY SLOW); and 2D density plots examining the bivariate distribution of core-forming propensity and acidic charge accessibility across different degradation rate categories. This multi-faceted analysis provides insights into how structural and electrostatic properties of ubiquitination sites correlate with protein turnover rates, proteasome processing efficiency, and degradation kinetics, revealing potential structural determinants of ubiquitin-mediated protein degradation.

library(dplyr)
library(ggpubr)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggdensity)
library(RColorBrewer)
library(ggrepel)
library(earth)
options(scipen=999)
options(stringsAsFactors = FALSE)
setwd("~/Other_project/AA_topology/reproducibility/f.PTM&AD")

source("../0.script/_plot_fig.R")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

load(paste0(data_folder,"AF2_AA_TopoTrend_dedup.rdata"))
topo_feature$AA<-bio3d::aa321(topo_feature$AA)

topo_feature$id<-paste0(topo_feature$AA,topo_feature$Pos)

load(paste0(data_folder,"AF2_Human_charged_Dist.rdata"))


#####

ubl_hl<-read.table(paste0(data_folder,"PTM_and_AD/ubiquitylation_sites_all.csv"),quote="\"",header=T,sep=",",fill=TRUE)
summary(as.numeric(ubl_hl$half.life..min.))
ubl_hl$half.life..min.<-as.numeric(ubl_hl$half.life..min.)

ubl_hl$LogHL<-round(log2(ubl_hl$half.life..min.+0.1),3)

ubl_hl_reduce<-ubl_hl[which(!is.na(ubl_hl$LogHL)),]

######
protein_in<-unique(ubl_hl_reduce$Protein.IDs)

topo_feature<-topo_feature[which(topo_feature$Uniprot%in%protein_in),]

charged_aa_dist<-charged_aa_dist[which(charged_aa_dist$Uniprot%in%protein_in),]

######
ubl_hl_topo<-inner_join(ubl_hl_reduce,topo_feature,by=join_by("Protein.IDs"=="Uniprot","Position"=="id"))

ubl_hl_topo<-left_join(ubl_hl_topo,charged_aa_dist[,c("Uniprot","id","logBasic","logAcidic","AcidicQS","BasicQS")],by=join_by("Protein.IDs"=="Uniprot","Position"=="id"))


ubl_hl_topo$ub_acc<-ubl_hl_topo$CF10QS+ubl_hl_topo$AcidicQS


ubl_hl_topo$occupancy_percentage<-as.numeric(gsub("%","",ubl_hl_topo$PCdG.occupancy))
ubl_hl_topo$LogOccup<-log2(0.001+ubl_hl_topo$occupancy_percentage)

write.table(ubl_hl_topo,
            paste0("UBL_AF2_HL.csv"),
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)


length(unique(ubl_hl_topo$Protein.IDs))

##########
# 
# plot_x<-"half.life..min."
plot_x<-"LogHL"
p<-1
plot_list<-list()

tri_colsets<-matrix(c("#377EB8","#B3CDE3","#4DAF4A","#CCEBC5","#E41A1C","#FBB4AE"),byrow = T,ncol=2)

for(feature in c("CF10QS","AcidicQS","ub_acc")){
  
  plot_list[[p]]<-plot_bar_interval_10(ubl_hl_topo,plot_x,feature,feature,tri_colsets[p,])

  p<-p+1
}

pdf_file<-paste0(output_folder,"UBL_AF2_",plot_x,"_bar_intervalue.pdf")
pcomb<-ggarrange(plotlist=plot_list,ncol=3,nrow=1)

ggsave(pdf_file,pcomb,width=12,height=3,limitsize = FALSE)

#######

plot_x<-"LogHL"
p<-1
plot_box<-list()
for(feature in c("CF10QS","AcidicQS","ub_acc")){
  plot_box[[p]]<-qu_boxplot(ubl_hl_topo,0.1,plot_x,feature,feature)
  
  p<-p+1
}

pdf_file<-paste0(output_folder,"UBL_AF2_",plot_x,"_box.pdf")
p<-ggarrange(plotlist=plot_box,ncol=3,nrow=1,common.legend=T)
ggsave(pdf_file,p,width=6,height=4,limitsize = FALSE)

#############

plot_df<-ubl_hl_topo[which(ubl_hl_topo$occupancy_percentage>0),]
plot_x<-"LogOccup"
p<-1
plot_box<-list()

for(feature in c("CF10QS","AcidicQS","ub_acc")){
  plot_box[[p]]<-qu_boxplot(plot_df,0.1,plot_x,feature,feature)
  
  p<-p+1
}

pdf_file<-paste0(output_folder,"UBL_AF2_",plot_x,"_box.pdf")
p<-ggarrange(plotlist=plot_box,ncol=3,nrow=1,common.legend=T)
ggsave(pdf_file,p,width=6,height=4,limitsize = FALSE)

######



plot_box<-list()
p<-1
for(plotVar in c("CF10QS","AcidicQS","ub_acc")){
  
  plot_data<-ubl_hl_topo[which(!is.na(ubl_hl_topo$ProteasomeCluster)),]
  
  my_comparisons<-list(c("DOWN","NC"),c("NC","UP"),c("DOWN","UP"))
  
  color_set<-brewer.pal(4, 'Dark2')
  
  plot_box[[p]]<-ggplot(plot_data,aes(x=factor(ProteasomeCluster),y= .data[[plotVar]]))+
    labs(title=stringr::str_wrap(plotVar,width=40),size=1)+
    xlab("")+
    # ylim(c(min(plot_data[,plotVar])-0.1,max(plot_data[,plotVar])+0.3))+
    geom_boxplot(aes(color=ProteasomeCluster),
                 width = 0.6,
                 alpha = 0.85,
                 lwd = 0.8, 
                 fatten = 1.2)+
    scale_fill_manual(values = palette_box) +
    scale_color_manual(values= palette_box) +
    stat_compare_means(aes(label="p.signif"),comparisons = my_comparisons,method="wilcox.test")+
    # theme_bw() +
    theme_minimal() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size=7))
  
  p<-p+1
  
  
}

pdf_file<-paste0(output_folder,"UBL_AF2_ProteasomeCluster_box.pdf")
p<-ggarrange(plotlist=plot_box,ncol=3,nrow=1,common.legend=T)
ggsave(pdf_file,p,width=8,height=4,limitsize = FALSE)


#########


plot_box<-list()
p<-1
for(plotVar in c("CF10QS","AcidicQS","ub_acc")){
  
  plot_data<-ubl_hl_topo[which(!is.na(ubl_hl_topo$HL.Cluster)),]
  
  plot_data$HL.Cluster<-factor(plot_data$HL.Cluster,levels=c("VERY FAST","FAST","SLOW","VERY SLOW"))
  
  my_comparisons<-list(c("VERY FAST","SLOW"),c("FAST","VERY SLOW"),c("VERY FAST","VERY SLOW"))
  
  color_set<-brewer.pal(4, 'Dark2')
  
  plot_box[[p]]<-ggplot(plot_data,aes(x=factor(HL.Cluster),y= .data[[plotVar]]))+
    labs(title=stringr::str_wrap(plotVar,width=40),size=1)+
    xlab("")+
    # ylim(c(min(plot_data[,plotVar])-0.1,max(plot_data[,plotVar])+0.3))+
    geom_boxplot(aes(color=HL.Cluster),                   width = 0.6,
                 alpha = 0.85,
                 lwd = 0.8, 
                 fatten = 1.2)+
    stat_compare_means(aes(label="p.signif"),comparisons = my_comparisons,method="wilcox.test")+
    theme_minimal() +
    scale_color_manual(values=palette_box)+ 
    scale_fill_manual(values = palette_box) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size=7))
  
  p<-p+1
  
  
}

pdf_file<-paste0(output_folder,"UBL_AF2_HL.Cluster_box.pdf")
p<-ggarrange(plotlist=plot_box,ncol=3,nrow=1,common.legend=T)
ggsave(pdf_file,p,width=8,height=3.5,limitsize = FALSE)



#######

target_1<-"CF10QS"
target_2<-"AcidicQS"
p_list<-list()
k<-1

for(ubl_type in c("VERY SLOW","SLOW","FAST","VERY FAST")){
  
  plot_data<-ubl_hl_topo[which(ubl_hl_topo$HL.Cluster==ubl_type),]
  plot_title<-ubl_type
  p_list[[k]]<-plot_dens2(plot_data,plot_title,target_1,target_2)
  k<-k+1
}


p_comb<-ggarrange(plotlist=p_list,nrow=1,ncol=4)
pdf_file<-paste0(output_folder,"UBL_AF2_HL.Cluster_dens.pdf")
ggsave(pdf_file,p_comb,width=12,height=3)
