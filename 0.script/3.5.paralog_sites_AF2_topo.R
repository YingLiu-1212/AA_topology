# This R script performs comparative structural analysis of paralogous proteins to investigate how amino acid conservation patterns correlate with protein topological features. The code integrates multiple data sources including paralog sequence alignments and AlphaFold2-predicted structural metrics to analyze structural differences between conserved positions, missense variants, and indel sites in human paralog pairs with >70% sequence identity. Key analytical steps include: filtering paralog pairs by sequence identity threshold; merging topological feature data (CF10QS, LD15QS, CF10QS.Trend4bi) with paralog conservation annotations; classifying positions into conserved, missense, or indel categories. The script generates comprehensive visualizations comparing topological feature distributions across different conservation categories, including both global analyses across all amino acids and focused examinations of specific residue types (S, A, V, L, T). Output includes multi-panel PDF figures that systematically reveal how evolutionary conservation correlates with structural environments, providing insights into structural constraints on protein evolution. The analysis leverages custom plotting functions from external R scripts for consistent visualization formatting.

library(ggrepel)
# library(msa)
# library(Biobase)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(bio3d)
library(tidyr)
options(scipen=999)
options(stringsAsFactors = FALSE)

setwd("~/Other_project/AA_topology/reproducibility/c.Global_analysis")

source("../0.script/_plot_fig.R")
source("../0.script/_DN_function.R")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

paralog_nonConsv<-read.table(paste0(data_folder,"UP000005640_paralog_nonConsv_AA.txt"),quote="\"",header=T,sep="\t",fill=TRUE)

idcut<-70

idty_cut<-paralog_nonConsv[which(paralog_nonConsv$identity>idcut),]
idty_cut_prot<-unique(idty_cut$target)
length(idty_cut_prot)

# data_folder<-"~/Published_data/AlphaFold/Topology/"
load(paste0(data_folder,"AF2_AA_TopoTrend_dedup.rdata"))
topo_feature$AA<-bio3d::aa321(topo_feature$AA)

topo_df_para<-topo_feature[which(topo_feature$Uniprot%in%idty_cut_prot),]

topo_df_para$id<-paste0(topo_df_para$Pos,topo_df_para$AA)

para_info_comb<-left_join(topo_df_para,idty_cut,by=join_by("Uniprot"=="target","Pos"=="aa_pos","AA"=="aa_from"))

para_info_comb$AA_paralog<-"Conserved"
para_info_comb$AA_paralog[which(para_info_comb$aa_to=="-")]<-"Indel"
para_info_comb$AA_paralog[which((!is.na(para_info_comb$aa_to))&(para_info_comb$aa_to!="-"))]<-"Missense"

para_info_comb$id<-paste0(para_info_comb$Uniprot,"_",para_info_comb$id)

para_info_comb<-para_info_comb[order(para_info_comb$AA_paralog,decreasing = T),]
para_info_comb<-para_info_comb[!duplicated(para_info_comb$id),]

para_info_comb$id<-paste0(para_info_comb$Pos,para_info_comb$AA)

#######


load(paste0(data_folder,"AF2_WG_charged_DN_Norm.rdata"))

charged_aa_dist$id<-paste0(charged_aa_dist$Pos,charged_aa_dist$AA)

para_info_comb<-left_join(para_info_comb,charged_aa_dist[,c("Uniprot","id","logBasic","logAcidic","logCharged","BasicQS","AcidicQS","ChargedQS")],by=join_by("Uniprot","id"))
#######


# source("~/Other_project/AA_topology/script/_plot_fig.R")

plot_df<-para_info_comb

table(para_info_comb$AA_paralog)

p<-1
plot_list<-list()
for(plot_x in c("CF10QS","LD15QS","CF10QS.Trend4bi")){
  
  plot_y<-"AA_paralog"
  
  
  x_min<-floor(min(plot_df[,plot_x],na.rm = T))
  x_max<-ceiling(max(plot_df[,plot_x],na.rm = T))
  
  plot_list[[p]]<-plot_bar_even_10_logy(plot_df,plot_x,plot_y,x_min,x_max)
  p<-p+1
  
}

pdf_file<-paste0(output_folder,"Paralog70_AF2_prop_dist.pdf")
plots<-ggarrange(plotlist=plot_list,ncol=3,nrow=1,common.legend = T)

ggsave(pdf_file,plots,width=10,height=3,limitsize = FALSE)

######

p<-1
plot_list<-list()
AA_type<-c("S","A","V","L","T")
for(i in 1:length(AA_type)){
  
  
  plot_y<-"AA_paralog"
  
  
  for(plot_x in c("CF10QS","LD15QS","CF10QS.Trend4bi")){
    para_AA<-plot_df[which(plot_df$AA==AA_type[i]),]
    
    plot_list[[p]]<-plot_bar_even_10_logy(para_AA,plot_x,plot_y,x_min,x_max)
    p<-p+1
    
  }
  
  
}

pdf_file<-paste0(output_folder,"Paralog70_AF2_AA_prop_dist.pdf")
plots<-ggarrange(plotlist=plot_list,ncol=3,nrow=5,common.legend = T,labels = rep(AA_type,each=3))

ggsave(pdf_file,plots,width=10,height=12,limitsize = FALSE)
