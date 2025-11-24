# This R script generates 2D density visualization plots for protein structural topology features across the human proteome. The code loads AlphaFold2-predicted topological data for human proteins and creates comparative density plots between core topological metrics (CF10QS) and their derived trend features (CF10QS.Trend4bi). Using custom visualization functions from the _plot_fig.R script, it produces bivariate density contour plots that reveal the distribution patterns and correlations between different structural parameters at the proteome scale. The analysis encompasses the entire human proteome dataset, visualizing the relationship between core packing density metrics and their spatial trend variations across all amino acid positions. Output includes multi-panel PDF figures that provide comprehensive overviews of structural feature distributions, enabling identification of common patterns and outliers in protein structural space. This visualization supports large-scale structural bioinformatics analyses by revealing global trends in protein architecture and facilitating quality assessment of structural predictions.

library(dplyr)
library(ggpubr)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggdensity)
library(RColorBrewer)
library(ggrepel)
options(scipen=999)
options(stringsAsFactors = FALSE)

setwd("~/Other_project/AA_topology/reproducibility/d.SFI")
source("../0.script/_plot_fig.R")

topo_used<-"CF10QS"

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

load(paste0(data_folder,"AF2_AA_TopoTrend_dedup.rdata"))

topo_additional<-paste0(topo_used,c(".Trend4bi"))
topo_all<-c("CF10QS",topo_additional)
topo_pair<-combn(topo_all,2)


p_list<-list()
i<-1
plot_data<-topo_feature

for(j in 1:ncol(topo_pair)){
  
  target_1<-topo_pair[1,j]
  target_2<-topo_pair[2,j]
  
  plot_title<-paste0("Human Whole Proteomics (",nrow(plot_data)," AA)")
  
  p_list[[i]]<-plot_dens(plot_data,plot_title,target_1,target_2)
  
  i<-i+1
  
}

p<-ggarrange(plotlist=p_list,ncol=ncol(topo_pair),nrow=1)

pdf_file<-paste0(output_folder,"AF2_Human_Proteome_Dens.pdf")

ggsave(pdf_file,p,width=ncol(topo_pair)*3,height=3,limitsize = FALSE)


