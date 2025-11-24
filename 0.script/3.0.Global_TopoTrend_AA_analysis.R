# This R script performs comprehensive analysis of amino acid distribution patterns across protein structural features using topological data from human proteins. The code processes AlphaFold2-predicted structural metrics to investigate how different amino acid types are distributed along specific structural parameters. Key functionalities include: classification of amino acids into seven biochemical groups (hydrophobic, polar, basic, acidic, aromatic, sulfhydryl, and unique) through a custom mapping function; processing of large-scale topological data from PDB structures; and generation of detailed distribution profiles across structural feature intervals. The analysis calculates proportional distributions of amino acids across binned structural values and produces two types of visualizations: comprehensive plots showing all amino acid types, and focused plots highlighting polar residues. The script generates publication-ready line plots with labeled data points using ggplot2, with enhanced visualization through color-coding by biochemical groups and intelligent label placement via ggrepel to minimize overlaps. Output includes both numerical proportion tables and PDF figures, enabling systematic investigation of amino acid preference patterns in protein structural contexts.

convert_aa_to_group<-function(aa_list){
  
  aa_group<-aa_list
  aa_group[which(aa_group%in%c("A","V","L","I","M"))]<-"Hydrophobic"
  aa_group[which(aa_group%in%c("S","T","N","Q"))]<-"Polar"
  aa_group[which(aa_group%in%c("K","R","H"))]<-"Basic"
  aa_group[which(aa_group%in%c("D","E"))]<-"Acidic"
  aa_group[which(aa_group%in%c("F","W","Y"))]<-"Aromatic"  
  aa_group[which(aa_group=="C")]<-"Sulfhydryl"  
  aa_group[which(aa_group%in%c("G","P"))]<-"Unique"  
  
  aa_group
  
}
library(bio3d)
library(dplyr)
library(ggpubr)
library(stringr)
library(reshape2)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggrepel)
options(scipen=999)
options(stringsAsFactors = FALSE)

setwd("~/Other_project/AA_topology/reproducibility/c.Global_analysis")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

load(paste0(data_folder,"uniprotkb_pdb_human_Topology_dedup.rdata"))

topo_trend$AA<-aa321(topo_trend$AA)
topo_trend<-topo_trend[which(!is.na(topo_trend$AA)),]

# AA_property<-read.table("../0.Data/AA_property2.txt",quote="",header=T,sep="\t",fill=TRUE)


# length(unique(topo_trend$pdb_id)) 
# length(unique(topo_trend$pdb_chain)) 
# length(unique(topo_trend$Uniprot)) 

AA_count<-table(topo_trend$AA)
AA_count<-cbind(AA=names(AA_count),count=AA_count,percent=round(AA_count/sum(AA_count),4))

for(feature in c("LD15QS")){
  
  aa_feat<-topo_trend[,c(feature,"AA")]
  aa_feat<-aa_feat[which(!is.na(aa_feat[,feature])),]
  values<-aa_feat[,feature]
  
  if(min(values)<0){
    breaks <- seq(-1, 1, length.out = 11) 
  }else{
    breaks <- seq(0, 1, length.out = 11) 
  }
  
  
  groups <- cut(values, breaks = breaks, include.lowest = TRUE)  
  aa_feat$interval<-groups
  
  count_feat<-table(aa_feat[,c(2,3)])
  # write.table(count_feat,
  #             paste0("AA_count_",feature,".csv"),
  #             append = F, quote = T, sep = ",",row.names = T, col.names = NA)
  
  # prop_feat<-round(count_feat/rowSums(count_feat),4)*100
  
  
  prop_feat<-round(count_feat/nrow(aa_feat),4)*100
  write.table(prop_feat,
              paste0("AA_prop_",feature,".csv"),
              append = F, quote = T, sep = ",",row.names = T, col.names = NA)
  
  #######
  
  plot_data<-melt(prop_feat)
  plot_data$AA<-as.character(plot_data$AA)
  
  color_sets<-brewer.pal(7,"Dark2")
  #######
  
  
  plot_data$label<-as.character(plot_data$AA)
  plot_data$label[which(!plot_data$interval%in%c("[-1,-0.8]","(0.8,1]","[0,0.1]","(0.9,1]"))]<-""
  
  
  plot_data$group<-convert_aa_to_group(plot_data$AA)
  
  plot_data$group<-factor(plot_data$group,levels=c("Hydrophobic","Aromatic","Polar","Basic","Acidic","Sulfhydryl","Unique"))
  
  p1<-ggplot(plot_data, aes(x = factor(interval),y=value,colour=group,group=AA)) +
    geom_line() +
    geom_point(size=2)+
    guides(color=guide_legend(nrow=10, byrow=TRUE)) +
    scale_x_discrete(expand = c(0.2, 0.2))+
    geom_label_repel(aes(label=label,colour=group),
                     show.legend = FALSE,
                     min.segment.length = 0,
                     max.overlaps=100,
                     max.iter = 1e5, max.time = 1)+
    scale_colour_manual(values=color_sets)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y = element_text()) +
    labs(title = paste0("AA type distribution along ",feature),x = "", y = "Proportion (%)")
  
  ggsave(paste0(output_folder,"PDB_AA_prop_",feature,".pdf"),p1,width=8,height=3)
  
  ########
  
  plot_data2<-plot_data[which(plot_data$AA%in%c("S","T","N","Q","C")),]
  
  
  p2<-ggplot(plot_data2, aes(x = factor(interval),y=value,colour=group,group=AA)) +
    geom_line() +
    geom_point(size=2)+
    guides(color=guide_legend(nrow=10, byrow=TRUE)) +
    scale_x_discrete(expand = c(0.2, 0.2))+
    geom_label_repel(aes(label=label,colour=group),
                     show.legend = FALSE,
                     min.segment.length = 0,
                     max.overlaps=100,
                     max.iter = 1e5, max.time = 1)+
    scale_colour_manual(values=color_sets)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y = element_text()) +
    labs(title = paste0("AA type distribution along ",feature),x = "", y = "Proportion (%)")
  
  ggsave(paste0(output_folder,"PDB_Polar_AA_prop_",feature,".pdf"),p2,width=8,height=3)
  
  
}







