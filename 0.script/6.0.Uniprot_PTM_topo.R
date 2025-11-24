# This R script analyzes the relationship between post-translational modifications (PTMs) and protein structural topology using AlphaFold2-predicted structural features. The code integrates comprehensive PTM annotation data from UniProt with structural topology metrics to investigate how different types of modifications are distributed across structural environments. Key analytical steps include: processing and filtering PTM data for major modification types including phosphorylation (serine, threonine, tyrosine) and lysine modifications (acetylation, succinylation); integrating PTM site information with structural topology features (CF10QS, LD15QS, CF10QS.Trend4bi); performing quality control to ensure amino acid compatibility between modification types and residue contexts; and generating comparative distribution plots that show how PTM sites are distributed across structural feature intervals compared to non-modified residues. The analysis produces multi-panel visualizations that systematically compare the structural contexts of different PTM types, providing insights into how protein structural environment influences or is influenced by post-translational modifications, with potential implications for understanding PTM regulation, functional consequences, and structural determinants of modification sites.

library(tidyr)
library(dplyr)
library(ggpubr)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggdensity)
library(RColorBrewer)
library(ggrepel)
# library(Biostrings)
options(scipen=999)
options(stringsAsFactors = FALSE)

setwd("~/Other_project/AA_topology/reproducibility/f.PTM&AD")

source("../0.script/_plot_fig.R")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

PTM_all<-read.table(paste0(data_folder,"PTM_and_AD/uniprotkb_modified_residue_info.txt"),quote="\"",header=T,sep=",",fill=TRUE)

PTM_prot<-unique(PTM_all$Entry)
#######

load(paste0(data_folder,"AF2_AA_TopoTrend_dedup.rdata"))
topo_feature$AA<-bio3d::aa321(topo_feature$AA)
topo_feature$id<-paste0(topo_feature$Pos,topo_feature$AA)

topo_feature<-topo_feature[which(topo_feature$Uniprot%in%PTM_prot),]

topo_feature$Pos<-as.character(topo_feature$Pos)
############

topo_PTM<-left_join(topo_feature,PTM_all,by=join_by("Uniprot"=="Entry","Pos"))

topo_PTM$PTM<-topo_PTM$Type
length(which(!is.na(topo_PTM$PTM)))
topo_PTM$PTM[is.na(topo_PTM$PTM)]<-"Unknown"

topo_PTM<-topo_PTM[!duplicated(topo_PTM[,c("pdb_chain","id","PTM")]),]

PTM_use<-c("N6-acetyllysine","N6-succinyllysine","Phosphoserine","Phosphothreonine","Phosphotyrosine")

topo_PTM$PTM[which(!topo_PTM$PTM%in%c(PTM_use,"Unknown"))]<-"Others"


err_anno<-which((topo_PTM$Type=="N6-acetyllysine" & topo_PTM$AA!="K") |
                  (topo_PTM$Type=="N6-succinyllysine" & topo_PTM$AA!="K") |
                  (topo_PTM$Type=="Phosphoserine" & topo_PTM$AA!="S") |
                  (topo_PTM$Type=="Phosphothreonine" & topo_PTM$AA!="T") |
                  (topo_PTM$Type=="Phosphotyrosine" & topo_PTM$AA!="Y"))

print(length(err_anno))

topo_PTM<-topo_PTM[-err_anno,]
topo_PTM$PTM<-topo_PTM$Type
topo_PTM$PTM[is.na(topo_PTM$PTM)]<-"Other"
####

PTM_type<-c("N6-acetyllysine","N6-succinyllysine","Phosphoserine","Phosphothreonine","Phosphotyrosine")
AA_type<-c("K","K","S","T","Y")
# plot_x<-"CF10QS"

length(unique(topo_PTM$Uniprot))
nrow(topo_PTM)

# topo_PTM$phos_acc<-topo_PTM$CF10QS+topo_PTM$BasicQS
#c("CF10QS","LD15QS","CF10QS.Trend4bi","logBasic","logAcidic","BasicQS","AcidicQS")
for(plot_x in c("CF10QS","LD15QS","CF10QS.Trend4bi")){
  p<-1
  plot_list<-list()
  
  PTM_single<-topo_PTM[which(topo_PTM$AA%in%AA_type),]
  
  x_min<-floor(min(PTM_single[,plot_x],na.rm = T))
  x_max<-ceiling(max(PTM_single[,plot_x],na.rm = T))
  
  PTM_single$Modification<-"Unknown"
  PTM_single$Modification[!is.na(PTM_single$Type)]<-"Modification"
  PTM_single$Modification<-factor(PTM_single$Modification,levels=c("Modification","Unknown"))
  
  plot_y<-"Modification"
  
  plot_list[[p]]<-plot_bar_even_10(PTM_single,plot_x,plot_y,x_min,x_max)
  p<-p+1
  
  plot_y<-"PTM"
  for(i in 1:length(PTM_type)){
    
    PTM_single<-topo_PTM[which(topo_PTM$PTM%in%c(PTM_type[i],"Other") & topo_PTM$AA==AA_type[i]),]
    
    PTM_single$PTM<-factor(PTM_single$PTM,levels=c(PTM_type[i],"Other"))
    
    plot_list[[p]]<-plot_bar_even_10(PTM_single,plot_x,plot_y,x_min,x_max)
    p<-p+1
    
  }
  
  
  
  pdf_file<-paste0(output_folder,"PTM_AF2_",plot_x,"_prop_dist.pdf")
  p<-ggarrange(plotlist=plot_list,ncol=3,nrow=2)
  
  ggsave(pdf_file,p,width=9,height=5,limitsize = FALSE)
  ####
  
  
  
  
}
