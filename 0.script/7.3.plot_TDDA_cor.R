# This R script generates a comprehensive scatter plot visualization to analyze the relationship between structural topology features and radiotherapy treatment outcomes across cancer genes from TCGA data. The code integrates correlation coefficients calculated from interval-based analysis of progression-free survival with gene family. It specifically focuses on four key gene families (ADAM, SCNA, PCDHB, DOCK) involved in cell adhesion, signaling, and extracellular matrix interactions. The visualization highlights genes that show significant negative correlations in either CF10QS or LD15QS structural features with treatment outcomes, using color-coding to distinguish gene families and intelligent label placement to avoid overlaps. The resulting plot provides a global overview of how different structural topology metrics correlate with radiotherapy response across diverse gene families, enabling identification of potential structural biomarkers for treatment resistance and insights into family-specific patterns of structural determinants in cancer therapy response.

library(matrixStats)
library(reshape2)
library(ggplot2)
library(viridis)
#library(pROC)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)
options(scipen=999)
options(stringsAsFactors = FALSE)

setwd("~/Other_project/AA_topology/reproducibility/g.Cancer")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

###########

response_info_all<-read.table(paste0("TCGA_RD_topo_intervals_response.csv"),quote="",header=T,sep=",",fill=TRUE)

family_file<-c("ADAM","SCNA","PCDHB","DOCK")
all_gene<-NULL
for(group in family_file){
  
  group_gene<-read.table(paste0(data_folder,"ASPD/",group,".tsv"),quote="",header=T,sep="\t",fill=TRUE)
  
  group_gene$family<-group
  
  all_gene<-rbind(all_gene,group_gene)
}

# response<-response_info_all[which(response_info_all$ncases>100),]
response<-response_info_all
print(nrow(response))

response<-left_join(response,all_gene[,c("Entry","family")],by=join_by("prot_id"=="Entry"))

col_sets<-data.frame(ZZ="#00000050",PCDHB="#377EB8",DOCK="#4DAF4A",ADAM="#984EA3",SCNA="#A65628")


hot_gene<-which(response$prot_id%in%all_gene$Entry & ((response$LD15QS_cor<0 & response$CF10QS_cor<0.5)| (response$CF10QS_cor<0 & response$LD15QS_cor<0.5)))

length(hot_gene)

response$label<-""
response$label[hot_gene]<-response$gene_name[hot_gene]

response$type<-"ZZ"
response$type[hot_gene]<-response$family[hot_gene]

response$size<-0.01
response$size[which(response$type!="ZZ")]<-10

response<-response[order(response$type,decreasing = T),]

plot3<-ggplot(data=response,aes(LD15QS_cor,CF10QS_cor,colour=type)) +
  xlim(-1.1,1.1)+ylim(-1.1,1.1)+
  #xlab("Radiation Therapy")+ylab("Chemotherapy")+
  #labs(title="CF10QS_cor")+ 
  geom_point(aes(size=size))+ 
  scale_size_area(max_size = 2) +
  #geom_abline(lwd=0.1)+
  geom_hline(yintercept = 0,linetype=3,lwd=0.1)+
  geom_vline(xintercept = 0,linetype=3,lwd=0.1)+
  geom_text_repel(aes(LD15QS_cor,CF10QS_cor,label = label),show.legend = FALSE,max.overlaps =10000,min.segment.length = 0)+
  scale_colour_manual(values=col_sets)+
  theme_minimal()+
  theme(legend.position = "none")


ggsave(paste0(output_folder,"TCGA_RD_Cor_ECM_highlight.pdf"),plot3,width=4,height=4)

