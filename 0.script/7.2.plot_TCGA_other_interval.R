# This R script performs focused analysis of specific cancer-related genes in TCGA data to investigate the relationship between protein structural topology and radiotherapy treatment outcomes. The code analyzes a curated set of seven key cancer genes (including KRAS, PIK3CA, PTEN, and TP53) to examine how structural features correlate with progression-free survival in patients receiving radiation therapy. Using an interval-based analytical approach, the script divides structural feature values (CF10QS and LD15QS) into multiple bins and calculates the distribution of disease progression versus censored cases within each interval. The analysis specifically identifies genes with positive correlations between structural feature intervals and clinical outcomes, generating multi-panel visualization plots that show how mutation distribution patterns across different structural environments relate to treatment response. This targeted approach enables detailed investigation of structural determinants of treatment outcomes in well-characterized cancer genes, providing insights into how specific protein structural contexts may influence radiotherapy efficacy in molecularly defined patient subgroups.

plot_intervals<-function(data_CF,feature,nintervals,plot_title){
  
  df<-data_CF
  breaks <- seq(-1, 1, length.out = nintervals+1) 
  groups <- cut(df[,feature], breaks = breaks, include.lowest = TRUE)  
  df$interval<-groups
  
  count_feat<-table(df[,c("interval","PFS_STATUS")])
  
  plot_data<-melt(count_feat,value.name ="nCases")
  
  plot_data$Proportion<-0
  
  pfs0<-which(plot_data$PFS_STATUS=="0:CENSORED")
  plot_data$Proportion[pfs0]<-round(plot_data$nCases[pfs0]/sum(plot_data$nCases[pfs0]),4)*100
  pfs1<-which(plot_data$PFS_STATUS=="1:PROGRESSION")
  plot_data$Proportion[pfs1]<-round(plot_data$nCases[pfs1]/sum(plot_data$nCases[pfs1]),4)*100
  
  #print(plot_data)
  
  Proportion0<-plot_data[pfs0,c("interval","Proportion")]
  Proportion1<-plot_data[pfs1,c("interval","Proportion")]
  Proportion_comb<-left_join(Proportion0,Proportion1,by="interval")
  
  #print(Proportion_comb)
  
  diff<-round(cor(Proportion_comb$Proportion.x,Proportion_comb$Proportion.y),3)
  
  palette_box <-c("#E41A1C","#377EB8")
  
  p1<-ggplot(plot_data, aes(x = factor(interval), y = Proportion,group=PFS_STATUS,colour=PFS_STATUS)) +  
    ggtitle(paste0(plot_title))+
    xlab(feature)+
    geom_line(alpha=0.8) +
    geom_point(aes(size=nCases),alpha=0.8)+
    scale_fill_manual(values = palette_box) +
    scale_color_manual(values= palette_box) +
    annotate('text',x=4,y=13,
             label=paste0("Correlation = ",diff))+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
  
  
  list(p1,Proportion_comb)
  
}


library(matrixStats)
library(reshape2)
library(ggplot2)
library(viridis)
library(stringr)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
options(scipen=999)
options(stringsAsFactors = FALSE)

setwd("~/Other_project/AA_topology/reproducibility/g.Cancer")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"


load("tcga_atlas_2018_mut_topo.rdata")

tcga_topo<-tcga_treatment_topo[which(!is.na(tcga_treatment_topo$PFS_MONTHS)),]

tcga_topo<-tcga_treatment_topo[which(tcga_treatment_topo$TREATMENT_TYPE=="Radiation Therapy"),]


tcga_topo<-tcga_topo[!duplicated(tcga_topo[,c("uniprot_id","HGVSp_Short","PATIENT_ID")]),]

gene_list<-c("P01116","O75874","P60484","Q14315","P42336","P01111","P04637")

tcga_topo<-tcga_topo[which(tcga_topo$uniprot_id%in%gene_list),]

length(unique(tcga_topo$uniprot_id))


mut_count<-table(tcga_topo$Hugo_Symbol)
mut_count<-mut_count[order(mut_count,decreasing = T)]


nintervals<-5
ncut<-50
subCut<-20
color_sets<-brewer.pal(9,"Set1")

gene_list<-names(mut_count[mut_count>ncut])
length(gene_list)

response_info_all<-NULL

plot_CF<-list()
plot_LD<-list()
j_1<-1
j_2<-1

for(g in 1:length(gene_list)){
  gene_name<-gene_list[g]
  patient_treatment_CF<-tcga_topo[which(tcga_topo$Hugo_Symbol==gene_name),]
  prot_id<-patient_treatment_CF$uniprot_id[1]
  
  ##progression-free survival,PFS
  ##disease-free survival,DFS
  
  patient_treatment_CF<-patient_treatment_CF[!is.na(patient_treatment_CF$CF10QS),]
  
  patient_treatment_CF<-patient_treatment_CF[!is.na(patient_treatment_CF$PFS_MONTHS),]
  
  patient_treatment_CF$PFS_STATUS<-factor(patient_treatment_CF$PFS_STATUS,levels=c("1:PROGRESSION","0:CENSORED"))
  
  ncases<-length(unique(patient_treatment_CF$PATIENT_ID))
  
  if(ncases>ncut & min(table(patient_treatment_CF$PFS_STATUS))>subCut & length(table(patient_treatment_CF$PFS_STATUS))>1){
    
    plot_title<-paste0(gene_name," ",ncases," cases")
    
    info1<-plot_intervals(patient_treatment_CF,"CF10QS",nintervals,plot_title)
    info2<-plot_intervals(patient_treatment_CF,"LD15QS",nintervals,plot_title)
    
    
    CF10QS_prop<-info1[[2]]
    LD15QS_prop<-info2[[2]]
    
    CF10QS_cor<-round(cor(CF10QS_prop$Proportion.x,CF10QS_prop$Proportion.y),3)
    
    LD15QS_cor<-round(cor(LD15QS_prop$Proportion.x,LD15QS_prop$Proportion.y),3)
    
    if(CF10QS_cor>0){
      
      plot_CF[[j_1]]<-info1[[1]]
      
      j_1<-j_1+1
    }
    if(LD15QS_cor>0){
      
      plot_LD[[j_2]]<-info2[[1]]

      j_2<-j_2+1
    }
    
    
  } ## if(ncases>10)

}

#main_title<-"TCGA S123 progression by CF10QS"
length(plot_CF)

p_cf<-ggarrange(plotlist=plot_CF,nrow=3,ncol=3,common.legend = T)

ggsave(paste0(output_folder,"TCGA_other_gene_CFinterval_RD.pdf"),p_cf,width=9,height=9)

p_ld<-ggarrange(plotlist=plot_LD,nrow=3,ncol=3,common.legend = T)

ggsave(paste0(output_folder,"TCGA_other_gene_LDinterval_RD.pdf"),p_ld,width=9,height=9)


