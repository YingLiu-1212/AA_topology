# This R script performs comprehensive analysis of cancer mutation structural topology in relation to radiotherapy treatment outcomes using TCGA (The Cancer Genome Atlas) data. The code investigates how protein structural features influence progression-free survival (PFS) in cancer patients receiving radiation therapy. Key analytical components include: processing TCGA mutation data with structural topology annotations; implementing an interval-based analysis method that divides structural feature values into multiple bins and calculates the distribution of disease progression versus censored cases within each interval; computing correlation coefficients between progression and censored case proportions across structural feature intervals to identify genes where structural context correlates with treatment outcomes; conducting both pan-cancer analysis and focused analysis of genes from the ASPD collection; generating multi-panel visualization plots that show the relationship between structural feature intervals and clinical outcomes for individual genes; and identifying genes with negative correlations where specific structural environments are associated with differential treatment responses. This analysis provides insights into how protein structural context may influence radiotherapy efficacy and identifies potential structural biomarkers for treatment response prediction in cancer therapy.

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
#library(pROC)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
options(scipen=999)
options(stringsAsFactors = FALSE)
setwd("~/Other_project/AA_topology/reproducibility/g.Cancer")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

load("tcga_atlas_2018_mut_topo.rdata")

tcga_topo<-tcga_treatment_topo[which(tcga_treatment_topo$TREATMENT_TYPE=="Radiation Therapy"),]

tcga_topo<-tcga_topo[which(!is.na(tcga_treatment_topo$PFS_MONTHS)),]

tcga_topo<-tcga_topo[!duplicated(tcga_topo[,c("uniprot_id","HGVSp_Short","PATIENT_ID")]),]

mut_count<-table(tcga_topo$Hugo_Symbol)
mut_count<-mut_count[order(mut_count,decreasing = T)]
# rownames(prot_gene)<-prot_gene$uniprot_id

nintervals<-5
ncut<-50
subCut<-20
color_sets<-brewer.pal(9,"Set1")

# prot_list<-names(mut_count[mut_count>ncut])
gene_list<-names(mut_count[mut_count>ncut])
length(gene_list)

# gene_list<-gene_list[1:1000]

response_info_all<-NULL

for(g in 1:length(gene_list)){
  # if(g%%100==0)print(g)
  gene_name<-gene_list[g]
  
  
  prot_topo_CF<-tcga_topo[which(tcga_topo$Hugo_Symbol==gene_name),]
  prot_id<-prot_topo_CF$uniprot_id[1]
  
  ##progression-free survival,PFS
  ##disease-free survival,DFS
  
  prot_topo_CF<-prot_topo_CF[!is.na(prot_topo_CF$CF10QS),]
  
  prot_topo_CF<-prot_topo_CF[!is.na(prot_topo_CF$PFS_MONTHS),]
  
  prot_topo_CF$PFS_STATUS<-factor(prot_topo_CF$PFS_STATUS,levels=c("1:PROGRESSION","0:CENSORED"))
  
  ncases<-length(unique(prot_topo_CF$PATIENT_ID))
  
  if(ncases>ncut & min(table(prot_topo_CF$PFS_STATUS))>subCut & length(table(prot_topo_CF$PFS_STATUS))>1){
    
    plot_title<-paste0(gene_name)
    
    info1<-plot_intervals(prot_topo_CF,"CF10QS",nintervals,plot_title)
    info2<-plot_intervals(prot_topo_CF,"LD15QS",nintervals,plot_title)  
    info3<-plot_intervals(prot_topo_CF,"CF10QS.Trend4bi",nintervals,plot_title)  
    
    CF10QS_prop<-info1[[2]]
    LD15QS_prop<-info2[[2]]
    
    CF10QS_cor<-round(cor(CF10QS_prop$Proportion.x,CF10QS_prop$Proportion.y),3)
    LD15QS_cor<-round(cor(LD15QS_prop$Proportion.x,LD15QS_prop$Proportion.y),3)
    
    response_info<-data.frame(Cancer="All",ncases=ncases,gene_name,prot_id,CF10QS_cor,LD15QS_cor)
    
    response_info_all<-rbind(response_info_all,response_info)
    
  } ## if(ncases>10)
  
  
  
}


write.table(response_info_all,
            paste0("TCGA_RD_topo_intervals_response.csv"),
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)

############
gene_list<-read.table(paste0(data_folder,"ASPD/All.tsv"),quote="",header=T,sep="\t",fill=TRUE)


tcga_topo<-tcga_topo[which(tcga_topo$uniprot_id%in%gene_list$Entry),]

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
    
    if(CF10QS_cor<0){
      
      plot_CF[[j_1]]<-info1[[1]]
      
      j_1<-j_1+1
    }
    if(LD15QS_cor<0){
      
      plot_LD[[j_2]]<-info2[[1]]

      j_2<-j_2+1
    }
  } ## if(ncases>10)

}



p_cf<-ggarrange(plotlist=plot_CF,nrow=7,ncol=4,common.legend = T)

ggsave(paste0(output_folder,"TCGA_ASPD_gene_CFinterval_RD.pdf"),p_cf,width=12,height=21)


p_ld<-ggarrange(plotlist=plot_LD,nrow=6,ncol=4,common.legend = T)

ggsave(paste0(output_folder,"TCGA_ASPD_gene_LDinterval_RD.pdf"),p_ld,width=12,height=18)

############
