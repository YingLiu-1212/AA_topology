# This R script performs comprehensive ROC analysis to evaluate the predictive power of mutation counts in specific gene families for cancer treatment outcomes across different tumor stages. The code analyzes TCGA data to assess how missense mutation counts in four key gene families (ADAM, SCNA, PCDHB, DOCK) can predict progression-free survival in patients receiving radiation therapy or chemotherapy. The analysis is stratified by tumor stage (Stage I, Stages I-III combined, and Stage IV) to investigate stage-specific predictive patterns. For each gene family and tumor stage combination, the script calculates ROC curves using missense mutation counts as predictors and progression-free survival status as the outcome, then generates comparative ROC plots with AUC values and case numbers. 

library(matrixStats)
library(reshape2)
library(ggplot2)
library(viridis)
library(pROC)
library(dplyr)
library(ggpubr)
library(stringr)
library(RColorBrewer)
options(scipen=999)
options(stringsAsFactors = FALSE)

setwd("~/Other_project/AA_topology/reproducibility/g.Cancer")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

load("tcga_atlas_2018_mut_topo.rdata")

tcga_topo<-tcga_treatment_topo[which(!is.na(tcga_treatment_topo$PFS_MONTHS)),]

tcga_treatment_topo<-tcga_treatment_topo[tcga_treatment_topo$TREATMENT_TYPE%in%c("Radiation Therapy","Chemotherapy"),]


tcga_treatment_topo<-tcga_treatment_topo[!duplicated(tcga_treatment_topo[,c("uniprot_id","HGVSp_Short","PATIENT_ID","TREATMENT_TYPE")]),]


tcga_treatment_topo_s123<-tcga_treatment_topo[which(tcga_treatment_topo$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("STAGE I","STAGE I/II (NOS)","STAGE IA","STAGE IB","STAGE II","STAGE IIA","STAGE IIB","STAGE IIC","STAGE III","STAGE IIIA","STAGE IIIB","STAGE IIIC")),]


tcga_treatment_topo_s1<-tcga_treatment_topo[which(tcga_treatment_topo$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("STAGE I","STAGE I/II (NOS)","STAGE IA","STAGE IB")),]

tcga_treatment_topo_s4<-tcga_treatment_topo[which(tcga_treatment_topo$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("STAGE IV","STAGE IVA","STAGE IVB","STAGE IVC")),]


calc_ROC_compare<-function(group_treatment_CF,treatment_type){
  
  group_treatment_CF<-group_treatment_CF[!is.na(group_treatment_CF$CF10QS),]
  
  if(treatment_type!=""){
    group_treatment_type_CF<-group_treatment_CF[which(group_treatment_CF$TREATMENT_TYPE==treatment_type),]
  }else{
    group_treatment_type_CF<-group_treatment_CF
  }
  
  group_treatment_type_CF<-group_treatment_type_CF[!duplicated(group_treatment_type_CF[,c("uniprot_id","HGVSp_Short","PATIENT_ID")]),]
  
  group_treatment_type_CF<-group_treatment_type_CF[!is.na(group_treatment_type_CF$PFS_STATUS),]
  group_treatment_type_CF<-group_treatment_type_CF[group_treatment_type_CF$PFS_STATUS!="",]
  
  
  group_sum<- group_treatment_type_CF %>%
    group_by(PATIENT_ID) %>%
    summarise(
      # CF10_sum = sum(CF10QS),
      #         revLD15_sum = sum(1-LD15QS),
      Missense_count = length(CF10QS),
      PFS_STATUS=PFS_STATUS[1],
      PFS_MONTHS=PFS_MONTHS[1]
    )
  
  # print(head(group_sum))
  
  roc_count <- roc(response = group_sum$PFS_STATUS,predictor = group_sum$Missense_count)  
  
  print("count")
  print(roc_count$auc)
  
  roc_count
}

##########

plot_ROC_group<-function(treatment_type,tcga_treatment_topo,pdf_prefix){
  
  family_file<-c("ADAM","SCNA","PCDHB","DOCK","All")
  color_set<-c(
    "#984EA3","#A65628","#377EB8","#4DAF4A","#ED1C24")
  # color_set<-brewer.pal(8, 'Dark2')
  
  #treatment_type<-"Chemotherapy"
  
  count_auc<-rep(0,length(family_file))
  ncases<-rep(0,length(family_file))
  pdf(file=paste0(pdf_prefix,".pdf"), width=5, height=5)
  
  for(i in 1:length(family_file)){
    
    group<-family_file[i]
    
    group_gene<-read.table(paste0(data_folder,"ASPD/",group,".tsv"),quote="",header=T,sep="\t",fill=TRUE)
    
    group_tcga_treatment_topo<-tcga_treatment_topo[which(tcga_treatment_topo$uniprot_id%in%c(group_gene$Entry)),]
    
    
    if(treatment_type!=""){
      ncases[i]<-length(unique(group_tcga_treatment_topo[which(group_tcga_treatment_topo$TREATMENT_TYPE==treatment_type),"PATIENT_ID"]))
    }else{
      ncases[i]<-length(unique(group_tcga_treatment_topo[,"PATIENT_ID"]))
    }
    
    
    
    roc<-calc_ROC_compare(group_tcga_treatment_topo,treatment_type)
    
    
    if(i==1){
      plot_count_roc<-plot(roc, col=color_set[i], main=treatment_type,legacy.axes=T,lwd=1) 
      count_auc[i]<-round(plot_count_roc$auc,2)
      
    }else{
      plot_count_roc<-plot(roc, add=TRUE, col=color_set[i],lwd=1)
      count_auc[i]<-round(plot_count_roc$auc,2)
    }
  }
  
  legend("bottomright",box.lwd = 0,
         legend=c(paste0(family_file,"(",ncases," Cases) AUC=",round(count_auc,2))),cex = 0.7, col=color_set,lty=1)
  
  dev.off()
  
  
  
}

####

#

treatment_type<-"Radiation Therapy"
pdf_prefix<-paste0(output_folder,"TCGA_ASPD_ROC_S123_",treatment_type)
plot_ROC_group(treatment_type,tcga_treatment_topo_s123,pdf_prefix)

pdf_prefix<-paste0(output_folder,"TCGA_ASPD_ROC_S4_",treatment_type)
plot_ROC_group(treatment_type,tcga_treatment_topo_s4,pdf_prefix)

pdf_prefix<-paste0(output_folder,"TCGA_ASPD_ROC_S1_",treatment_type)
plot_ROC_group(treatment_type,tcga_treatment_topo_s1,pdf_prefix)




####




