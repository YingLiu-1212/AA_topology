# This R script performs comparative survival analysis between patients carrying mutations in ASPD genes versus non-carriers across different cancer treatment modalities. The code implements a binary classification approach to compare progression-free survival (PFS) and overall survival (OS) between two patient groups: those with at least one mutation in ASPD genes (ASPD+) and those with no ASPD mutations (ASPD-). Key analytical steps include: processing TCGA patient data with treatment information and clinical outcomes; creating binary ASPD carrier status based on the presence or absence of mutations in ASPD genes; stratifying patients by treatment type (radiation therapy, chemotherapy, or all treatments combined) and tumor stage (early-stage I-III vs late-stage IV); generating Kaplan-Meier survival curves with confidence intervals and log-rank p-values; and producing comparative survival plots that clearly display case numbers for each group. 

plot_surv_pfs<-function(data,var,cut,pdf_prefix){
  
  data<-data[!is.na(data[,var]),]
  data<-data[!is.na(data$PFS_MONTHS),]

  data$CF_type<-""
  data$CF_type[which(data[,var]==0)]<-paste0(var,"- ",length(which(data[,var]==0))," cases")
  
  data$CF_type[which(data[,var]>0)]<-paste0(var,"+ ",length(which(data[,var]>0))," cases")
  
  plot_data<-data[which(data$CF_type!=""),]
  
  
  plot_data$PFS_STATUS<-as.numeric(str_split_fixed(plot_data$PFS_STATUS,":",n=2)[,1])
  
  
  fit <- survfit(Surv(PFS_MONTHS,PFS_STATUS) ~ CF_type, data = plot_data)
  
  plot<-ggsurvplot(fit, data = plot_data,conf.int = TRUE, pval = TRUE,title=paste0(pdf_prefix," ",nrow(plot_data)," patients"),palette = c("#377EB8","#E41A1C"), xlab="PFS_MONTHS")
  
  ggsave(paste0(output_folder,"Surv_PFS_",pdf_prefix,".pdf"), print(plot$plot), width = 6, height = 5)
}

plot_surv_os<-function(data,var,cut,pdf_prefix){
  
  data<-data[!is.na(data[,var]),]
  data<-data[!is.na(data$OS_MONTHS),]
  # 
  data$CF_type<-""
  data$CF_type[which(data[,var]==0)]<-paste0(var,"- ",length(which(data[,var]==0))," cases")
  
  data$CF_type[which(data[,var]>0)]<-paste0(var,"+ ",length(which(data[,var]>0))," cases")
  
  plot_data<-data[which(data$CF_type!=""),]
  
  plot_data<-plot_data[!is.na(plot_data$OS_MONTHS),]
  
  plot_data$OS_STATUS<-as.numeric(str_split_fixed(plot_data$OS_STATUS,":",n=2)[,1])
  
  fit <- survfit(Surv(OS_MONTHS,OS_STATUS) ~ CF_type, data = plot_data)
  
  plot<-ggsurvplot(fit, data = plot_data,conf.int = TRUE, pval = TRUE,title=paste0(pdf_prefix," ",nrow(plot_data)," patients"),palette = c("#377EB8","#E41A1C"),xlab="OS_MONTHS")
  
  ggsave(paste0(output_folder,"Surv_OS_",pdf_prefix,".pdf"), print(plot$plot), width = 6, height = 5)
}

library(survival)
library(survminer)
library(stringr)

library(dplyr)
setwd("~/Other_project/AA_topology/reproducibility/g.Cancer")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

load("tcga_atlas_2018_mut_topo.rdata")

gene_list<-read.table(paste0(data_folder,"ASPD/All.tsv"),quote="",header=T,sep="\t",fill=TRUE)

tcga_treatment_topo$ASPD<-0
tcga_treatment_topo$ASPD[which(tcga_treatment_topo$uniprot_id%in%gene_list$Entry)]<-1

tcga_treatment_topo<-tcga_treatment_topo[!duplicated(tcga_treatment_topo[,c("uniprot_id","HGVSp_Short","PATIENT_ID","TREATMENT_TYPE")]),]

patient_sum<- tcga_treatment_topo %>%
  group_by(PATIENT_ID,TREATMENT_TYPE) %>%
  summarise(CANCER_TYPE_ACRONYM=CANCER_TYPE_ACRONYM[1],
            ASPD = sum(ASPD),
            AJCC_PATHOLOGIC_TUMOR_STAGE = unique(AJCC_PATHOLOGIC_TUMOR_STAGE),
            OS_STATUS=OS_STATUS[1],
            OS_MONTHS=OS_MONTHS[1],
            PFS_STATUS=PFS_STATUS[1],
            PFS_MONTHS=PFS_MONTHS[1]
  )




###

cf_cut<-0.5

for(treatment in c("Radiation Therapy","Chemotherapy")){
  patient_sum_treatment<-patient_sum[which(patient_sum$TREATMENT_TYPE==treatment),]
  
  patient_sum_treatment_s123<-patient_sum_treatment[which(patient_sum_treatment$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("STAGE I","STAGE I/II (NOS)","STAGE IA","STAGE IB","STAGE II","STAGE IIA","STAGE IIB","STAGE IIC","STAGE III","STAGE IIIA","STAGE IIIB","STAGE IIIC")),]
  
  patient_sum_treatment_s4<-patient_sum_treatment[which(patient_sum_treatment$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("STAGE IV","STAGE IVA","STAGE IVB","STAGE IVC")),]
  
  cancer_type<-"S123"
  
  plot_surv_pfs(patient_sum_treatment_s123,"ASPD",cf_cut,paste0("TCGA_",treatment,"_ASPDcarr_",cancer_type))
  
  plot_surv_os(patient_sum_treatment_s123,"ASPD",cf_cut,paste0("TCGA_",treatment,"_ASPDcarr_",cancer_type))
  # 
  # cancer_type<-"S4"
  # 
  # plot_surv_pfs(patient_sum_treatment_s4,"ASPD",cf_cut,paste0("TCGA_",treatment,"_ASPDcarr_",cancer_type))
  # plot_surv_os(patient_sum_treatment_s4,"ASPD",cf_cut,paste0("TCGA_",treatment,"_ASPDcarr_",cancer_type))
  
}


treatment<-"All"
patient_sum_treatment<-patient_sum

patient_sum_treatment_s123<-patient_sum_treatment[which(patient_sum_treatment$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("STAGE I","STAGE I/II (NOS)","STAGE IA","STAGE IB","STAGE II","STAGE IIA","STAGE IIB","STAGE IIC","STAGE III","STAGE IIIA","STAGE IIIB","STAGE IIIC")),]

patient_sum_treatment_s4<-patient_sum_treatment[which(patient_sum_treatment$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("STAGE IV","STAGE IVA","STAGE IVB","STAGE IVC")),]

cancer_type<-"S123"

plot_surv_pfs(patient_sum_treatment_s123,"ASPD",cf_cut,paste0("TCGA_",treatment,"_ASPDcarr_",cancer_type))

plot_surv_os(patient_sum_treatment_s123,"ASPD",cf_cut,paste0("TCGA_",treatment,"_ASPDcarr_",cancer_type))

cancer_type<-"S4"

plot_surv_pfs(patient_sum_treatment_s4,"ASPD",cf_cut,paste0("TCGA_",treatment,"_ASPDcarr_",cancer_type))


plot_surv_os(patient_sum_treatment_s4,"ASPD",cf_cut,paste0("TCGA_",treatment,"_ASPDcarr_",cancer_type))
