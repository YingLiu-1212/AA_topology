# This R script performs comprehensive survival analysis to evaluate the prognostic significance of ASPD gene mutation counts in cancer patients receiving different treatment modalities. The code implements Kaplan-Meier survival analysis to compare progression-free survival (PFS) and overall survival (OS) between patient groups with high versus low ASPD mutation counts. Key analytical features include: processing TCGA patient data with treatment information and clinical outcomes; calculating ASPD mutation counts per patient; stratifying patients by treatment type (radiation therapy, chemotherapy, or all treatments combined) and tumor stage (early-stage I-III vs late-stage IV); generating survival curves with confidence intervals and log-rank p-values; and producing publication-ready survival plots.

plot_surv_pfs<-function(data,var,cut,pdf_prefix){
  
  data<-data[!is.na(data[,var]),]
  data<-data[!is.na(data$PFS_MONTHS),]
  
  low_cut<-quantile(data[,var],cut)
  high_cut<-quantile(data[,var],1-cut)
  
  print(low_cut)
  
  data$CF_type<-""
  data$CF_type[which(data[,var]<=low_cut)]<-paste0(var,"_Low")
  
  data$CF_type[which(data[,var]>=high_cut)]<-paste0(var,"_High")
  
  plot_data<-data[which(data$CF_type!=""),]
  
  
  plot_data$PFS_STATUS<-as.numeric(str_split_fixed(plot_data$PFS_STATUS,":",n=2)[,1])
  
  
  fit <- survfit(Surv(PFS_MONTHS,PFS_STATUS) ~ CF_type, data = plot_data)
  
  plot<-ggsurvplot(fit, data = plot_data,conf.int = TRUE, pval = TRUE,title=paste0(pdf_prefix," ",nrow(plot_data)," patients"),palette = c("#E41A1C","#377EB8"), xlab="PFS_MONTHS")
  
  ggsave(paste0(output_folder,"Surv_PFS_",pdf_prefix,".pdf"), print(plot$plot), width = 6, height = 5)
}

plot_surv_os<-function(data,var,cut,pdf_prefix){
  
  data<-data[!is.na(data[,var]),]
  data<-data[!is.na(data$PFS_MONTHS),]
  
  low_cut<-quantile(data[,var],cut)
  high_cut<-quantile(data[,var],1-cut)
  
  data$CF_type<-""
  data$CF_type[which(data[,var]<=low_cut)]<-paste0(var,"_Low")
  
  data$CF_type[which(data[,var]>=high_cut)]<-paste0(var,"_High")
  
  plot_data<-data[which(data$CF_type!=""),]
  
  plot_data<-plot_data[!is.na(plot_data$OS_MONTHS),]
  
  plot_data$OS_STATUS<-as.numeric(str_split_fixed(plot_data$OS_STATUS,":",n=2)[,1])
  
  fit <- survfit(Surv(OS_MONTHS,OS_STATUS) ~ CF_type, data = plot_data)
  
  plot<-ggsurvplot(fit, data = plot_data,conf.int = TRUE, pval = TRUE,title=paste0(pdf_prefix," ",nrow(plot_data)," patients"),palette = c("#E41A1C","#377EB8"), xlab="OS_MONTHS")
  
  ggsave(paste0(output_folder,"Surv_OS_",pdf_prefix,".pdf"), print(plot$plot), width = 6, height = 5)
}

library(survival)
library(survminer)
library(stringr)
library(stringr)
library(dplyr)
setwd("~/Other_project/AA_topology/reproducibility/g.Cancer")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

load("tcga_atlas_2018_mut_topo.rdata")

gene_list<-read.table(paste0(data_folder,"ASPD/All.tsv"),quote="",header=T,sep="\t",fill=TRUE)
tcga_treatment_topo<-tcga_treatment_topo[which(tcga_treatment_topo$uniprot_id%in%gene_list$Entry),]
tcga_treatment_topo$ASPD<-1

tcga_treatment_topo<-tcga_treatment_topo[!duplicated(tcga_treatment_topo[,c("uniprot_id","HGVSp_Short","PATIENT_ID","TREATMENT_TYPE")]),]


patient_sum<- tcga_treatment_topo %>%
  group_by(PATIENT_ID,TREATMENT_TYPE) %>%
  summarise(CANCER_TYPE_ACRONYM=CANCER_TYPE_ACRONYM[1],
            ASPD_count = sum(ASPD),
            AJCC_PATHOLOGIC_TUMOR_STAGE = unique(AJCC_PATHOLOGIC_TUMOR_STAGE),
            OS_STATUS=OS_STATUS[1],
            OS_MONTHS=OS_MONTHS[1],
            PFS_STATUS=PFS_STATUS[1],
            PFS_MONTHS=PFS_MONTHS[1]
  )

patient_sum<-data.frame(patient_sum)

# write.table(patient_sum,
#             "TCGA_patient_ADPS_sum.csv",
#             append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)



###

## c("Radiation Therapy","Chemotherapy")
cf_cut<-0.5


for(treatment in c("Radiation Therapy","Chemotherapy")){
  patient_sum_treatment<-patient_sum[which(patient_sum$TREATMENT_TYPE==treatment),]
  
  
  patient_sum_treatment_s123<-patient_sum_treatment[which(patient_sum_treatment$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("STAGE I","STAGE I/II (NOS)","STAGE IA","STAGE IB","STAGE II","STAGE IIA","STAGE IIB","STAGE IIC","STAGE III","STAGE IIIA","STAGE IIIB","STAGE IIIC")),]
  
  
  patient_sum_treatment_s4<-patient_sum_treatment[which(patient_sum_treatment$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("STAGE IV","STAGE IVA","STAGE IVB","STAGE IVC")),]
  
  cancer_type<-"S123"
  
  plot_surv_pfs(patient_sum_treatment_s123,"ASPD_count",cf_cut,paste0("TCGA_",treatment,"_Count_",cancer_type))
  
  plot_surv_os(patient_sum_treatment_s123,"ASPD_count",cf_cut,paste0("TCGA_",treatment,"_Count_",cancer_type))
  
  # cancer_type<-"S4"
  # 
  # plot_surv_pfs(patient_sum_treatment_s4,"ASPD_count",cf_cut,paste0("TCGA_",treatment,"_Count_",cancer_type))
  # 
  # plot_surv_os(patient_sum_treatment_s4,"ASPD_count",cf_cut,paste0("TCGA_",treatment,"_Count_",cancer_type))
  
  
}

treatment<-"All"
patient_sum_treatment<-patient_sum

patient_sum_treatment_s123<-patient_sum_treatment[which(patient_sum_treatment$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("STAGE I","STAGE I/II (NOS)","STAGE IA","STAGE IB","STAGE II","STAGE IIA","STAGE IIB","STAGE IIC","STAGE III","STAGE IIIA","STAGE IIIB","STAGE IIIC")),]


patient_sum_treatment_s4<-patient_sum_treatment[which(patient_sum_treatment$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("STAGE IV","STAGE IVA","STAGE IVB","STAGE IVC")),]

cancer_type<-"S123"

plot_surv_pfs(patient_sum_treatment_s123,"ASPD_count",cf_cut,paste0("TCGA_",treatment,"_Count_",cancer_type))

plot_surv_os(patient_sum_treatment_s123,"ASPD_count",cf_cut,paste0("TCGA_",treatment,"_Count_",cancer_type))

cancer_type<-"S4"

plot_surv_pfs(patient_sum_treatment_s4,"ASPD_count",cf_cut,paste0("TCGA_",treatment,"_Count_",cancer_type))

plot_surv_os(patient_sum_treatment_s4,"ASPD_count",cf_cut,paste0("TCGA_",treatment,"_Count_",cancer_type))

