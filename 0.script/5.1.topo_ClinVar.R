# This R script performs comprehensive analysis of ClinVar missense variants in relation to protein structural topology features predicted by AlphaFold2. The code integrates clinical variant annotations from ClinVar with structural topology metrics to investigate how protein structural context influences clinical variant interpretation. Key analytical steps include: processing and standardizing ClinVar missense mutation data; classifying variants into clinical significance categories (Pathogenic/Likely Pathogenic, Benign/Likely Benign, Uncertain Significance); integrating structural topology features (CF10QS, LD15QS, CF10QS.Trend4bi) for variant positions; generating distribution plots of structural features across clinical categories; performing receiver operating characteristic (ROC) analysis to evaluate the predictive power of topological features for clinical classification; and conducting stratified analysis by protein structural fragility index (SFI) to examine how predictive performance varies across proteins with different structural properties. The analysis produces multiple visualization outputs including distribution plots, density comparisons, ROC curves, and fragility-stratified performance plots, providing systematic insights into the relationship between protein structural environment and clinical variant pathogenicity assessment.

plot_classification2<-function(plot_data,use_feature,use_label,pdf_prefix){
  
  my_comparisons <- list(c("Likely Benign & Benign","Likely Pathogenic & Pathogenic"))
  
  p1<-ggplot(plot_data,aes(x=factor(Classification),y={{ use_feature }}))+
    labs(title=stringr::str_wrap(use_label,width=40),size=3)+ 
    ylim(c(0,1.3))+
    xlab("")+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(shape=16, position = position_jitter(0.2),alpha=0.5)+
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)),comparisons = my_comparisons,label.y=1.2)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10))
  
  ggsave(paste0(pdf_prefix,"_",use_label,".pdf"),p1,width=4,height=6)
}

library(ggplot2)
library(dplyr)
library(ggpubr)
library(stringr)
# library(Biostrings)
library(pROC)
library(ggdensity)
library(RColorBrewer)
options(scipen=999)
options(stringsAsFactors = FALSE)
setwd("~/Other_project/AA_topology/reproducibility/e.Clinical_SNVs")
source("../0.script/_plot_fig.R")

data_folder<-"../a.resouce/"
SFI_folder<-"../d.SFI/"
output_folder<-"../o.output_figures/"


load(paste0(data_folder,"AF2_AA_TopoTrend_dedup.rdata"))
topo_feature$AA<-bio3d::aa321(topo_feature$AA)

topo_feature$id<-paste0(topo_feature$AA,topo_feature$Pos)

# source("~/Other_project/AA_topology/script/_DN_function.R")
# load(paste0(af2_folder,"AF2_WG_charged_DN_Norm.rdata"))
####

ClinVar_data<-read.table(paste0(data_folder,"UP000005640_clinvar_missense.txt"),quote="\"",header=T,sep=",",row.names = NULL,fill=TRUE)
use_prot<-unique(ClinVar_data$uniprot_id)

topo_feature<-topo_feature[which(topo_feature$Uniprot%in%use_prot),]

ClinVar_CF<-inner_join(ClinVar_data,topo_feature,by=join_by("uniprot_id"=="Uniprot","aa_pos"=="Pos","aa_from"=="AA"))
# 
# ClinVar_CF<-inner_join(ClinVar_CF,charged_aa_dist[,c("Uniprot","id","BasicQS","AcidicQS")],by=join_by("uniprot_id"=="Uniprot","id"=="id"))

ClinVar_CF<-ClinVar_CF[!is.na(ClinVar_CF$CF10),]

ClinVar_CF$ClinicalSignificance<-tolower(ClinVar_CF$ClinicalSignificance)

clinical_class<-table(ClinVar_CF$ClinicalSignificance)

clinical_class<-clinical_class[order(clinical_class,decreasing = T)]

ClinVar_CF$ClinicalClass<-"Uncertain"
ClinVar_CF$ClinicalClass[grepl("benign",ClinVar_CF$ClinicalSignificance)]<-"Benign or Likely benign"
ClinVar_CF$ClinicalClass[grepl("pathogenic",ClinVar_CF$ClinicalSignificance)]<-"Pathogenic or Likely pathogenic"



ClinVar_CF$ClinicalClass[grepl("conflicting",ClinVar_CF$ClinicalSignificance)]<-"Uncertain"

write.table(ClinVar_CF,
            paste0("ClinVar_AF2.csv"),
            append = FALSE, quote = T, sep = ",",row.names = F, col.names = T)

table(ClinVar_CF$ClinicalClass)
nrow(ClinVar_CF)


#####

plot_var<-"LD15QS"
p<-plot_bar_even_10(ClinVar_CF,plot_var,"ClinicalClass",-1,1)

pdf_file<-paste0(output_folder,"ClinVar_",plot_var,"_ClinicalClass_distr.pdf")

ggsave(pdf_file,p,width=6,height=3.2)

plot_var<-"CF10QS"
p<-plot_bar_even_10(ClinVar_CF,plot_var,"ClinicalClass",-1,1)

pdf_file<-paste0(output_folder,"ClinVar_",plot_var,"_ClinicalClass_distr.pdf")

ggsave(pdf_file,p,width=6,height=3.2)


#######

# tmp<-ClinVar_CF[,c("ClinicalSignificance","ClinicalClass")]

predict_clinvar<-ClinVar_CF
predict_clinvar<-predict_clinvar[which(predict_clinvar$ClinicalClass!="Uncertain"),]

predict_clinvar<-predict_clinvar[order(predict_clinvar$ClinicalClass,decreasing = T),]

predict_clinvar<-predict_clinvar[!duplicated(predict_clinvar[,c("uniprot_id","aa_pos","ClinicalClass")]),]


# plot_var<-"CF10QS"
# plot_var<-"LD15QS"

pdf_file<-paste0(output_folder,"ClinVar_AF2_dens")

type_order<-c("Pathogenic or Likely pathogenic","Benign or Likely benign")
plot_dens_combn(predict_clinvar,type_var="ClinicalClass",type_order,"CF10QS",2,pdf_file)


#######
# color_set<-brewer.pal(8, 'Dark2')
color_set<-c(
  "#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00"
  ,"#FFFF33","#A65628","#F781BF","#999999")

plot_df<-predict_clinvar


plot_feature<-c("CF10QS","LD15QS","CF10QS.Trend4bi")
auc_value<-rep(0,length(plot_feature))

pdf(file=paste0(output_folder,"ClinVar_AF2_topo_ROC_all.pdf"), width=5, height=5)

for(i in 1:length(plot_feature)){
  
  roc_data<-roc(response = predict_clinvar$ClinicalClass,predictor = predict_clinvar[,plot_feature[i]])
  
  if(i==1){
    
    plot_roc<-plot(roc_data, col=color_set[i], main="AA  topo estimate clinical effects",legacy.axes=T,lwd=0.6)
    auc_value[i]<-round(plot_roc$auc,2)
  }else{
    plot_roc<-plot(roc_data, add=TRUE, col=color_set[i],lwd=0.6)
    auc_value[i]<-round(plot_roc$auc,2)
  }
  
}

legend("bottomright",box.lwd = 0,
       legend=c(paste0(plot_feature," AUC=",round(auc_value,2))),cex = 0.7, col=color_set,lty=1)

dev.off()

print(auc_value)
########


##################
protein_FI<-read.table(paste0(SFI_folder,"AF2_WG_Proteins_Fragile.csv"),header=T,sep=",",quote = "",fill=TRUE)

ncut<-10
protein_FI$FI_group <-paste0("FI_",ntile(protein_FI$LogG1G3_prop, ncut))
quantile_labels<-paste0("FI_",c(1:ncut))

predict_clinvar<-inner_join(predict_clinvar,protein_FI,by=join_by("uniprot_id"=="Uniprot"))

# plot_x<-"LogG1G3_prop"

table(predict_clinvar$FI_group)


auc_list<-NULL
for(pred_var in plot_feature){
  
  auc<-rep(0,ncut)
  for(i in 1:length(quantile_labels)){
    
    group_df<-predict_clinvar[which(predict_clinvar$FI_group==quantile_labels[i]),]
    
    roc_data<-roc(response = group_df$ClinicalClass,predictor = group_df[,pred_var])
    
    auc[i]<-roc_data$auc
  }
  
  auc_list<-cbind(auc_list,auc)
  
}

colnames(auc_list)<-plot_feature

plot_data<-data.frame(auc_list)
plot_data$FI_group<-quantile_labels

plot_data$FI_group<-factor(plot_data$FI_group,levels=quantile_labels)

plot_data<-melt(plot_data,id.vars = "FI_group")


colnames(plot_data)<-c("FI_group","Predictor","AUC")

p<-ggplot(plot_data, aes(x = FI_group,fill=Predictor)) +
  geom_point(aes(x = FI_group,y = AUC,color=Predictor), size = 3) +
  geom_line(aes(x = FI_group,y = AUC, group = Predictor,color=Predictor), linewidth = 1.2) +
  scale_color_manual(values=color_set)+ 
  labs(title = "") +
  theme_minimal() +
  scale_x_discrete(limits = quantile_labels) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

pdf_file<-paste0(output_folder,"ClinVar_AF2_SFI_ROC.pdf")

ggsave(pdf_file,p,width=6,height=4,limitsize = FALSE)



