# This R script performs comprehensive analysis of clinical genetic variants from the ClinGen database in relation to protein structural topology features predicted by AlphaFold2. The code systematically processes clinical variant annotations across multiple gene panels to investigate how structural properties differentiate pathogenic from benign mutations. Key analytical components include: parsing and standardizing ClinGen variant data with amino acid code conversion from three-letter to one-letter codes; integrating clinical classifications with AlphaFold2-predicted structural metrics (CF10QS, LD15QS, CF10QS.Trend4bi); generating comparative boxplot visualizations of structural features across clinical significance categories (Pathogenic, Benign, Uncertain Significance) for multiple gene panels; and performing receiver operating characteristic (ROC) analysis to evaluate the predictive power of topological features for clinical variant classification. The analysis produces multi-panel comparative plots across different genes and comprehensive ROC curves assessing the diagnostic value of structural metrics. This work provides insights into how protein structural context influences clinical variant interpretation and demonstrates the potential of structural topology features as complementary evidence for variant pathogenicity assessment in clinical genetics.

convert_aa_code <- function(three_letter_code) {  
  
  amino_acid_map <- c(  
    "Ala" = "A", "Arg" = "R", "Asn" = "N", "Asp" = "D",  
    "Cys" = "C", "Glu" = "E", "Gln" = "Q", "Gly" = "G",  
    "His" = "H", "Ile" = "I", "Leu" = "L", "Lys" = "K",  
    "Met" = "M", "Phe" = "F", "Pro" = "P", "Ser" = "S",  
    "Thr" = "T", "Trp" = "W", "Tyr" = "Y", "Val" = "V",  
    "Sec" = "U", "Pyl" = "O"
  )  
  
  one_letter_code <- amino_acid_map[three_letter_code] 
  
  one_letter_code
}

plot_classification2<-function(plot_data,use_feature,use_label,pdf_prefix){
  
  # color_set<-brewer.pal(8, 'Dark2')  
  color_set<-c(
    "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00"
    ,"#FFFF33","#A65628","#F781BF","#999999")
  
  my_comparisons <- list(c("Likely Benign & Benign","Likely Pathogenic & Pathogenic"))
  
  p1<-ggplot(plot_data,aes(x=factor(Classification),y={{ use_feature }}))+
    labs(title=stringr::str_wrap(use_label,width=40),size=2)+ 
    ylim(c(-1.1,1.3))+
    xlab("")+
    ylab(pdf_prefix)+
    geom_boxplot(aes(color=Classification),outlier.shape = NA,                 width = 0.6,
                 alpha = 0.85,
                 lwd = 0.8, 
                 fatten = 1.2 )+
    geom_jitter(aes(color=Classification),shape=16, position = position_jitter(0.2),alpha=0.5)+
    scale_color_manual(values=color_set)+ 
    scale_fill_manual(values=color_set)+ 
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)),comparisons = my_comparisons,label.y=1.2)+
    theme_minimal() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size=15))
  
  p1
  
  # ggsave(paste0(pdf_prefix,"_",use_label,".pdf"),p1,width=4,height=6)
}
library(pROC)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(stringr)
# library(Biostrings)
library(RColorBrewer)
options(scipen=999)
options(stringsAsFactors = FALSE)
setwd("~/Other_project/AA_topology/reproducibility/e.Clinical_SNVs")
source("../0.script/_plot_fig.R")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

load(paste0(data_folder,"AF2_AA_TopoTrend_dedup.rdata"))
topo_feature$AA<-bio3d::aa321(topo_feature$AA)
topo_feature$id<-paste0(topo_feature$Pos,topo_feature$AA)

prot_gene<-read.table(paste0(data_folder,"UP000005640_9606_prot_gene.txt"),quote="\"",header=T,sep="\t",fill=TRUE)


pannel_info<-as.matrix(read.table(paste0(data_folder,"ClinGen/pannel_list.txt"),quote="\"",header=F,sep="\t",fill=TRUE))

ClinGen_all<-NULL

plot_cf<-list()
plot_ld<-list()
m<-1
for(i in 1:length(pannel_info)){
  
  
  pannel<-pannel_info[i]
  
  Clin_data<-read.table(paste0(data_folder,"ClinGen/",pannel,".csv"),quote="\"",header=T,sep=",",fill=TRUE)
  
  ncol(Clin_data)
  
  if(ncol(Clin_data)==1){
    Clin_data<-read.table(paste0(data_folder,"ClinGen/",pannel,".csv"),quote="\"",header=T,sep="\t")
  }  ## if(ncol(Clin_data)==1)
  
  tmp<-str_split_fixed(Clin_data[,1],"\\(",n=3)
  tmp2<-str_split_fixed(tmp[,2],"\\):",n=2)
  tmp[,3]<-gsub("\\)","",tmp[,3])
  
  
  Clin_aa<-data.frame(gene_name=tmp2[,1],DNA_var=tmp2[,2],
                      aa_from=substr(tmp[,3],3,5),
                      aa_to=substr(tmp[,3],nchar(tmp[,3])-2,nchar(tmp[,3])),aa_pos=substr(tmp[,3],6,nchar(tmp[,3])-3),Classification=Clin_data$Classification)
  
  
  Clin_aa<-Clin_aa[which(Clin_aa$aa_pos!=""),]
  Clin_aa$aa_from<-convert_aa_code(Clin_aa$aa_from)
  Clin_aa$aa_to<-convert_aa_code(Clin_aa$aa_to)
  Clin_aa<-Clin_aa[!is.na(Clin_aa$aa_to),]
  
  Clin_aa<-inner_join(Clin_aa,prot_gene,by="gene_name")
  Clin_aa$id<-paste0(Clin_aa$aa_pos,Clin_aa$aa_from)
  
  length(unique(Clin_aa$id))
  
  target_CF<-topo_feature[which(topo_feature$Uniprot%in%unique(Clin_aa$uniprot_id)),]
  
  Clin_aa_CF<-inner_join(Clin_aa,target_CF,by="id")
  
  if(nrow(Clin_aa_CF)>0){
    
    Clin_aa_CF$Category<-Clin_aa_CF$Classification
    Clin_aa_CF$Classification[which(Clin_aa_CF$Classification%in%c("Likely Benign","Benign"))]<-"Likely Benign & Benign"
    Clin_aa_CF$Classification[which(Clin_aa_CF$Classification%in%c("Likely Pathogenic","Pathogenic"))]<-"Likely Pathogenic & Pathogenic"
    
    ClinGen_all<-rbind(ClinGen_all,cbind(Clin_aa_CF,Pannel=pannel))
    
    plot_CF<-Clin_aa_CF
    
    plot_CF$Classification<-factor(plot_CF$Classification,levels=c("Likely Pathogenic & Pathogenic","Likely Benign & Benign","Uncertain Significance"))
    
    if(nrow(plot_CF)>50 & length(which(plot_CF$Classification=="Likely Pathogenic & Pathogenic"))>4 & length(which(plot_CF$Classification=="Likely Benign & Benign"))>4 & length(which(plot_CF$Classification=="Uncertain Significance"))>4){
      
      plot_CF<-plot_CF[!duplicated(plot_CF[,c("Classification","uniprot_id","id")]),]
      plot_cf[[m]]<-plot_classification2(plot_CF,CF10QS,pannel,pdf_prefix="AF2_CF10QS")
      
      plot_ld[[m]]<-plot_classification2(plot_CF,LD15QS,pannel,pdf_prefix="AF2_LD15QS")
      m<-m+1
      
    } ## if(nrow(plot_CF)>20)
    
    
  } ### if(nrow(Clin_aa_CF)>0)
  
  
}  ##  for(i in 1:length(pannel_info))

table(ClinGen_all$Classification)
length(unique(ClinGen_all$Pannel))

comb_cf<-ggarrange(plotlist=plot_cf,ncol=5,nrow=4,common.legend = T)
ggsave(paste0(output_folder,"ClinGen_boxplot_AF2_CF10QS.pdf"),comb_cf,width=15,height=20)

comb_ld<-ggarrange(plotlist=plot_ld,ncol=5,nrow=4,common.legend = T)
ggsave(paste0(output_folder,"ClinGen_boxplot_AF2_LD15QS.pdf"),comb_ld,width=15,height=20)

write.table(ClinGen_all,
            paste0("ClinGen_AF2CF.csv"),
            append = FALSE, quote = T, sep = ",",row.names = F, col.names = T)



predict_clingene<-ClinGen_all[which(ClinGen_all$Classification%in%c("Likely Benign & Benign","Likely Pathogenic & Pathogenic")),]


########
color_set<-c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00"
  ,"#FFFF33","#A65628","#F781BF","#999999")
plot_feature<-c("CF10QS","LD15QS","CF10QS.Trend4bi")
auc_value<-rep(0,length(plot_feature))

pdf(file=paste0(output_folder,"ClinGene_AF2CF_ROC_all.pdf"), width=5, height=5)

for(i in 1:length(plot_feature)){
  
  roc_data<-roc(response = predict_clingene$Classification,predictor = predict_clingene[,plot_feature[i]]) 
  
  if(i==1){
    
    plot_roc<-plot(roc_data, col=color_set[i], main="AA topo estimate clinical effects",legacy.axes=T,lwd=0.6) 
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

##########
# protein_FI<-read.table(paste0("../4.Protein_category/AF2_WG_Proteins_Fragile.txt"),header=T,sep=",",quote = "",fill=TRUE)
# 
# ncut<-2
# protein_FI$FI_group <-paste0("FI_",ntile(protein_FI$LogG1G3_prop, ncut))
# quantile_labels<-paste0("FI_",c(1:ncut))
# 
# predict_clingene<-inner_join(predict_clingene,protein_FI,by=join_by("Uniprot"=="Uniprot"))
# 
# # plot_x<-"LogG1G3_prop"
# 
# table(predict_clingene$FI_group)
# 
# auc<-rep(0,ncut)
# for(i in 1:length(quantile_labels)){
#   
#   group_df<-predict_clingene[which(predict_clingene$FI_group==quantile_labels[i]),]
#   
#   roc_data<-roc(response = group_df$Classification,predictor = group_df[,"CF10QS"])
#   
#   auc[i]<-roc_data$auc
# }
# 
# 
# plot_data<-data.frame(FI_group=quantile_labels,AUC=auc)
# plot_data$FI_group<-factor(plot_data$FI_group,levels=quantile_labels)
# 
# p<-ggplot(plot_data, aes(x = FI_group)) +
#   geom_point(aes(x = FI_group,y = AUC), size = 3) +
#   geom_line(aes(x = FI_group,y = AUC, group = AUC), linewidth = 1.2) +
#   labs(title = "") +
#   theme_minimal() +
#   scale_x_discrete(limits = quantile_labels) +
#   theme(legend.position = "top",
#         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
# 
# p



########

pdf_file<-paste0(output_folder,"ClinGen_AF2_Trend_dens.pdf")


type_order<-c("Likely Pathogenic & Pathogenic","Likely Benign & Benign")
plot_dens_combn_cf(predict_clingene,type_var="Classification",type_order,topo_used,pdf_file)

