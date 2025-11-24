# This R script performs integrated analysis of phosphorylation sites in mouse tissues, combining phosphoproteomics data with AlphaFold2-predicted structural and charge distribution features. The code systematically investigates how structural context influences phosphorylation site properties, including phosphorylation abundance, turnover rates (T50), and their relationship with protein-level stability metrics. Key analytical steps include: processing phosphorylation site data across multiple mouse tissues and cell types; integrating structural topology metrics (CF10QS) with basic charge distribution features (BasicQS) to create composite phosphorylation accessibility scores (phos_acc); generating tissue-specific analyses comparing brain tissues with other organs; creating comprehensive visualizations including interval bar plots showing relationships between phosphorylation turnover rates and structural features across serine, threonine, and tyrosine residues; comparative box plots analyzing phosphorylation accessibility across different turnover rate categories; and scatter plots examining correlations between phosphorylation site turnover and protein-level stability across different cell types. This multi-dimensional analysis provides insights into how structural and electrostatic environments influence phosphorylation site dynamics, stability, and regulatory functions across different biological contexts in mouse models.

library(tidyr)
library(dplyr)
library(ggpubr)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggdensity)
library(RColorBrewer)
library(ggrepel)
options(scipen=9)
# options(digits = 3)
options(stringsAsFactors = FALSE)
setwd("~/Other_project/AA_topology/reproducibility/f.PTM&AD")

source("../0.script/_plot_fig.R")
source("../0.script/_DN_function.R")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

load(paste0(data_folder,"AF2_Mouse_AA_TopoTrend_dedup.rdata"))
load(paste0(data_folder,"AF2_Mouse_charged_Dist.rdata"))

phos_info<-read.table(paste0(data_folder,"PTM_and_AD/final_phosphorylation_abu_t50.txt"),header=T,sep="\t",quote = "",fill=TRUE)

# cor(log2(0.1+phos_info$phos_abundance),log2(0.1+phos_info$phos_T50),use="na.or.complete")

phos_abu_t50<-separate_longer_delim(phos_info,c(PTMsite), delim = ",")
phos_abu_t50$PTMsite<-gsub("\\(","",phos_abu_t50$PTMsite)
phos_abu_t50$PTMsite<-gsub("\\)","",phos_abu_t50$PTMsite)
phos_abu_t50<-unique(phos_abu_t50)

cell_order<-unique(phos_abu_t50$Cell)

#####
protein_info<-read.table(paste0(data_folder,"PTM_and_AD/final_protein_abu_t50.txt"),header=T,sep="\t",quote = "",fill=TRUE)

phos_abu_t50_protein<-left_join(phos_abu_t50,protein_info,by=c("Proteins","Cell"))

phos_abu_t50_protein$logProteinT50<-log2(0.5+phos_abu_t50_protein$Protein_T50)
phos_abu_t50_protein$logProteinT50[which(phos_abu_t50_protein$logProteinT50>10)]<-10


phos_abu_t50_protein$logPhosT50<-log2(0.5+phos_abu_t50_protein$phos_T50)
phos_abu_t50_protein$logPhosT50[which(phos_abu_t50_protein$logPhosT50>10)]<-10






######

protein_in<-unique(phos_abu_t50_protein$Proteins)

topo_feature<-topo_feature[which(topo_feature$Uniprot%in%protein_in),]

topo_feature$AA<-bio3d::aa321(topo_feature$AA)
topo_feature$id<-paste0(topo_feature$AA,topo_feature$Pos)

phos_topo<-inner_join(phos_abu_t50_protein,topo_feature,by=join_by("Proteins"=="Uniprot","PTMsite"=="id"))

phos_topo<-phos_topo[which(phos_topo$AA%in%c("S","T","Y")),]

######
charged_aa_dist<-charged_aa_dist[which(charged_aa_dist$Uniprot%in%protein_in),]


phos_protein_comb<-inner_join(phos_topo,charged_aa_dist[,c("Uniprot","id","BasicQS","AcidicQS")],by=join_by("Proteins"=="Uniprot","PTMsite"=="id"))


phos_protein_comb$phos_acc<-phos_protein_comb$CF10QS+phos_protein_comb$BasicQS

phos_protein_comb$phos_logAb<-log2(0.1+phos_protein_comb$phos_abundance)

phos_protein_comb$tissue<-"Brain"
phos_protein_comb$tissue[which(phos_protein_comb$Cell%in%c("Liver","Kidney","Spleen","Lung","Heart","Gut","Plasma"))]<-"Other"

write.table(phos_protein_comb,
            paste0("Mouse_phos_topo.csv"),
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)



#######
plot_phos<-phos_protein_comb

# plot_x<-"logProteinT50"
plot_x<-"logPhosT50"


for(tissue in c("Brain","Other")){
  
  tri_colsets<-matrix(rep(c("#377EB8","#B3CDE3","#4DAF4A","#CCEBC5","#E41A1C","#FBB4AE"),3),byrow = T,ncol=2)
  p<-1
  plot_list<-list()
  for(aa in c("S","T","Y")){
    
    for(feature in c("CF10QS","BasicQS","phos_acc")){
      
      plot_df<-plot_phos[which(plot_phos$AA==aa & plot_phos$tissue==tissue),]
      
      plot_list[[p]]<-plot_bar_interval_10(plot_df,plot_x,feature,paste0(feature," of ",aa),tri_colsets[p,])
      p<-p+1
    }
  }
  
  
  pdf_file<-paste0(output_folder,"Mouse_Phos_AA_",plot_x,"_",tissue,"_bar_intervalue.pdf")
  fig<-ggarrange(plotlist=plot_list,ncol=3,nrow=3)
  ggsave(pdf_file,fig,width=10,height=8,limitsize = FALSE)
}
######

plot_phos<-phos_protein_comb
for(tissue in c("Brain","Other")){
  p<-1
  plot_list<-list()
  for(plot_x in c("phos_T50","Protein_T50")){
    for(aa in c("S","T","Y")){
      
      plot_df<-plot_phos[which(plot_phos$AA==aa & plot_phos$tissue==tissue),]
      
      plot_list[[p]]<-qu_boxplot(plot_df,0.1,plot_x,"phos_acc",paste0(plot_x," of ",aa))
      p<-p+1
      
    }
  }
  
  pdf_file<-paste0(output_folder,"Mouse_Phos_AA_",tissue,"_box.pdf")
  p<-ggarrange(plotlist=plot_list,ncol=6,nrow=1,common.legend = T)
  
  ggsave(pdf_file,p,width=10,height=4,limitsize = FALSE)
  
}

########
scatter_pair_feature<-function(plot_data,feature1,feature2,title_prefix){
  
  plot_data<-plot_data[which(!is.na(plot_data[,feature1])&!is.na(plot_data[,feature2])),]
  
  plot_label<-paste0(title_prefix," ",nrow(plot_data)," sites")
  
  # cor<-round(cor(plot_data[[feature1]],plot_data[[feature2]]),3)
  
  plot1<-ggplot(data=plot_data,aes(.data[[feature1]],.data[[feature2]])) +
    labs(title=paste0(plot_label))+ 
    #xlab("protein_rna_ratio")+ylab(feature)+
    # geom_point(alpha=0.2,size=1,color="#377EB8")+
    # geom_density_2d(colour="#E41A1C",linewidth=0.5)+
    geom_point(alpha=0.2,size=0.5,stroke=0.5)+
    geom_smooth(method="glm",colour="#FF2600",linewidth=0.5,fill="#377EB8")+
    stat_cor(r.digits = 2,p.digits =2,color="#FF2600")+
    theme_bw()+
    theme(plot.title = element_text(size = 10))
  
  
  plot1
  
}

plot_phos<-phos_protein_comb

p<-1
plot_list<-list()
for(cell in cell_order){
  
  plot_df<-plot_phos[which(plot_phos$Cell==cell),]
  
  plot_list[[p]]<-scatter_pair_feature(plot_df,"logPhosT50","logProteinT50",cell)
  p<-p+1
}

pdf_file<-paste0(output_folder,"Mouse_Phos_scatter_phos_protein_HL.pdf")
p<-ggarrange(plotlist=plot_list,ncol=3,nrow=5)

ggsave(pdf_file,p,width=7.5,height=13,limitsize = FALSE)
