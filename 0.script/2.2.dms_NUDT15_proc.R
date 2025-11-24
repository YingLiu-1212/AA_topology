# This R script performs integrated analysis of deep mutational scanning (DMS) data with protein structural topology features for NUDT15. The code processes missense mutation data per amino acid position, then integrates these functional measurements with structural metrics (CF10 and LD15) derived from AlphaFold2 predictions. Key functionalities include: data filtering and aggregation using dplyr, amino acid code conversion via bio3d, and comprehensive correlation analysis between functional scores and structural features. The script generates a 2×2 panel of scatter plots with smoothed trend lines and Spearman correlation statistics, visualizing relationships between activity/abundance scores and structural parameters. All results are exported as both processed data tables and publication-quality PDF figures. The analysis enables systematic investigation of structure-function relationships in protein variants and supports reproducibility through modular data processing pipelines.

avg_min_n<-function(numbers,n){
  
  numbers<-numbers[!is.na(numbers)]
  if(length(numbers)>0){
    sorted_numbers <- sort(numbers,decreasing = F)
    smallest_n <- sorted_numbers[1:n]
    avg_of_smallest <- mean(smallest_n)
  }else{
    avg_of_smallest<-NA
  }
  
  
  avg_of_smallest
}

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(bio3d)
library(dplyr)
library(stringr)
options(scipen=0)
options(digits = 2)
options(stringsAsFactors = FALSE)

setwd("~/Other_project/AA_topology/reproducibility/b.Protein_expression")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

dms_gene<-read.table(paste0(data_folder,"NUDT15/DMS_NUDT15_abundance.txt"),quote="",header=T,sep="\t",fill=TRUE)

tmp<-str_split_fixed(dms_gene$Variant,">",n=2)

dms_gene$aa_from<-tmp[,1]
dms_gene$aa_to<-tmp[,2]

dms_gene<-dms_gene[which(dms_gene$Info=="missense"),]

summary(dms_gene$Average_abundance_score)
summary(dms_gene$Average_drug_sensitivity_score)
summary(dms_gene$Final.NUDT15.activity.Score)

dms_aa<-dms_gene %>%
  group_by(aa_from,Variant_position) %>%
  summarise(
    Average_abundance_score=mean(Average_abundance_score,na.rm=T),
    Average_drug_sensitivity_score=mean(Average_drug_sensitivity_score,na.rm=T),
    Final.NUDT15.activity.Score = mean(Final.NUDT15.activity.Score,na.rm=T)
    
  )

topo_trend<-read.table(paste0(data_folder,"NUDT15/AA_topology_AF-Q9NV35-F1-model_v4_chains.csv"),quote="",header=T,sep=",",fill=TRUE)
topo_trend$AA<-aa321(topo_trend$AA)

dms_topo_trend<-inner_join(dms_aa,topo_trend,by=join_by("Variant_position"=="Pos","aa_from"=="AA"))

data_df<-data.frame(dms_topo_trend)

write.table(data_df,
            paste0("dms_topo_NUDT15.csv"),
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)

gene_name<-"NUDT15"

pdf_prefix<-paste0(output_folder,"DMS_",gene_name)

p_list<-list()
i<-1
for(var_2 in c("CF10","LD15")){
  
  # cor<-round(cor(data_df[,var_1],data_df[,var_2],method="spearman",use="complete.obs"),2)
  # print(cor)
  
  for(var_1 in c("Average_abundance_score","Average_drug_sensitivity_score","Final.NUDT15.activity.Score")){
    
    p_list[[i]]<-ggplot(data=data_df,aes(.data[[var_1]],.data[[var_2]])) +
      labs(title=paste0(gene_name))+
      xlab(var_1)+
      ylab(var_2)+
      geom_point(alpha=0.3,size=1,stroke=0.5)+
      geom_smooth(colour="#FF2600",linewidth=0.5)+
      stat_cor(r.digits = 2,p.digits =2,color="#FF2600")+
      theme_minimal()
    
    i<-i+1
    
  }
  

}

p<-ggarrange(plotlist=p_list,ncol=3,nrow=2,common.legend = T)
ggsave(paste0(pdf_prefix,".pdf"),p,width = 9, height = 6)



