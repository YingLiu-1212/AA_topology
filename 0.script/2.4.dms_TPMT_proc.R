# This R script performs integrated analysis of deep mutational scanning (DMS) data with protein structural topology features for TPMT. The code processes missense mutation data per amino acid position, then integrates these functional measurements with structural metrics (CF10 and LD15) derived from AlphaFold2 predictions. Key functionalities include: data filtering and aggregation using dplyr, amino acid code conversion via bio3d, and comprehensive correlation analysis between functional scores and structural features. The script generates a 2×2 panel of scatter plots with smoothed trend lines and Spearman correlation statistics, visualizing relationships between activity/abundance scores and structural parameters. All results are exported as both processed data tables and publication-quality PDF figures. The analysis enables systematic investigation of structure-function relationships in protein variants and supports reproducibility through modular data processing pipelines.

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

dms_gene<-read.table(paste0(data_folder,"PTEN&TPMT/DMS_TPMT_abundance.txt"),quote="",header=T,sep="\t",fill=TRUE)


dms_gene<-dms_gene[which(dms_gene$class=="missense"),]

summary(dms_gene$score)

dms_aa<-dms_gene %>%
  group_by(start,position) %>%
  summarise(
    # activity_score=avg_min_n(activity_score,3),
    #         abundance_score = avg_min_n(abundance_score,3)
    Average_abundance_score=mean(score,na.rm=T)
    
  )



# topo_trend<-read.table(paste0(data_folder,"PTEN&TPMT/AA_Trend_AF-P51580-F1-model_v4_chains.csv"),quote="",header=T,sep="\t",fill=TRUE)

topo_trend<-read.table(paste0(data_folder,"PTEN&TPMT/AA_topology_AF-P51580-F1-model_v4_chains.csv"),quote="",header=T,sep=",",fill=TRUE)
topo_trend$AA<-aa321(topo_trend$AA)

# cor(topo_trend$B,topo_trend$LD15)
# cor(topo_trend$B,topo_trend$LD15QS)


dms_topo_trend<-inner_join(dms_aa,topo_trend,by=join_by("position"=="Pos","start"=="AA"))

data_df<-data.frame(dms_topo_trend)


write.table(data_df,
            paste0("dms_topo_TPMT.csv"),
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)


cor(data_df$Average_abundance_score,data_df$CF10,use="na.or.complete")
cor(data_df$Average_abundance_score,data_df$LD15,use="na.or.complete")

gene_name<-"TPMT"
var_1<-"Average_abundance_score"


pdf_prefix<-paste0(output_folder,"DMS_",gene_name)

p_list<-list()
i<-1
for(var_2 in c("CF10","LD15")){
  
  # cor<-round(cor(data_df[,var_1],data_df[,var_2],method="spearman",use="complete.obs"),2)
  # print(cor)
  
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

p<-ggarrange(plotlist=p_list,ncol=2,nrow=1,common.legend = T)
ggsave(paste0(pdf_prefix,".pdf"),p,width = 6, height = 3)



