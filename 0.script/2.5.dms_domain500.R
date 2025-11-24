# This R script performs integrated analysis of deep mutational scanning (DMS) data with protein structural topology features for 500 human domains. The code processes missense mutation data per amino acid position, then integrates these functional measurements with structural metrics (CF10 and LD15) derived from AlphaFold2 predictions. Key functionalities include: data filtering and aggregation using dplyr, amino acid code conversion via bio3d, and comprehensive correlation analysis between functional scores and structural features. The script generates a 2×2 panel of scatter plots with smoothed trend lines and Spearman correlation statistics, visualizing relationships between fitness/abundance scores and structural parameters. All results are exported as both processed data tables and publication-quality PDF figures. The analysis enables systematic investigation of structure-function relationships in protein variants and supports reproducibility through modular data processing pipelines.

library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
options(scipen=0)
options(digits = 2)
options(stringsAsFactors = FALSE)


setwd("~/Other_project/AA_topology/reproducibility/b.Protein_expression")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

AA_fit_topology<-read.table(paste0(data_folder,"domain500/AA_fitness_AF2_topology.txt"),header=T,quote="",sep="\t",fill=TRUE)

prot_gene<-read.table(paste0(data_folder,"UP000005640_9606_prot_gene_seq.txt"),header=T,quote="",sep=",",fill=TRUE)
prot_gene$full_len<-nchar(prot_gene$AA_seq)

prot_gene<-prot_gene[,-3]

AA_fit_topology<-left_join(AA_fit_topology,prot_gene,by=join_by("uniprot_ID"=="uniprot_id"))



plot_pair_feature<-function(plot_data,feature1,feature2,color_fill){
  
  plot_label<-paste0(plot_data$uniprot_ID," (",plot_data$aa_len,"/",plot_data$full_len,")")

  plot1<-ggplot(data=plot_data,aes(.data[[feature1]],.data[[feature2]])) +
    # labs(title=paste0(plot_label,"; Cor.= ",cor))+ 
    labs(title=paste0(plot_label))+ 
    #xlab("protein_rna_ratio")+ylab(feature)+
    geom_point(alpha=0.3,size=1,stroke=0.5)+
    geom_smooth(colour="#FF2600",linewidth=0.5)+
    stat_cor(r.digits = 2,p.digits =2,color="#FF2600")+
    theme_minimal()
  
}


plot_df<-AA_fit_topology[which(AA_fit_topology$aa_len>=90),]
plot_df<-plot_df[order(plot_df$aa_len,decreasing = T),]

write.table(plot_df,
            paste0("dms_topo_domains.csv"),
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)

fit_prot<-unique(plot_df$uniprot_ID)


# plot_var<-"CF10"
# plot_var<-"LD15"

for(plot_var in c("CF10","LD15")){
  
  p<-list()
  
  for(j in 1:length(fit_prot)){
    
    plot_data<-plot_df[which(plot_df$uniprot_ID==fit_prot[j]),]
    
    p[[j]]<-plot_pair_feature(plot_data,plot_var,"fitness","wt_aa")
    
    
  }
  comb<-ggarrange(plotlist=p,ncol=4,nrow=4,common.legend = T)
  ggsave(paste0(output_folder,"DMS_domain500_AF2_",plot_var,".pdf"),comb,limitsize = FALSE,width=10,height=9)
  
  
}

