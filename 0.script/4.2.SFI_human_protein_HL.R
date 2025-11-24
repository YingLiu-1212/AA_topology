# This R script investigates the relationship between protein structural topology and cellular half-life across different human cell types. The analysis integrates experimental protein half-life measurements from multiple cellular contexts with AlphaFold2-predicted structural features to examine how topological composition influences protein turnover rates. Key computational steps include: processing protein half-life data from multiple replicates and calculating log2-transformed values; classifying amino acid positions into four topological groups (G1-G4) based on core-forming propensity (CF10QS) and spatial trend metrics; computing proportional distributions of these topological groups for each protein; and performing comprehensive correlation analysis between topological compositions and half-life measurements across different cell types. The script generates multiple visualization outputs including scatter plots showing relationships between half-life and topological group proportions for individual cell types, combined multi-panel figures for comparative analysis, and ROC curves evaluating the predictive power of G1G3 topological features for classifying proteins with extreme half-lives. This systematic analysis provides insights into how structural architecture contributes to protein stability and degradation kinetics in cellular environments, revealing potential structural determinants of protein turnover regulation.

plot_pair_feature<-function(plot_data,feature1,feature2){
  
  plot_label<-paste0(plot_data$Cell[1],"(",nrow(plot_data)," Proteins)")
  
  # cor<-round(cor(plot_data[[feature1]],plot_data[[feature2]]),3)
  
  plot1<-ggplot(data=plot_data,aes(.data[[feature1]],.data[[feature2]])) +
    labs(title=paste0(plot_label))+ 
    #xlab("protein_rna_ratio")+ylab(feature)+
    # geom_point(alpha=0.2,size=1,color="#377EB8")+
    # geom_density_2d(colour="#E41A1C",linewidth=0.5)+
    geom_point(alpha=0.2,size=1,stroke=0.5)+
    geom_smooth(method="glm",colour="#FF2600",linewidth=0.5,fill="#377EB8")+
    stat_cor(r.digits = 2,p.digits =2,color="#FF2600")+
    theme_bw()
  
}


library(dplyr)
library(ggpubr)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggdensity)
library(RColorBrewer)
library(ggrepel)
options(scipen=0)
options(digits = 2)
options(stringsAsFactors = FALSE)
setwd("~/Other_project/AA_topology/reproducibility/d.SFI")
source("../0.script/_plot_fig.R")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

load(paste0(data_folder,"AF2_AA_TopoTrend_dedup.rdata"))

prot_gene<-read.table(paste0(data_folder,"UP000005640_9606_prot_gene.txt"),header=T,sep="\t",quote = "",fill=TRUE)

protein_hl<-read.table(paste0(data_folder,"Protein_stability&half_life/human_protein_half_lives_nc2018.csv"),header=T,sep=",",quote = "",fill=TRUE)

protein_hl$HL<-(protein_hl$Rep_1+protein_hl$Rep_2)/2

protein_hl$logHL<-round(log2(protein_hl$HL),3)

summary(protein_hl$logHL)

protein_hl<-inner_join(protein_hl,prot_gene,by="gene_name")

protein_in<-unique(protein_hl$uniprot_id)


topo_feature<-topo_feature[which(topo_feature$Uniprot%in%protein_in),]

length(unique(topo_feature$Uniprot))

cut<-0

G_1<-topo_feature$CF10QS+topo_feature$CF10QS.Trend4bi+1
G_2<-topo_feature$CF10QS+(1-topo_feature$CF10QS.Trend4bi)
G_3<-(1-topo_feature$CF10QS)+(0-topo_feature$CF10QS.Trend4bi)
G_4<-(1-topo_feature$CF10QS)+topo_feature$CF10QS.Trend4bi


topo_feature$group<-"-"
topo_feature$group[which(G_1<cut)]<-"G1"
topo_feature$group[which(G_2<cut)]<-"G2"
topo_feature$group[which(G_3<cut)]<-"G3"
topo_feature$group[which(G_4<cut)]<-"G4"
table(topo_feature$group)

prot_topo<-topo_feature %>%
  group_by(Uniprot) %>%
  summarise(
    G1_prop = length(which(group=="G1"))/length(group),
    G2_prop = length(which(group=="G2"))/length(group),
    G3_prop = length(which(group=="G3"))/length(group),
    G4_prop = length(which(group=="G4"))/length(group),
    G1G3_proption = length(which(group%in%c("G1","G3")))/length(group),
    G2G4_proption = length(which(group%in%c("G2","G4")))/length(group),
    length = length(group)
  )


protein_hl_topo<-inner_join(protein_hl,prot_topo,by=join_by("uniprot_id"=="Uniprot"))

protein_hl_topo$LogG1_prop<-log2(0.01+protein_hl_topo$G1_prop)
protein_hl_topo$LogG2_prop<-log2(0.01+protein_hl_topo$G2_prop)
protein_hl_topo$LogG3_prop<-log2(0.01+protein_hl_topo$G3_prop)
protein_hl_topo$LogG4_prop<-log2(0.01+protein_hl_topo$G4_prop)


####
plot_var<-paste0("LogG",c(1:4),"_prop")

cells<-unique(protein_hl_topo$Cell)

p<-list()
i<-1
for(j in 1:length(cells)){
  
  for(var in plot_var){
    
    plot_data<-protein_hl_topo[which(protein_hl_topo$Cell==cells[j]),]
    
    p[[i]]<-plot_pair_feature(plot_data,"logHL",var)
    i<-i+1
  }
  
  
}
comb<-ggarrange(plotlist=p,ncol=4,nrow=4,common.legend = T)
ggsave(paste0(output_folder,"Protein_HL_Human_AF2_Topo_G1-G4_scatter.pdf"),comb,limitsize = FALSE,width=12,height=12)


protein_hl_topo$LogG1G3_prop<-log2(0.01+protein_hl_topo$G1_prop+protein_hl_topo$G3_prop)
protein_hl_topo$LogG2G4_prop<-log2(0.01+protein_hl_topo$G2_prop+protein_hl_topo$G4_prop)


write.table(protein_hl_topo,
            paste0("Protein_HL_Human_AF2_Topo.csv"),
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)

# protein_hl_topo$G1G2_ratio<-(0.01+protein_hl_topo$G1_prop)/(0.01+protein_hl_topo$G2_prop)
# protein_hl_topo$G3G4_ratio<-(0.01+protein_hl_topo$G3_prop)/(0.01+protein_hl_topo$G4_prop)
# 
# cor(protein_hl_topo$logHL,protein_hl_topo$G1G2_ratio)
# cor(protein_hl_topo$logHL,protein_hl_topo$G3G4_ratio)
# 
# protein_hl_topo$Fragile_index<-protein_hl_topo$G1G2_ratio+protein_hl_topo$G3G4_ratio

# cor(protein_hl_topo$logHL,protein_hl_topo$Fragile_index)

####
cells<-unique(protein_hl_topo$Cell)

plots<-list()
k<-1
# candi_var<-c("LogG1G3_prop","G1G2_ratio","LogG2G4_prop","G3G4_ratio","Fragile_index")

candi_var<-c("LogG1G3_prop","LogG2G4_prop")

for(plot_var in candi_var){
  
  for(j in 1:length(cells)){
    
    plot_data<-protein_hl_topo[which(protein_hl_topo$Cell==cells[j]),]
    
    plots[[k]]<-plot_pair_feature(plot_data,"logHL",plot_var)
    
    k<-k+1
  }
  
}

comb<-ggarrange(plotlist=plots,ncol=4,nrow=length(candi_var),common.legend = T)
ggsave(paste0(output_folder,"Protein_HL_Human_AF2_Topo_comb.pdf"),comb,limitsize = FALSE,width=16,height=4*length(candi_var))


#####

sub_hl_topo<-protein_hl_topo[,c("uniprot_id","Cell","HL","G1G3_proption","G2G4_proption")]

short_cut<-quantile(sub_hl_topo$HL,0.25)
long_cut<-quantile(sub_hl_topo$HL,0.75)

########
pdf(file=paste0(output_folder,"Protein_HL_Human_AF2_Topo_G1G3_ROC.pdf"), width=5, height=5)

all_cells<-unique(sub_hl_topo$Cell)
auc_value<-rep(0,length(all_cells))
color_set<-brewer.pal(4, 'Dark2')
for(i in 1:length(all_cells)){
  
  print(all_cells[i])
  
  data_df<-sub_hl_topo[which(sub_hl_topo$Cell==all_cells[i]),]
  
  data_df$Type<-"-"
  data_df$Type[which(data_df$HL<short_cut)]<-"Short"
  data_df$Type[which(data_df$HL>long_cut)]<-"Long"
  
  data_df<-data_df[which(data_df$Type!="-"),]
  
  roc_data<-roc(response = data_df$Type,predictor = data_df$G1G3_proption) 
  
  if(i==1){
    
    plot_roc<-plot(roc_data, col=color_set[i], main="G1G3_proption estimate protein half lives of Human",legacy.axes=T,lwd=0.6) 
    auc_value[i]<-round(plot_roc$auc,2)
  }else{
    plot_roc<-plot(roc_data, add=TRUE, col=color_set[i],lwd=0.6)
    auc_value[i]<-round(plot_roc$auc,2)
  }
  
}
legend("bottomright",box.lwd = 0,
       legend=c(paste0(all_cells," AUC=",round(auc_value,2))),cex = 0.7, col=color_set,lty=1)

dev.off()  

#######



