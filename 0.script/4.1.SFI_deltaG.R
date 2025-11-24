# This R script analyzes the relationship between protein structural topology features and thermodynamic stability (ΔG) using experimental data from deep mutational scanning studies. The code processes protein stability measurements from the Nature 2023 dataset, integrating them with AlphaFold2-derived structural topology metrics to investigate how different topological environments correlate with protein folding stability. Key analytical steps include: filtering and processing mutation data to extract wild-type stability measurements; calculating four distinct topological groups (G1-G4) based on combinations of core-forming propensity (CF10QS) and spatial trend features (CF10QS.Trend4bi); computing proportional distributions of these topological groups across protein structures; and performing correlation analysis between topological compositions and thermodynamic stability. The script generates multiple visualization outputs including scatter plots showing relationships between ΔG and topological group proportions, and ROC curves evaluating the predictive power of topological features for classifying protein stability extremes. This analysis provides insights into how structural architecture influences protein stability and enables quantitative assessment of topological metrics as predictors of folding energetics.

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

options(stringsAsFactors=FALSE)
library(bio3d)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(pROC)
options(scipen=0)
options(digits = 2)
setwd("~/Other_project/AA_topology/reproducibility/d.SFI")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

deltaG<-read.table(paste0(data_folder,"Protein_stability&half_life/DeltaG_pdb_nature2023.txt"),quote="\"",header=T,sep="\t",fill=TRUE)

deltaG<-deltaG[which(!grepl("ins",deltaG$mut_type) & !grepl("del",deltaG$mut_type) & !grepl(":",deltaG$mut_type) & !grepl("_",deltaG$WT_name)),]

deltaG$Ref<-substr(deltaG$mut_type,1,1)
deltaG$Pos<-substr(deltaG$mut_type,2,nchar(deltaG$mut_type)-1)
deltaG$Alt<-substr(deltaG$mut_type,nchar(deltaG$mut_type),nchar(deltaG$mut_type))

deltaG$Ref[which(deltaG$mut_type=="wt")]<-"-"
deltaG$Pos[which(deltaG$mut_type=="wt")]<-0
deltaG$Alt[which(deltaG$mut_type=="wt")]<-"-"


deltaG$Pos<-as.numeric(deltaG$Pos)

deltaG$id<-paste0(deltaG$Pos,deltaG$Ref)
deltaG$id[which(deltaG$mut_type=="wt")]<-"-"
deltaG$Pos[which(deltaG$mut_type=="wt")]<-NA

AA_deltaG<-deltaG %>%
  group_by(pdb,id) %>%
  reframe(
    Pos=as.numeric(unique(Pos)),
    AA=unique(Ref),
    WT_name=unique(WT_name),
    WT_cluster=unique(WT_cluster),
    deltaG = mean(deltaG)
  )

AA_deltaG_WT<-AA_deltaG[which(AA_deltaG$id=="-"),]
AA_deltaG_MUT<-AA_deltaG[which(AA_deltaG$id!="-"),]

topo_feature<-read.table(paste0(data_folder,"DeltaG_Topology.txt"),quote="\"",header=T,sep="\t",fill=TRUE)



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
  group_by(pdb_id) %>%
  summarise(
    G1_prop = length(which(group=="G1"))/length(group),
    G2_prop = length(which(group=="G2"))/length(group),
    G3_prop = length(which(group=="G3"))/length(group),
    G4_prop = length(which(group=="G4"))/length(group),
    length = length(group)
  )
# 
length(unique(prot_topo$pdb_id))

protein_deltaG_topo<-inner_join(AA_deltaG_WT,prot_topo,by=join_by("pdb"=="pdb_id"))

plot_var<-paste0("G",c(1:4),"_prop")

p<-list()
i<-1
for(var in plot_var){
  plot_data<-protein_deltaG_topo
  p[[i]]<-plot_pair_feature(plot_data,"deltaG",var)
  i<-i+1
}
comb<-ggarrange(plotlist=p,ncol=4,nrow=1,common.legend = T)
ggsave(paste0(output_folder,"Protein_deltaG_Human_Topo_G1-G4_scatter.pdf"),comb,limitsize = FALSE,width=12,height=3)


######


protein_deltaG_topo$G1G3_prop<-plot_data$G1_prop+plot_data$G3_prop
protein_deltaG_topo$G2G4_prop<-plot_data$G2_prop+plot_data$G4_prop


write.table(protein_deltaG_topo,
            paste0("Protein_deltaG_Human_AF2_Topo.csv"),
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)


plot_data<-protein_deltaG_topo

p<-plot_pair_feature(plot_data,"deltaG","G1G3_prop")
ggsave(paste0(output_folder,"Protein_deltaG_Topo_SFI.pdf"),p,limitsize = FALSE,width=4,height=4)

p<-plot_pair_feature(plot_data,"deltaG","G2G4_prop")
ggsave(paste0(output_folder,"Protein_deltaG_Topo_SSI.pdf"),p,limitsize = FALSE,width=4,height=4)

########
sub_deltaG_topo<-protein_deltaG_topo[,c("pdb","deltaG","G1G3_prop","G2G4_prop")]

short_cut<-quantile(sub_deltaG_topo$deltaG,0.25)
long_cut<-quantile(sub_deltaG_topo$deltaG,0.75)

sub_deltaG_topo$Type<-"-"
sub_deltaG_topo$Type[which(sub_deltaG_topo$deltaG<short_cut)]<-"Short"
sub_deltaG_topo$Type[which(sub_deltaG_topo$deltaG>long_cut)]<-"Long"

sub_deltaG_topo<-sub_deltaG_topo[which(sub_deltaG_topo$Type!="-"),]

roc_data<-roc(response = sub_deltaG_topo$Type,predictor = sub_deltaG_topo$G1G3_prop)

pdf(file=paste0(output_folder,"Protein_deltaG_Topo_SFI_ROC.pdf"), width=5, height=5)

plot_roc<-plot(roc_data, main="G1G3_proportion estimate protein stability",legacy.axes=T,lwd=0.6)
auc_value<-round(plot_roc$auc,2)
legend("bottomright",box.lwd = 0,
       legend=c(paste0("AUC=",round(auc_value,2))),cex = 0.7,lty=1)

dev.off()



