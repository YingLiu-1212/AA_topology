# This R script performs cross-species analysis of protein structural topology and half-life relationships using mouse proteome data. Building on the human protein analysis framework, this code extends the investigation to murine proteins to examine evolutionary conservation of structure-degradation relationships. The analysis processes experimental protein half-life measurements from mouse tissues and integrates them with AlphaFold2-predicted structural features specific to the mouse proteome. Key computational steps include: processing and normalizing mouse protein half-life data; handling multiple UniProt identifiers for individual proteins; classifying amino acid positions into four topological groups (G1-G4) using the same structural classification scheme applied to human proteins; computing proportional distributions of topological features; and performing correlation analysis between topological compositions and protein turnover rates. The script generates comparative visualizations including scatter plots showing relationships between half-life and individual topological group proportions, combined G1G3 and G2G4 composite metrics, and ROC curves evaluating the predictive power of topological features for classifying proteins with extreme half-lives. This cross-species analysis enables investigation of evolutionary conservation in structure-degradation relationships and provides validation of topological metrics as general predictors of protein stability across mammalian systems.

plot_pair_feature<-function(plot_data,feature1,feature2){
  
  plot_label<-paste0(plot_data$Cell[1],"(",nrow(plot_data)," Proteins)")
  
  # cor<-round(cor(plot_data[[feature1]],plot_data[[feature2]]),3)
  
  plot1<-ggplot(data=plot_data,aes(.data[[feature1]],.data[[feature2]])) +
    labs(title=paste0(plot_label))+ 
    #xlab("protein_rna_ratio")+ylab(feature)+
    # geom_point(alpha=0.2,size=1,color="#377EB8")+
    # geom_density_2d(colour="#E41A1C",linewidth=0.5)+
    geom_point(alpha=0.2,size=0.5,stroke=0.5)+
    geom_smooth(method="glm",colour="#FF2600",linewidth=0.5,fill="#377EB8")+
    stat_cor(r.digits = 2,p.digits =2,color="#FF2600")+
    theme_bw()
  
}

library(tidyr)
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

load(paste0(data_folder,"AF2Mouse_AA_TopoTrend_dedup.rdata"))

protein_hl<-read.table(paste0(data_folder,"Protein_stability&half_life/mouse_protein_half_lives_nature2011.csv"),header=T,sep=",",quote = "",fill=TRUE)

protein_hl<-separate_longer_delim(protein_hl,c(Uniprot.IDs), delim = ";")

# tmp<-str_split_fixed(protein_hl$Uniprot.IDs,";",n=2)
# protein_hl$Uniprot.IDs<-tmp[,1]

protein_hl<-protein_hl[which(protein_hl$Uniprot.IDs!=""),]
summary(protein_hl$half.life_avg)
protein_hl$HL<-protein_hl$half.life_avg
protein_hl$logHL<-round(log2(protein_hl$HL),3)
summary(protein_hl$logHL)
protein_in<-unique(protein_hl$Uniprot.IDs)

#######

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


# pdf_file<-paste0("Protein_HL_Mouse_Topo_group.pdf")
# plot_group<-topo_feature[which(topo_feature$group!="-"),]
# 
# plot_hdr_combn(plot_group,type_var="group",c("G1_G3","G2_G4"),topo_used,pdf_file)


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


protein_hl_topo<-inner_join(protein_hl,prot_topo,by=join_by("Uniprot.IDs"=="Uniprot"))

protein_hl_topo$LogG1_prop<-log2(0.01+protein_hl_topo$G1_prop)
protein_hl_topo$LogG2_prop<-log2(0.01+protein_hl_topo$G2_prop)
protein_hl_topo$LogG3_prop<-log2(0.01+protein_hl_topo$G3_prop)
protein_hl_topo$LogG4_prop<-log2(0.01+protein_hl_topo$G4_prop)


plot_data<-protein_hl_topo

plot_data$Cell<-"Mouse"

plot_var<-paste0("LogG",c(1:4),"_prop")

p<-list()
i<-1

for(var in plot_var){
  
  p[[i]]<-plot_pair_feature(plot_data,"logHL",var)
  i<-i+1
}


comb<-ggarrange(plotlist=p,ncol=4,nrow=1,common.legend = T)
ggsave(paste0(output_folder,"Protein_HL_Mouse_Topo_G1-G4_scatter.pdf"),comb,limitsize = FALSE,width=12,height=3)

plot_data<-protein_hl_topo

plot_data$LogG1G3_prop<-log2(0.01+plot_data$G1_prop+plot_data$G3_prop)
plot_data$LogG2G4_prop<-log(0.01+plot_data$G2_prop+plot_data$G4_prop)


write.table(plot_data,
            paste0("Protein_HL_Mouse_AF2_Topo.csv"),
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)


p<-plot_pair_feature(plot_data,"logHL","LogG1G3_prop")
ggsave(paste0(output_folder,"Protein_HL_Mouse_Topo_G1G3.pdf"),p,limitsize = FALSE,width=4,height=4)

p<-plot_pair_feature(plot_data,"logHL","LogG2G4_prop")
ggsave(paste0(output_folder,"Protein_HL_Mouse_Topo_G2G4.pdf"),p,limitsize = FALSE,width=4,height=4)



# p<-plot_pair_feature(plot_data,"logHL","G2G4_proption")
# 
# ggsave(paste0("Protein_HL_Mouse_Topo_G2G4.pdf"),p,limitsize = FALSE,width=4,height=4)



protein_hl_topo<-inner_join(protein_hl,prot_topo,by=join_by("Uniprot.IDs"=="Uniprot"))

########
sub_hl_topo<-protein_hl_topo[,c("Uniprot.IDs","HL","G1G3_proption","G2G4_proption")]

short_cut<-quantile(sub_hl_topo$HL,0.25)
long_cut<-quantile(sub_hl_topo$HL,0.75)

sub_hl_topo$Type<-"-"
sub_hl_topo$Type[which(sub_hl_topo$HL<short_cut)]<-"Short"
sub_hl_topo$Type[which(sub_hl_topo$HL>long_cut)]<-"Long"

sub_hl_topo<-sub_hl_topo[which(sub_hl_topo$Type!="-"),]

roc_data<-roc(response = sub_hl_topo$Type,predictor = sub_hl_topo$G1G3_proption)

roc_data$auc

pdf(file=paste0(output_folder,"Protein_HL_Mouse_Topo_G1G3_ROC.pdf"), width=5, height=5)

plot_roc<-plot(roc_data, main="G1G3_proption estimate protein half lives of mouse",legacy.axes=T,lwd=0.6) 
auc_value<-round(plot_roc$auc,2)
legend("bottomright",box.lwd = 0,
       legend=c(paste0("AUC=",round(auc_value,2))),cex = 0.7,lty=1)

dev.off()



