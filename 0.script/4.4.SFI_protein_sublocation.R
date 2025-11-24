# This R script analyzes the relationship between protein structural topology and subcellular localization using human proteome data from UniProt. The code integrates AlphaFold2-predicted structural features with comprehensive subcellular localization annotations to investigate how topological composition varies across different cellular compartments. Key analytical steps include: processing UniProt knowledgebase annotations to extract and categorize proteins into 22 distinct subcellular locations including nucleus, cytoplasm, organelles, and membranes; classifying amino acid positions into four topological groups (G1-G4) based on core-forming propensity and spatial trend metrics; calculating proportional distributions of structural features for each protein; and performing comparative analysis of topological compositions across subcellular locations. The script generates specialized visualizations including box plots that systematically compare G1G3 structural fragility index distributions across different cellular compartments, ordered by median values to reveal location-specific structural patterns. This analysis provides insights into how protein structural architecture correlates with cellular localization, potentially revealing compartment-specific structural constraints and adaptive features that may influence protein function, stability, and interactions within specific cellular environments.

library(ggdensity)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(ggpubr)
options(scipen=999)
options(digits = 9)
options(stringsAsFactors = FALSE)
setwd("~/Other_project/AA_topology/reproducibility/d.SFI")
source("../0.script/_plot_fig.R")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"

prot_anno<-read.table(paste0(data_folder,"uniprotkb_AND_model_organism_9606_AND_r_2025_04_29.tsv"),header=T,sep="\t",quote = "",fill=TRUE)

prot_anno<-prot_anno[,c(1,3,4,9)]

head(prot_anno$Subcellular.location..CC.)

cell_loc<-c("Nucleus","Chromosome","Nucleus, nucleoplasm","Nucleus speckle","Cytoplasm","Cytoplasm, cytoskeleton","Endoplasmic reticulum","Endoplasmic reticulum membrane","Early endosome","Late endosome","Golgi apparatus","Golgi apparatus membrane","Mitochondrion","Mitochondrion matrix","Mitochondrion inner membrane","Mitochondrion outer membrane","Lysosome","Lysosome membrane","Peroxisome","Cell membrane","Membrane","Secreted")

prot_anno$Loc<-"All"
for(loc in cell_loc){
  
  tmp<-which(grepl(loc,prot_anno$Subcellular.location..CC.))
  prot_anno$Loc[tmp]<-paste0(prot_anno$Loc[tmp],";",loc)
}

prot_anno_loc<-separate_longer_delim(prot_anno,c(Loc), delim = ";")

table(prot_anno_loc$Loc)

prot_anno_loc<-prot_anno_loc[,-4]

prot_anno_loc<-unique(prot_anno_loc)

write.table(prot_anno_loc,
            paste0("Human_protein_Uniprot_Loc.txt"),
            append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)




#######

load(paste0(data_folder,"AF2_AA_TopoTrend_dedup.rdata"))

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
    G1_prop = round(length(which(group=="G1"))/length(group),4),
    G2_prop = round(length(which(group=="G2"))/length(group),4),
    G3_prop = round(length(which(group=="G3"))/length(group),4),
    G4_prop = round(length(which(group=="G4"))/length(group),4),
    G1G3_proption = round(length(which(group%in%c("G1","G3")))/length(group),4),
    G2G4_proption = round(length(which(group%in%c("G2","G4")))/length(group),4),
    
    # G1_prop = length(which(group=="G1"))/length(group),
    # G2_prop = length(which(group=="G2"))/length(group),
    # G3_prop = length(which(group=="G3"))/length(group),
    # G4_prop = length(which(group=="G4"))/length(group),
    # G1G3_proption = length(which(group%in%c("G1","G3")))/length(group),
    # G2G4_proption = length(which(group%in%c("G2","G4")))/length(group),
    # length = length(group)
  )
length(unique(prot_topo$Uniprot))

prot_topo$LogG1G3_prop<-round(log2(0.01+prot_topo$G1G3_proption),4)
prot_topo$LogG2G4_prop<-round(log2(0.01+prot_topo$G2G4_proption),4)

# prot_topo$LogG1G3_prop<-log2(0.01+prot_topo$G1G3_proption)
# prot_topo$LogG2G4_prop<-log2(0.01+prot_topo$G2G4_proption)

write.table(prot_topo,
            "AF2_WG_Proteins_Fragile.csv",
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)



prot_loc<-read.table("Human_protein_Uniprot_Loc.txt",quote="\"",header=T,sep="\t",fill=TRUE)

prot_loc_topo<-inner_join(prot_loc,prot_topo,by=join_by("Entry"=="Uniprot"))


cell_loc<-c("All","Nucleus","Chromosome","Nucleus, nucleoplasm","Nucleus speckle","Cytoplasm","Cytoplasm, cytoskeleton","Endoplasmic reticulum","Endoplasmic reticulum membrane","Early endosome","Late endosome","Golgi apparatus","Golgi apparatus membrane","Mitochondrion","Mitochondrion matrix","Mitochondrion inner membrane","Mitochondrion outer membrane","Lysosome","Lysosome membrane","Peroxisome","Cell membrane","Membrane","Secreted")

prot_loc_topo$Loc<-factor(prot_loc_topo$Loc,levels=cell_loc)


p1<-ggplot(prot_loc_topo,aes(x=Loc,y=G1G3_proption))+
  #labs(title=stringr::str_wrap(use_feature,width=40),size=1)+
  #ylim(c(0,1.3))+
  #xlab("G1G3_proption")+
  #ylab(pdf_prefix)+
  geom_boxplot(aes(x = reorder(Loc, G1G3_proption, FUN = median), y = G1G3_proption),outlier.shape = NA)+
  theme_minimal() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=7),legend.position = "none")

ggsave(paste0(output_folder,"Uniprot_human_protein_Loc_SFI.pdf"),p1,width=10,height=4)



# 
# p2<-ggplot(prot_loc_topo,aes(x=factor(Loc),y= G2G4_proption))+
#   #labs(title=stringr::str_wrap(use_feature,width=40),size=1)+
#   #ylim(c(0,1.3))+
#   xlab("G2G4_proption")+
#   #ylab(pdf_prefix)+
#   geom_boxplot(aes(x = reorder(Loc, G2G4_proption, FUN = median), y = G2G4_proption),outlier.shape = NA)+
#   theme_minimal() +
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 15))+
#   theme(axis.text.x = element_text(angle = 60, hjust = 1,size=7),legend.position = "none")
# 
# ggsave(paste0("boxplot_Uniprot_human_Loc_G2G4.pdf"),p2,width=10,height=4)
# 
# 
