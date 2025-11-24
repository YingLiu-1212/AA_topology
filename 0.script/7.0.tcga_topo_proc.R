# This R script processes TCGA (The Cancer Genome Atlas) pan-cancer mutation data and integrates it with AlphaFold2-predicted protein structural topology features. The code performs systematic mapping of cancer-associated mutations to their structural contexts by combining TCGA clinical and genomic data with protein structural information. Key processing steps include: loading and processing TCGA pan-cancer atlas data from 2018; mapping gene symbols to UniProt identifiers using protein-gene correspondence tables; integrating mutation data with structural topology features (CF10QS and other metrics) based on protein position and amino acid changes; filtering to retain only mutations with available structural topology information; and exporting the integrated dataset for subsequent cancer structural biology analyses. 

#library(pROC)
library(dplyr)
options(scipen=999)
options(stringsAsFactors = FALSE)

setwd("~/Other_project/AA_topology/reproducibility/g.Cancer")

# source("../0.script/_plot_fig.R")

data_folder<-"../a.resouce/"
output_folder<-"../o.output_figures/"


load(paste0(data_folder,"AF2_AA_TopoTrend_dedup.rdata"))
topo_feature$AA<-bio3d::aa321(topo_feature$AA)

prot_gene<-read.table(paste0(data_folder,"UP000005640_9606_prot_gene.txt"),quote="",header=T,sep="\t",fill=TRUE)

tcga_treatment<-read.table(paste0(data_folder,"all_tcga_pan_can_atlas_2018.csv"),quote="",header=T,sep=",",fill=TRUE)

tcga_treatment<-left_join(tcga_treatment,prot_gene,by=join_by("Hugo_Symbol"=="gene_name"))

tcga_treatment<-tcga_treatment[!is.na(tcga_treatment$uniprot_id),]

topo_feature$Pos<-as.character(topo_feature$Pos)

tcga_treatment_topo<-left_join(tcga_treatment,topo_feature,by=join_by("uniprot_id"=="Uniprot","Protein_position"=="Pos","aa_from"=="AA"))

tcga_treatment_topo<-tcga_treatment_topo[which(!is.na(tcga_treatment_topo$CF10QS)),]

write.table(tcga_treatment_topo,
            "tcga_atlas_2018_mut_topo.txt",
            append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)

save(tcga_treatment_topo,file="tcga_atlas_2018_mut_topo.rdata")

