## _DN_function.R
library(dplyr)
library(stringr)

DN2_sep<-function(DN2_df,na_assign){
  
  # Processes dual-neighbor distance matrices by separating concatenated distance values 
  # for each amino acid type. This function takes a data frame containing combined distance 
  # strings (format: "distance1@neighbor1;distance2@neighbor2") for all 20 standard amino 
  # acids and parses them into individual numeric distance columns. For each amino acid, 
  # it extracts both distance values and retains the maximum distance between the two 
  # neighbors. The function handles missing values by optionally assigning a default 
  # distance of 100, making it suitable for computational analyses where complete matrices 
  # are required. Outputs a structured data frame with separate columns for each amino 
  # acid type's maximum neighbor distance.
  
  AA_order<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  
  DN2_sep<-DN2_df[,c("Uniprot","pdb_chain","id")]
  
  for(i in 1:length(AA_order)){
    #print(i)
    dist2_aa<-str_split_fixed(DN2_df[,AA_order[i]],";",n=2) 
    
    tmp1<-str_split_fixed(dist2_aa[,1],"@",n=2) 
    tmp2<-str_split_fixed(dist2_aa[,2],"@",n=2) 
    #dist<-as.numeric(tmp1[,1])+as.numeric(tmp2[,1])
    dist<-matrixStats::rowMaxs(cbind(as.numeric(tmp1[,1]),as.numeric(tmp2[,1])))
    if(na_assign==TRUE){
      dist[is.na(dist)]<-100
    }
    
    DN2_sep<-cbind(DN2_sep,dist)
  }
  
  colnames(DN2_sep)<-c(c("Uniprot","pdb_chain","id"),AA_order)
  
  DN2_sep
}


DN_sep2<-function(DN_df,na_assign,AA_order){
  
  # Processes single-neighbor distance matrices with customizable amino acid ordering.
  # This function parses distance strings (format: "distance@neighbor") for specified 
  # amino acid types and converts them into numeric distance values. Unlike DN_sep, 
  # it accepts a custom amino acid order parameter, providing flexibility for non-standard 
  # residue sets or specific analytical needs. Missing values can be automatically 
  # replaced with a default distance of 100 to ensure complete distance matrices. 
  # Returns a data frame with individual columns for each amino acid's distance 
  # measurements, preserving original protein identifiers and positional information.
  
  # AA_order<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  
  DN_sep<-DN_df[,c("Uniprot","chain_len","Pos","AA")]
  
  for(i in 1:length(AA_order)){
    #print(i)
    tmp<-str_split_fixed(DN_df[,AA_order[i]],"@",n=2) 
    dist<-as.numeric(tmp[,1])
    if(na_assign==TRUE){
      dist[is.na(dist)]<-100
    }
    
    DN_sep<-cbind(DN_sep,dist)
  }
  
  colnames(DN_sep)<-c(c("Uniprot","chain_len","Pos","AA"),AA_order)
  
  DN_sep
}

DN_sep<-function(DN_df,na_assign){
  
  # Processes standard single-neighbor distance matrices for all 20 amino acids.
  # This function parses concatenated distance-neighbor strings (format: "distance@neighbor") 
  # into separate numeric distance columns for each standard amino acid type. It 
  # systematically processes all 20 canonical amino acids in fixed order, extracting 
  # the distance component from each string. The function provides optional NA handling 
  # by replacing missing distances with a value of 100, which is useful for maintaining 
  # consistent matrix dimensions in subsequent analyses. Outputs a structured data 
  # frame with protein identifiers and individual distance columns for each residue type.
  
  AA_order<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  
  DN_sep<-DN_df[,c("Uniprot","pdb_chain","id")]
  
  for(i in 1:length(AA_order)){
    #print(i)
    tmp<-str_split_fixed(DN_df[,AA_order[i]],"@",n=2) 
    dist<-as.numeric(tmp[,1])
    if(na_assign==TRUE){
      dist[is.na(dist)]<-100
    }
    
    DN_sep<-cbind(DN_sep,dist)
  }
  
  colnames(DN_sep)<-c(c("Uniprot","pdb_chain","id"),AA_order)
  
  DN_sep
}

AA_ld_weight<-function(dist,chain_dn,AA_contact,AA_use){
  
  # Calculates amino acid contact weights based on local density (LD15RK) stratification.
  # This function computes the proportion of amino acid contacts occurring in high 
  # versus low local density regions. It divides the protein structure into two 
  # compartments based on the LD15RK metric (threshold = 0.5) and calculates the 
  # fraction of distances below the contact threshold in each compartment. The 
  # weighting provides insights into how amino acid contacts are distributed between 
  # densely packed interior regions and less dense exterior regions of the protein 
  # structure. Returns a vector with two values: outer_weight and inner_weight.
  
  subs<-which(chain_dn[,"LD15RK"]>=0.5)
  AA_inner_weight<-round(length(which(dist[subs]<AA_contact))/AA_use,4)
  
  subs<-which(chain_dn[,"LD15RK"]<0.5)
  
  AA_outer_weight<-round(length(which(dist[subs]<AA_contact))/AA_use,4)
  
  c(AA_outer_weight,AA_inner_weight)
  
}

AA_cf_weight<-function(dist,chain_dn,AA_contact,AA_use){
  
  # Computes amino acid contact weights stratified by core-forming (CF10RK) propensity.
  # This function analyzes the distribution of amino acid contacts between core-forming 
  # and surface regions defined by the CF10RK metric (threshold = 0.5). It calculates 
  # the proportion of distances below the specified contact threshold in each structural 
  # compartment, providing quantitative measures of how amino acids preferentially 
  # contact residues in the protein core versus surface regions. The weighting scheme 
  # helps characterize the structural preferences of different amino acid types 
  # in protein folding and stability. Returns a vector with outer_weight and inner_weight.
  
  subs<-which(chain_dn[,"CF10RK"]>=0.5)
  AA_outer_weight<-round(length(which(dist[subs]<AA_contact))/AA_use,4)
  
  subs<-which(chain_dn[,"CF10RK"]<0.5)
  
  AA_inner_weight<-round(length(which(dist[subs]<AA_contact))/AA_use,4)
  
  c(AA_outer_weight,AA_inner_weight)
  
}

get_AA_weight<-function(prot_dn_topo,AA_type,AA_contact,AA_use){
  
  # Generates comprehensive amino acid contact weight profiles across multiple structural contexts.
  # This function calculates five different contact weight metrics for a specific amino 
  # acid type within a protein structure: full structural context, core-forming outer 
  # and inner regions, and local density outer and inner regions. It integrates distance 
  # measurements with topological features to provide a multi-dimensional characterization 
  # of amino acid contact preferences. The function serves as a wrapper that combines 
  # full-structure analysis with compartment-specific calculations, returning a 
  # structured data frame with all weight metrics for systematic comparison across 
  # different structural environments.
  
  
  AA_dist<-prot_dn_topo[,AA_type]
  #print(AA_dist)
  
  AA_weight_full<-round(length(which(AA_dist<AA_contact))/AA_use,4)
  
  # AA_cf_weight(AA_dist,prot_dn_topo,AA_contact)
  # AA_ld_weight(AA_dist,prot_dn_topo,AA_contact)
  
  AA_weight<-data.frame(weight_type=paste0(AA_type,c("_weight_full","_CFouter_weight","_CFinner_weight","_LDouter_weight","_LDinner_weight")),
                        value=c(AA_weight_full,
               AA_cf_weight(AA_dist,prot_dn_topo,AA_contact,AA_use),
               AA_ld_weight(AA_dist,prot_dn_topo,AA_contact,AA_use)))
  
  # names(AA_weight)<-paste0(AA_type,c("_weight_full","_CFouter_weight","_CFinner_weight","_LDouter_weight","_LDinner_weight"))
  
  AA_weight
}

get_protein_AA_weight<-function(DN_topo_sep,protein_list,AA_contact,output_file){
  
  # Performs large-scale amino acid contact weight analysis across multiple proteins.
  # This function processes entire datasets of protein structures to compute comprehensive 
  # contact weight profiles for all 20 standard amino acids. It iterates through a 
  # provided protein list, calculates multiple contact weight metrics for each amino 
  # acid type in each protein, and appends the results to an output file. The function 
  # includes progress tracking and generates a tab-separated output file with complete 
  # weight profiles, enabling large-scale comparative analyses of amino acid structural 
  # preferences across diverse protein families. Supports batch processing of thousands 
  # of proteins with efficient file writing to handle large datasets.
  
  
  AA_order<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  
  # prot_all<-unique(protein_list$Uniprot)
  
  #WG_weight_all<-NULL
  
  for(p in 1:nrow(protein_list)){
    
    prot_id<-protein_list[p,"Uniprot"]
    AA_use<-protein_list[p,"use_AA"]
    
    if(p%%1000==0) print(p)
    
    prot_dn<-DN_topo_sep[which(DN_topo_sep$Uniprot==prot_id),]
    
    #AA_use<-nrow(prot_dn)
    
    AA_weight_all<-NULL
    
    for(i in 1:length(AA_order)){
      
      AA_weight<-get_AA_weight(prot_dn,AA_order[i],AA_contact,AA_use)
      
      AA_weight_all<-rbind(AA_weight_all,AA_weight)
      
    } ##for(i in 1:length(AA_order))
    
    
    AA_weight_all$Uniprot<-prot_id
    AA_weight_all$AA_use<-AA_use
    AA_weight_all$AA_target<-nrow(prot_dn)
    
    if(p==1){
      write.table(AA_weight_all,
                  output_file,
                  append = FALSE, quote = T, sep = "\t",row.names = F, col.names = T)
    }else{
      write.table(AA_weight_all,
                  output_file,
                  append = TRUE, quote = T, sep = "\t",row.names = F, col.names = F)
    }

  } ##for(p in 1:nrow(protein_list))
  
}







