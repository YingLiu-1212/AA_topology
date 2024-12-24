### AA_topoloty.R
#install.packages("bio3d", dependencies=TRUE)
avg_min_n<-function(numbers,n){
  
  sorted_numbers <- sort(numbers,decreasing = F)
  smallest_n <- sorted_numbers[1:n]
  avg_of_smallest <- mean(smallest_n)
  avg_of_smallest
}

avg_max_n<-function(numbers,n){
  
  sorted_numbers <- sort(numbers,decreasing = T)
  largest_n <- sorted_numbers[1:n]
  avg_of_largest <- mean(largest_n)
  avg_of_largest
}

local_density<-function(numbers,space){
  pos<-which(numbers==0)
  rm_range<-ceiling(space/2)
  count_rm<-c(max(0,(pos-rm_range)):min(length(numbers),(pos+rm_range)))
  
  count_in<-which(numbers<space)
  count_in<-count_in[!count_in%in%count_rm]
  
  length(count_in)
}

library(bio3d)
library(dplyr)
options(scipen=999)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

pdb_file<-args[1]
chain_id<-args[2]

print(pdb_file)
print(chain_id)

# pdb_file<-"D:/weilab_other/CRC/Code/5w0u.pdb"
# chain_id<-"A"

pdb_id<-substr(tools::file_path_sans_ext(basename(pdb_file)), 1, nchar(basename(pdb_file))-4)

output_label<-paste0(pdb_id,"_",chain_id)
pdb <- read.pdb(pdb_file)
chains <- pdb$chain

ca.inds <- atom.select(pdb, "calpha",chain=chain_id)

if (length(ca.inds$atom) <= 50) {
  stop("The number of Cα atoms in the specified chain is less than 50", call. = FALSE)
}

distance<-dm(pdb, inds = ca.inds,mask.lower=FALSE)
diag(distance)<-0

select_atom<-pdb$atom[ca.inds$atom, ]

resid<-select_atom$resid
resno<-select_atom$resno
resid_no<-paste0(resno,aa321(resid))

colnames(distance)<-rownames(distance)<-resid_no

distance<-round(distance,4)

select_atom$resid<-aa321(select_atom$resid)

inter_contact<-distance
coord<-select_atom
coord$id<-paste0(coord$resno,coord$resid)

#identical(rownames(inter_contact), colnames(inter_contact))

filter<-identical(rownames(inter_contact),coord$id) & length(which(duplicated(rownames(inter_contact))))==0

if (!filter) {
  stop("Filtering criteria not met. Check the input PDB file and chain", call. = FALSE)
}

if(min(coord$resno)<=0){
  prot_in<-which(coord$resno>0)
  inter_contact<-inter_contact[prot_in,prot_in]
  coord<-coord[prot_in,]
}


# Prepare Topology data frame
Topology_df <- data.frame(
  pdb_chain = paste0(pdb_id, "_", chain_id),
  chain_len = nrow(inter_contact),
  id = rownames(inter_contact),
  Pos = coord$resno,
  AA = coord$resid
)

# Calculate Contact flexibility

for(nums in c(5,8,10)){

  cf<-round(unlist(apply(inter_contact,1,avg_min_n,n=nums+1)),4)
  
  Topology_df[[paste0("CF", nums)]] <- cf
  Topology_df[[paste0("CF", nums, "RK")]] <- round(rank(Topology_df[[paste0("CF", nums)]]) / nrow(Topology_df), 4)
  Topology_df[[paste0("CF", nums)]] <- round(scale(log10(Topology_df[[paste0("CF", nums)]]))[, 1], 4)
  
}

# Calculate Local Density (LD)
for (space in c(5, 8, 10, 12, 15, 20)) {
  ld <- round(unlist(apply(inter_contact, 1, local_density, space = space)), 4)
  Topology_df[[paste0("LD", space)]] <- ld
  Topology_df[[paste0("LD", space, "RK")]] <- round(rank(Topology_df[[paste0("LD", space)]]) / nrow(Topology_df), 4)
  Topology_df[[paste0("LD", space)]] <- round(scale(log10(Topology_df[[paste0("LD", space)]] + 1))[, 1], 4)
}

###########


write.table(Topology_df,
            paste0("AA_topology_",output_label,".csv"),
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)










