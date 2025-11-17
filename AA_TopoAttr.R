
library(bio3d)
library(dplyr)
library(stringr)


options(scipen = 999, stringsAsFactors = FALSE)


avg_min_n <- function(numbers, n) {
  mean(sort(numbers, decreasing = FALSE)[1:n])
}

avg_max_n <- function(numbers, n) {
  mean(sort(numbers, decreasing = TRUE)[1:n])
}

local_density <- function(numbers, space) {
  pos <- which(numbers == 0)
  rm_range <- ceiling(space / 2)
  count_rm <- max(0, (pos - rm_range)):min(length(numbers), (pos + rm_range))
  count_in <- which(numbers < space)
  length(count_in[!count_in %in% count_rm])
}

add_nb_pos <- function(original_df, pos_var) {
  original_df %>%
    mutate(
      Pos_L1 = .data[[pos_var]] - 1,
      Pos_L2 = .data[[pos_var]] - 2,
      Pos_L3 = .data[[pos_var]] - 3,
      Pos_L4 = .data[[pos_var]] - 4,
      Pos_L5 = .data[[pos_var]] - 5,
      Pos_L6 = .data[[pos_var]] - 6,
      Pos_R1 = .data[[pos_var]] + 1,
      Pos_R2 = .data[[pos_var]] + 2,
      Pos_R3 = .data[[pos_var]] + 3,
      Pos_R4 = .data[[pos_var]] + 4,
      Pos_R5 = .data[[pos_var]] + 5,
      Pos_R6 = .data[[pos_var]] + 6
    )
}

add_nb_CF <- function(target_df, CF_df, LCC_var, output_var) {
  positions <- c("L1", "L2", "L3", "L4", "L5", "L6", "R1", "R2", "R3", "R4", "R5", "R6")
  
  for (pos in positions) {
    pos_col <- paste0("Pos_", pos)
    target_df <- left_join(
      target_df, 
      CF_df[, c("pdb_chain", "Pos", LCC_var)],
      by = setNames(c("pdb_chain", "Pos"), c("pdb_chain", pos_col))
    )
    

    new_col_name <- paste0(LCC_var, ".", pos)
    colnames(target_df)[ncol(target_df)] <- new_col_name
  }
  

  target_df <- left_join(
    target_df,
    CF_df[, unique(c("pdb_chain", "Pos", "AA", output_var, LCC_var))],
    by = c("pdb_chain", "Pos")
  )
  
  target_df
}

safe_cor <- function(x, y, ...) {
  tryCatch(round(cor(x, y, ...), 4), error = function(e) NA)
}

calc_topo_bi_avg <- function(topo_df, topo_var, n) {
  var_c <- paste0(".R", 1:n)
  var_n <- paste0(".L", 1:n)
  
  topo_c <- topo_df[, paste0(topo_var, c("", var_c))]
  topo_n <- topo_df[, paste0(topo_var, c("", var_n))]
  
  p_c <- apply(topo_c, 1, function(x) safe_cor(as.numeric(x), 1:ncol(topo_c), use = "complete.obs"))
  p_n <- apply(topo_n, 1, function(x) safe_cor(as.numeric(x), 1:ncol(topo_n), use = "complete.obs"))
  
  data.frame((p_c + p_n) / 2)
}

multi_topo_trend <- function(topo_var, nb_CF) {
  # Trend2bi <- calc_topo_bi_avg(nb_CF, topo_var, 2)
  # Trend3bi <- calc_topo_bi_avg(nb_CF, topo_var, 3)
  Trend4bi <- calc_topo_bi_avg(nb_CF, topo_var, 4)

  
  TopoTrends <- Trend4bi
  colnames(TopoTrends) <- paste0(topo_var, ".FPI")
  TopoTrends
}

QS_calc <- function(values) {
  qu <- rank(values, na.last = "keep") / length(values)
  round(qu * 2 - 1, 4)
}

convert_to_QS <- function(topo, trend_feature) {
  TrendQS_all <- NULL
  all_chains <- unique(topo$chain)
  
  for (chain_id in all_chains) {
    topo_chain <- topo[which(topo$chain == chain_id), ]
    for (feat in trend_feature) {
      topo_chain[, feat] <- QS_calc(topo_chain[, feat])
    }
    TrendQS_all <- rbind(TrendQS_all, topo_chain)
  }
  TrendQS_all
}

# 主函数
main <- function() {
  # 获取命令行参数
  args <- commandArgs(trailingOnly = TRUE)
  pdb_file <- args[1]
  output_folder <- args[2]
  
  # 第一部分：计算拓扑特征
  pdb_id <- substr(tools::file_path_sans_ext(basename(pdb_file)), 1, nchar(basename(pdb_file)) - 4)
  output_label <- paste0(pdb_id)
  pdb <- read.pdb(pdb_file)
  chains <- unique(pdb$atom$chain)
  
  Topology_chains <- NULL
  
  for (chain_id in chains) {
    ca.inds <- atom.select(pdb, "calpha", insert = "", chain = chain_id)
    
    if (length(ca.inds$atom) <= 50) next
    
    distance <- dm(pdb, inds = ca.inds, mask.lower = FALSE)
    diag(distance) <- 0
    
    select_atom <- pdb$atom[ca.inds$atom, ]
    resid <- select_atom$resid
    resno <- select_atom$resno
    resid_no <- paste0(resno, resid)
    
    colnames(distance) <- rownames(distance) <- resid_no
    distance <- round(distance, 4)
    
    coord <- select_atom
    coord$id <- paste0(coord$resno, coord$resid)
    
    if (!identical(rownames(distance), coord$id) || any(duplicated(rownames(distance)))) next
    
    if (min(coord$resno) <= 0) {
      prot_in <- which(coord$resno > 0)
      distance <- distance[prot_in, prot_in]
      coord <- coord[prot_in, ]
    }
    

    Topology_df <- data.frame(
      pdb_chain = paste0(pdb_id, "_", chain_id),
      chain_len = nrow(distance),
      id = rownames(distance),
      Pos = coord$resno,
      AA = coord$resid
    )
    

    for (nums in c(8, 10, 15)) {
      cf <- round(apply(distance, 1, avg_min_n, n = nums + 1), 4)
      Topology_df[[paste0("CF", nums)]] <- cf
      qu <- rank(Topology_df[[paste0("CF", nums)]]) / nrow(Topology_df)
      Topology_df[[paste0("CF", nums, "QS")]] <- round(qu * 2 - 1, 4)
      Topology_df[[paste0("CF", nums)]] <- round(scale(log10(Topology_df[[paste0("CF", nums)]]))[, 1], 4)
    }
    

    for (space in c(10, 15, 20)) {
      ld <- round(apply(distance, 1, local_density, space = space), 4)
      Topology_df[[paste0("LD", space)]] <- ld
      qu <- rank(Topology_df[[paste0("LD", space)]]) / nrow(Topology_df)
      Topology_df[[paste0("LD", space, "QS")]] <- round(qu * 2 - 1, 4)
      Topology_df[[paste0("LD", space)]] <- round(scale(log10(Topology_df[[paste0("LD", space)]] + 1))[, 1], 4)
    }
    
    Topology_chains <- rbind(Topology_chains, Topology_df)
  }
  
  if (is.null(Topology_chains)) {
    stop("No valid chains found in the PDB file.")
  }
  

  calc_var <- "CF10QS"
  output_var <- c("CF10", "CF10QS", "LD15", "LD15QS")
  
  topo_df <- Topology_chains
  topo_df <- topo_df[!duplicated(topo_df[, c("pdb_chain", "Pos")]), ]
  

  topo_df$pdb_id <- substr(topo_df$pdb_chain, 1, nchar(topo_df$pdb_chain) - 2)
  topo_df$chain <- substr(topo_df$pdb_chain, nchar(topo_df$pdb_chain), nchar(topo_df$pdb_chain))
  
  topo_info <- unique(topo_df[, c("pdb_chain", "pdb_id", "chain", "chain_len", "id", "Pos")])
  aa_nb <- add_nb_pos(topo_info, "Pos")
  aa_nb_CF <- add_nb_CF(aa_nb, topo_df, calc_var, output_var)
  

  CF_Trend <- multi_topo_trend(calc_var, aa_nb_CF)
  trend_feature <- colnames(CF_Trend)
  
  topo_comb <- cbind(aa_nb_CF[, c("pdb_chain", "pdb_id", "chain", "chain_len", "id", "Pos", "AA", output_var)], CF_Trend)
  topo_Trend <- convert_to_QS(topo_comb, trend_feature)
  

  output_file <- paste0(output_folder, "AA_TopoAttr_", output_label, ".txt")
  write.table(topo_Trend, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  message("Analysis complete. Results saved to: ", output_file)
}


main()