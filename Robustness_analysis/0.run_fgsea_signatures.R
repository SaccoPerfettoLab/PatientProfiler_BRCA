

run_fgsea_signatures <- function(gene_signatures, 
                                 transcriptomic_data, 
                                 NES_filter = TRUE, 
                                 pval_cutoff = 0.05) {
  # Aggregate once
  gene_mat <- aggregate(transcriptomic_data[,-1], by = list(gene_name = transcriptomic_data$gene_name), FUN = mean)
  gene_ids <- gene_mat$gene_name
  gene_mat <- gene_mat[,-1]
  rownames(gene_mat) <- gene_ids
  
  # Prepare fgsea results list
  enrichment_list <- vector("list", ncol(gene_mat))
  colnames_mat <- colnames(gene_mat)
  
  message('Running enrichment analysis')
  for (i in seq_along(enrichment_list)) {
    stats <- gene_mat[, i]
    names(stats) <- rownames(gene_mat)
    stats <- stats[!is.na(stats)]
    stats <- sort(stats, decreasing = TRUE)
    
    fgsea_res <- fgsea(pathways = gene_signatures, stats = stats, scoreType = "std") # "pos"
    if (nrow(fgsea_res) > 0) {
      fgsea_res$patient_id <- colnames_mat[i]
      enrichment_list[[i]] <- fgsea_res
    }
  }
  
  enrichment_df <- do.call(rbind, enrichment_list)
  enrichment_df <- enrichment_df[, c("pathway", "NES", "padj", "patient_id")]
  
  # Best signature per patient
  # enrichment_sign <- enrichment_df %>%
  #   filter(!NES_filter | NES > 0) %>%
  #   group_by(patient_id) %>%
  #   slice_min(order_by = padj, n = 1) %>%
  #   slice_max(order_by = NES, n = 1) %>%
  #   ungroup() %>%
  #   filter(padj <= pval_cutoff)
  
  return(enrichment_df)
}

