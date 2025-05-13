library(dplyr)
process_results_pb <- function(res, sim) {
    gs <- rownames(sim)
    resTable <- res$table$B
    resTable <- do.call("rbind", resTable)
    if (!inherits(res$tbl, "error")) {
        # add metadata
        gi <- metadata(sim)$gene_info %>%
            dplyr::filter(gene %in% gs) %>%
            dplyr::mutate_at("cluster_id", as.character) %>%
            dplyr::select(-"logFC") %>%
            dplyr::mutate(., sim_lfc = eval(parse(text = ifelse(
                .$sim_mean.B == 0, "0", "log2(sim_mean.B/sim_mean.A)"
            ))))

        resTable <- left_join(gi, resTable, by = c("gene", "cluster_id"), copy = TRUE) %>%
            {
                if ("logFC" %in% names(.)) {
                    dplyr::rename(., est_lfc = logFC)
                } else {
                    .
                }
            } %>%
            dplyr::mutate(is_de = as.integer(!category %in% c("ee", "ep")))
    }
    resTable
}

process_results_mm <- function(res, sim) {
  gs <- rownames(sim)
  
  # Combine list of dataframes into one big dataframe
  resTable <- bind_rows(res, .id = "cluster_id")
  
  # Extract gene metadata
  gi <- metadata(sim)$gene_info %>%
    filter(gene %in% gs) %>%
    mutate(cluster_id = as.character(cluster_id)) %>%
    select(-"logFC") %>%
    mutate(sim_lfc = ifelse(sim_mean.B == 0, 0, log2(sim_mean.B / sim_mean.A)))
  
  # Ensure cluster_id is character and join metadata
  resTable <- resTable %>%
    mutate(cluster_id = as.character(cluster_id)) %>%
    left_join(gi, by = c("gene", "cluster_id"))
  
  # Rename logFC or beta to est_lfc (if either exists)
  if ("beta" %in% names(resTable)) {
    resTable <- resTable %>% rename(est_lfc = beta)
  } else if ("logFC" %in% names(resTable)) {
    resTable <- resTable %>% rename(est_lfc = logFC)
  }
  
  # Final processing
  resTable <- resTable %>%
    rename(
      p_val = p_val,
      p_adj.loc = p_adj.loc,
      p_adj.glb = p_adj.glb
    ) %>%
    mutate(
      contrast = "B",
      is_de = as.integer(!category %in% c("ee", "ep"))
    )
  
  return(resTable)
}

process_results_scdd <- function(res, sim) {
  gs <- rownames(sim)
  
  # Prepare metadata
  gi <- metadata(sim)$gene_info %>%
    filter(gene %in% gs) %>%
    mutate(cluster_id = as.character(cluster_id)) %>%
    select(-"logFC") %>%
    mutate(sim_lfc = ifelse(sim_mean.B == 0, 0, log2(sim_mean.B / sim_mean.A)))
  
  # Join metadata with scDD results
  resTable <- res %>%
    mutate(cluster_id = as.character(cluster_id)) %>%
    left_join(gi, by = c("gene", "cluster_id")) %>%
    mutate(
      contrast = "scDD",
      est_lfc = NA_real_,
      lfcSE = NA_real_,
      stat = NA_real_,
      is_de = as.integer(!category %in% c("ee", "ep"))
    )
  
  return(resTable)
}
