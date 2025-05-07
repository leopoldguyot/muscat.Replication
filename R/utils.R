library(dplyr)
process_results <- function(res, sim) {
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
        .$sim_mean.B == 0, "0", "log2(sim_mean.B/sim_mean.A)"))))
    
    resTable <- left_join(gi, resTable, by = c("gene", "cluster_id"), copy = TRUE) %>%
      {if ("logFC" %in% names(.))
        dplyr::rename(., est_lfc = logFC) else .} %>%
      dplyr::mutate(is_de = as.integer(!category %in% c("ee", "ep")))
  }
  resTable
}