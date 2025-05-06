if (!inherits(res$tbl, "error")) {
    # add metadata
    gi <- metadata(sim)$gene_info %>%
        dplyr::filter(gene %in% gs) %>%
        dplyr::mutate_at("cluster_id", as.character) %>%
        dplyr::select(-"logFC") %>%
        dplyr::mutate(., sim_lfc = eval(parse(text = ifelse(
            .$sim_mean.B == 0, "0", "log2(sim_mean.B/sim_mean.A)"))))

    res$tbl <- left_join(gi, res$tbl, by = c("gene", "cluster_id")) %>%
        {if ("logFC" %in% names(.))
            dplyr::rename(., est_lfc = logFC) else .} %>%
        dplyr::mutate(did = wcs$did, sid = wcs$sid, mid = wcs$mid,
                      i = wcs$i, j = wcs$j, g = wcs$g, c = wcs$c, k = wcs$k, s = wcs$s,
                      is_de = as.integer(!category %in% c("ee", "ep")))
}
