suppressMessages({
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
    library(ExperimentHub)
    library(ggplot2)
    library(ggrastr)
})
source(file.path("LPSR", "utils.R"))

eh <- ExperimentHub()
sce <- eh[["EH3297"]]

pbsSum <- aggregate_data(sce,
                      by = c("cluster_id", "sample_id"),
                      fun = "sum")

pbsMean <- aggregate_data(sce,
                         by = c("cluster_id", "sample_id"),
                         fun = "mean")
plotMDSMean <- pbMDS(pbsMean) + scale_shape_manual(values = c("Vehicle" = 17, "LPS" = 8))
plotMDSSum <- pbMDS(pbsSum) + scale_shape_manual(values = c("Vehicle" = 17, "LPS" = 8))

ggsave(filename = "figs/MDS_LPS_Mean.pdf", plotMDSMean)
ggsave(filename = "figs/MDS_LPS_Sum.pdf", plotMDSSum)


ei <- metadata(sce)$experiment_info

ids <- c("cluster_id", "group_id", "sample_id")
cd_df <- data.frame(cell_id = seq_len(ncol(sce)),
                    colData(sce), do.call(cbind, reducedDims(sce)))
cd_df <- cd_df[sample(nrow(cd_df), nrow(cd_df)), ]
# for ea. cluster, sample equal group sizes
cs <- split(setDT(cd_df),
            by = c("cluster_id", "group_id"),
            flatten = FALSE, sorted = TRUE) %>%
    map_depth(2, "cell_id") %>%
    map_depth(1, function(u)
        lapply(u, sample, min(sapply(u, length))))
cs <- unname(unlist(cs))
df <- cd_df[sample(cs, length(cs)), ] # randomize rows

# plot UMAP colored by cluster, sample & group ID
for (id in ids) {
    plot <- ggplot(df, aes_string(x = "UMAP_1", y = "UMAP_2", col = id)) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        geom_point_rast(size = 0.1, alpha = 0.1) +
        theme_void() + theme(aspect.ratio = 1)
    ggsave(paste0("figs/umap_LPS_", id, ".pdf"), plot)
}
