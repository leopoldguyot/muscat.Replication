suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(scDD)
    library(SingleCellExperiment)
    library(tidyverse)
})

apply_edgeR <- function(sce) {
    pbDS(pb = sce, method = "edgeR", filter = "none", verbose = FALSE)
}

apply_limma_trend <- function(sce) {
    pbDS(pb = sce, method = "limma-trend", filter = "none", verbose = FALSE)
}

apply_limma_voom <- function(sce) {
    pbDS(pb = sce, method = "limma-voom", filter = "none", verbose = FALSE)
}

apply_DESeq2 <- function(sce) {
    pbDS(pb = sce, method = "DESeq2", filter = "none", verbose = FALSE)
}

apply_MM_vst <- function(sce) {
    counts(sce) <- assay(sce)
    mmDS(x = sce, method = "vst", verbose = FALSE)
}

apply_MM_nbinom <- function(sce) {
    counts(sce) <- assay(sce)
    mmDS(x = sce, method = "nbinom", verbose = FALSE)
}

apply_MM_dream2 <- function(sce) {
    counts(sce) <- assay(sce)
    mmDS(x = sce, method = "dream2", verbose = FALSE)
}

apply_scdd <- function(sce) {
    kids <- levels(sce$cluster_id)
    cells_by_k <- split(colnames(sce), sce$cluster_id)
    normcounts(sce) <- assay(sce)
    suppressMessages(
        res <- lapply(kids, function(k) {
            res <- results(scDD(sce[, cells_by_k[[k]]],
                min.nonzero = 20, condition = "group_id",
                categorize = FALSE, testZeroes = FALSE,
                param = BiocParallel::MulticoreParam(workers = 1)
            ))
            data.frame(
                gene = rownames(sce),
                cluster_id = k,
                p_val = res$nonzero.pvalue,
                p_adj.loc = res$nonzero.pvalue.adj,
                row.names = NULL,
                stringsAsFactors = FALSE
            )
        })
    )
    df <- bind_rows(res)
    df$p_adj.glb <- p.adjust(df$p_val)
    return(df)
}
