suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
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
