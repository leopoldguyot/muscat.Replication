suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
})
prep_data <- function(data) {
    stopifnot(data %in% c("Kang", "LPS"))
    switch(data,
        "Kang" = prep_Kang_data(),
        "LPS" = prep_LPS_data()
    )
}

prep_Kang_data <- function() {
    eh <- ExperimentHub()
    sce <- eh[["EH2259"]]
    reducedDims(sce) <- NULL # remove dimensionality reductions
    sce <- sce[, sce$multiplets == "singlet"] # remove multiplets
    sce <- sce[, !is.na(sce$cell)] # remove unassigned cells
    sce <- sce[, sce$stim == "ctrl"] # keep control samples only
    sce$sample_id <- factor(paste0(sce$stim, sce$ind))

    sce <- prepSCE(sce, "cell", "sample_id", "stim", TRUE)

    sce <- prepSim(sce,
        verbose = TRUE,
        # keep genes w/ count > 1 in >= 10 cells
        min_count = 1, min_cells = 10,
        # keep cells w/ >= 100 detected genes & cluster w/ > 100 cells
        min_genes = 100, min_size = 100
    )

    sce
}

prep_LPS_data <- function() {
    eh <- ExperimentHub()
    sce <- eh[["EH3297"]]
    reducedDims(sce) <- NULL # remove dimension reductions
    assays(sce) <- SimpleList(counts = counts(sce)) # remove slots other than counts
    sce <- sce[, sce$group_id == "Vehicle"] # keep reference samples only
    sce <- prepSCE(sce, "cluster_id", "sample_id", "group_id", TRUE) # prep. SCE for 'muscat'
    # prep. SCE for simulation w/ 'muscat::simData'
    sce <- prepSim(sce,
        verbose = FALSE,
        # keep genes w/ count > 1 in >= 10 cells
        min_count = 1, min_cells = 10,
        # keep cells w/ >= 100 detected genes & cluster w/ > 100 cells
        min_genes = 100, min_size = 100
    )

    sce
}


simulate_data <- function(prepData,
                          ng,
                          nc,
                          ns,
                          nk,
                          p_dd,
                          probs) {
    set.seed(124)

    sim <- simData(prepData,
        paired = FALSE, lfc = 2,
        ng = nrow(prepData), nc = nc,
        ns = ns, nk = nk,
        p_dd = p_dd, probs = probs
    )

    sim <- sim[rowSums(counts(sim) > 0) >= 10, ]
    sim <- sim[sample(nrow(sim), min(nrow(sim), ng)), ]

    gi <- metadata(sim)$gene_info
    gi <- dplyr::filter(gi, gene %in% rownames(sim))
    metadata(sim)$gene_info <- gi

    sim <- computeLibraryFactors(sim)
    sim <- logNormCounts(sim)
    assays(sim)$cpm <- calculateCPM(sim)
    assays(sim)$vstresiduals <- suppressWarnings(
        vst(counts(sim), show_progress = FALSE)$y
    )
    sim
}

simulate_testis_data <- function(prepData) {
    library(SingleCellExperiment)
    library(CTdata)


    sce <- CTdata::testis_sce()
}
