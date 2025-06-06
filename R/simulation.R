suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
    library(ExperimentHub)
    library(CTdata)
})
prep_data <- function(data) {
    stopifnot(data %in% c("Kang", "LPS", "testis"))
    switch(data,
        "Kang" = prep_Kang_data(),
        "LPS" = prep_LPS_data(),
        "testis" = prep_testis_data()
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


simulate_data <- function(
        prepData,
        ng,
        nc,
        ns,
        nk,
        p_dd,
        probs) {
    set.seed(124)
    paired <- FALSE
    if(length(unique(colData(prepData)$sample_id)) <= 3)
        paired <- TRUE
    sim <- simData(prepData,
        paired = paired, lfc = 2,
        ng = nrow(prepData), nc = nc,
        ns = ns, nk = nk,
        p_dd = p_dd, probs = probs
    )

    sim <- sim[rowSums(counts(sim) > 0) >= 10, ]
    sim <- sim[sample(nrow(sim), min(nrow(sim), ng)), ]

    gi <- metadata(sim)$gene_info
    gi <- dplyr::filter(gi, gene %in% rownames(sim))
    metadata(sim)$gene_info <- gi
    sim
}

prep_testis_data <- function() {
    sce <- CTdata::testis_sce()
    reducedDims(sce) <- NULL # remove dimensionality reductions
    assays(sce) <- SimpleList(counts = counts(sce)) # remove slots other than counts
    colData(sce)[["group_id"]] <- "CTRL"
    sce <- sce[, !is.na(sce$clusters)]


    sce <- prepSCE(sce, kid = "clusters",
            sid = "Donor",
            gid = "group_id")
    sce <- prepSim(sce,
                   verbose = FALSE,
                   # keep genes w/ count > 1 in >= 10 cells
                   min_count = 1, min_cells = 10,
                   # keep cells w/ >= 100 detected genes & cluster w/ > 100 cells
                   min_genes = 100, min_size = 100
    )
    sce
}
