suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
    library(scuttle)
})

aggregate_assay <- function(data, assay, method) {
    assays(data) <- list("values" = switch(assay,
        "counts" = counts(data),
        "cpm" = calculateCPM(counts(data)),
        "logcounts" = normalizeCounts(computeLibraryFactors(data)),
        "vstresiduals" = vst(counts(data), show_progress = FALSE)$y
    ))
    switch(method,
        "None" = data,
        "Mean" = aggregate_wrapper(x = data, fun = "mean"),
        "Sum" = aggregate_wrapper(x = data, fun = "sum")
    )
}



aggregate_wrapper <- function(x,
                              by = c("cluster_id", "sample_id"),
                              fun = c("sum", "mean")) {
    scale = FALSE
    fun <- match.arg(fun)
    assay <- assayNames(x)[1]
    for (i in by) if (!is.factor(x[[i]]))
        x[[i]] <- factor(x[[i]])
    suppressWarnings(y <- summarizeAssayByGroup(x, assay.type = assay,
                                                ids = (ids <- colData(x)[by]),
                                                statistics = fun,
                                                BPPARAM = BPPARAM = SerialParam()))
    colnames(y) <- y[[by[length(by)]]]
    if (length(by) == 1)
        return(assay(y))
    if (is.factor(ids <- y[[by[1]]]))
        ids <- droplevels(ids)
    is <- split(seq_len(ncol(y)), ids)
    pb <- map(is, ~assay(y)[, ., drop = FALSE])
    for (i in seq_along(pb)) {
        fill <- setdiff(unique(y[[by[2]]]), colnames(pb[[i]]))
        if (length(fill != 0)) {
            foo <- matrix(0, nrow(x), length(fill))
            colnames(foo) <- fill
            foo <- cbind(pb[[i]], foo)
            o <- paste(sort(unique(y[[by[2]]])))
            pb[[i]] <- foo[, o]
        }
    }
    md <- metadata(x)
    md$agg_pars <- list(assay = assay, by = by, fun = fun, scale = scale)
    pb <- SingleCellExperiment(pb, rowData = rowData(x), metadata = md)
    cd <- data.frame(colData(x)[, by])
    for (i in names(cd)) if (is.factor(cd[[i]]))
        cd[[i]] <- droplevels(cd[[i]])
    ns <- table(cd)
    if (length(by) == 2) {
        ns <- asplit(ns, 2)
        ns <- map(ns, ~c(unclass(.)))
    }
    else ns <- c(unclass(ns))
    int_colData(pb)$n_cells <- ns
    if (length(by) == 2) {
        cd <- colData(x)
        ids <- colnames(pb)
        counts <- vapply(ids, function(u) {
            m <- as.logical(match(cd[, by[2]], u, nomatch = 0))
            vapply(cd[m, ], function(u) length(unique(u)), numeric(1))
        }, numeric(ncol(colData(x))))
        cd_keep <- apply(counts, 1, function(u) all(u == 1))
        cd_keep <- setdiff(names(which(cd_keep)), by)
        if (length(cd_keep) != 0) {
            m <- match(ids, cd[, by[2]], nomatch = 0)
            cd <- cd[m, cd_keep, drop = FALSE]
            rownames(cd) <- ids
            colData(pb) <- cd
        }
    }
    return(pb)
}
