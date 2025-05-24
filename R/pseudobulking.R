suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
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
                                assay = NULL,
                                by = c("cluster_id", "sample_id"),
                                fun = c("sum", "mean")) {
    fun <- match.arg(fun)
        if (is.null(assay)) {
        assay <- assayNames(x)[1]
    }
    for (i in by) {
        if (!is.factor(x[[i]])) {
            x[[i]] <- factor(x[[i]])
        }
    }
    mat <- assay(x, assay)
    groups <- interaction(lapply(by, function(f) x[[f]]), drop = TRUE)
    agg_fun <- switch(fun,
                      "sum" = function(m) rowsum(m, group = groups),
                      "mean" = function(m) rowsum(m, group = groups) / as.vector(table(groups)))

    agg_mat <- agg_fun(t(mat))  # transpose because rowsum aggregates by rows
    agg_mat <- t(agg_mat)       # transpose back

    new_sce <- SingleCellExperiment(list(agg = agg_mat),
                                    rowData = rowData(x))
    colnames(new_sce) <- levels(groups)
    return(new_sce)
}
