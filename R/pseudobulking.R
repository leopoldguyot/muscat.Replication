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
        "Mean" = aggregateData(x = data, fun = "mean"),
        "Sum" = aggregateData(x = data, fun = "sum")
    )
}
