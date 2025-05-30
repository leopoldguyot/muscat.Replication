library(stringr)
library(dplyr)
library(tidyr)
library(targets)

tarDF <- tar_manifest()
tarDF_filtered <- tarDF %>%
    filter(str_count(name, "_") >= 5) %>%
    separate(
        name,
        into = c("data", "NC", "props", "valueType", "aggregationType", "model"),
        sep = "_",
        extra = "merge",
        fill = "right",
        remove = FALSE
    )

paths <- sapply(tarDF_filtered$name, function(name) {
    object <- tar_read_raw(name)
    cat("\nsaving object: ", name)
    path <- file.path("dataOutput", "modelResults", paste0(name, ".rds"))
    saveRDS(object = object, path)
    path
})

tarDF_filtered[["path"]] <- paths

write.csv(tarDF_filtered, "dataOutput/modelResults/resultMetadata.csv")
