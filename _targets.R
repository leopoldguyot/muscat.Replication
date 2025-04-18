# _targets.R
library(targets)
library(tarchetypes)

tar_option_set(
    packages = c("muscat", "ExperimentHub", "muscData", "SingleCellExperiment")
)

tar_source()

target_analysis_pipeline <- function(data) {
    prep_sym <- as.symbol(paste0("prep", data))
    sim_sym <- as.symbol(paste0("sim", data))
    list(
      tar_target_raw(name = paste0("prep", data),
                     command = substitute(prep_data(data))),
      tar_target_raw(name = paste0("sim", data),
                     command = substitute(simulate_data(PREP), list(PREP = prep_sym))),
        tar_target_raw(
            name = paste0(data, "_aggregateSum"),
            command = substitute(aggregate_assay(data = SIM, method = "Sum"), list(SIM = sim_sym))
        ),
        tar_target_raw(
            name = paste0(data, "_aggregateMean"),
            command = substitute(aggregate_assay(data = SIM, method = "Mean"), list(SIM = sim_sym))
        )
    )
}

# Return the full target list:
c(
    target_analysis_pipeline("Kang"),
    target_analysis_pipeline("LPS")
)
