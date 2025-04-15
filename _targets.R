# _targets.R
library(targets)
library(tarchetypes)

tar_option_set(
    packages = c("muscat", "ExperimentHub", "muscData", "SingleCellExperiment")
)

tar_source()

target_analysis_pipeline <- function(data) {
    sim_sym <- as.symbol(paste0("sim", data))

    list(
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
    list(
        tar_target(name = prepKang, command = prep_data("Kang")),
        tar_target(name = prepLPS, command = prep_data("LPS")),
        tar_target(name = simKang, command = simulate_data(prepKang)),
        tar_target(name = simLPS, command = simulate_data(prepLPS))
    ),
    target_analysis_pipeline("Kang"),
    target_analysis_pipeline("LPS")
)
