# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)

# Set target options:
tar_option_set(
  packages = c("muscat",
               "ExperimentHub",
               "muscData",
               "SingleCellExperiment"
               ) # Packages that your targets need for their tasks.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

target_analysis_pipeline <- function(data) {
  c(
    tar_target_raw(
      paste0(data, "_aggregateSum"),
      quote(aggregate_assay(paste0("sim",data), sum))
    ),
    tar_target_raw(
      paste0(data, "_aggregateMean"),
      quote(aggregate_assay(paste0("sim",data), mean))
    )
  )
}

list(
  tar_target(
    name = prepKang,
    command = prep_data("Kang")
  ),
  tar_target(
    name = prepLPS,
    command = prep_data("LPS")
  ),
  tar_target(
    name = simKang,
    command = simulate_data(prepKang)
  ),
  tar_target(
    name = simLPS,
    command = simulate_data(prepLPS)
  ),
  target_analysis_pipeline("Kang"),
  target_analysis_pipeline("LPS")
)
