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

# Replace the target list below with your own:
list(
  tar_target(
    name = prepData,
    command = list(prep_Kang_data(),  prep_LPS_data())
  ),
  # Dynamically simulate on each dataset
  tar_target(
    simData,
    simulate_data(prepData),
    pattern = map(prepData)
  ),
  tar_target(
    aggMeanData,
    aggregate_assay(simData, "mean"),
    pattern = map(simData)
  ),
  tar_target(
    aggSumData,
    aggregate_assay(simData, "sum"),
    pattern = map(simData)
  )
)
