# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(tidyr)

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
valuesData <- list(
  dataset = c("Kang", "LPS")
)
valuesAgg <- list(
  aggregation = c("None","Mean", "Sum")
)

targets <- tar_map(
  values = valuesData,
  tar_map(
    values = valuesAgg,
    tar_target(
      name = prepData,
      command = prep_data(dataset)
    ),
    # Dynamically simulate on each dataset
    tar_target(
      simData,
      simulate_data(prepData),
    ),
    tar_target(
      name = aggregateData,
      command = aggregate_assay(simData, aggregation))
  )
)
  
list(targets)
