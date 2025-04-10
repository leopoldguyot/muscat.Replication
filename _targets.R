# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("muscat", "ExperimentHub") # Packages that your targets need for their tasks.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_target(
    name = prepKang,
    command = prep_Kang_data()
  ),
  tar_target(
    name = prepLPS,
    command = prep_LPS_data()
  ),
  # Combine the prepped datasets into a list
  tar_target(prepList, list(prepKang, prepLPS)),
  # Dynamically simulate on each dataset
  tar_target(
    simData,
    simulate_data(prepData),
    pattern = map(prepData = prepList)
  )
)
