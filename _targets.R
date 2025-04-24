# _targets.R
library(targets)
library(tarchetypes)

tar_option_set(
    packages = c(
        "muscat",
        "ExperimentHub",
        "muscData",
        "SingleCellExperiment",
        "scater",
        "sctransform",
        "tidyr"
    )
)

tar_source()

target_analysis_base <- function(data) {
    prep_sym <- as.symbol(paste0("prep", data))
    sim_sym <- as.symbol(paste0("sim", data))
    list(
        tar_target_raw(
            name = paste0("prep", data),
            command = substitute(prep_data(data))
        ),
        target_analysis_node1(data, list(name = "NC20_nill",
                                         ng = 4e3,
                                         nc = 16*20,
                                         ns = NULL,
                                         nk = NULL,
                                         p_dd = diag(6)[1, ],
                                         probs = NULL)),
        target_analysis_node1(data, list(name = "NC100_nill",
                                         ng = 4e3,
                                         nc = 16*100,
                                         ns = NULL,
                                         nk = NULL,
                                         p_dd = diag(6)[1, ],
                                         probs = NULL)),
        target_analysis_node1(data, list(name = "NC400_nill",
                                         ng = 4e3,
                                         nc = 16*400,
                                         ns = NULL,
                                         nk = NULL,
                                         p_dd = diag(6)[1, ],
                                         probs = NULL)),
        target_analysis_node1(data, list(name = "NC400_de10",
                                         ng = 4e3,
                                         nc = 16*400,
                                         ns = NULL,
                                         nk = NULL,
                                         p_dd = c(0.9, 0, 0.1, 0, 0, 0),
                                         probs = NULL)),
        target_analysis_node1(data, list(name = "NC400_dp10",
                                         ng = 4e3,
                                         nc = 16*400,
                                         ns = NULL,
                                         nk = NULL,
                                         p_dd = c(0.9, 0, 0, 0.1, 0, 0),
                                         probs = NULL)),
        target_analysis_node1(data, list(name = "NC400_dm10",
                                         ng = 4e3,
                                         nc = 16*400,
                                         ns = NULL,
                                         nk = NULL,
                                         p_dd = c(0.9, 0, 0, 0, 0.1, 0),
                                         probs = NULL)),
        target_analysis_node1(data, list(name = "NC400_db10",
                                         ng = 4e3,
                                         nc = 16*400,
                                         ns = NULL,
                                         nk = NULL,
                                         p_dd = c(0.9, 0, 0, 0, 0, 0.1),
                                         probs = NULL))
    )
}

target_analysis_node1 <- function(data, params) {
    prep_sym <- as.symbol(paste0("prep", data))
    sim_sym <- as.symbol(paste0("sim", data, "_", params$name))
    list(
        tar_target_raw(
            name = paste0("sim", data, "_", params$name),
            command = substitute(simulate_data(PREP,
                                               NG,
                                               NC,
                                               NS,
                                               NK,
                                               PDD,
                                               PROBS),
                                 list(PREP = prep_sym,
                                      NG = params$ng,
                                      NC = params$nc,
                                      NS = params$ns,
                                      NK = params$nk,
                                      PDD = params$p_dd,
                                      PROBS = params$probs))
        ),
        target_analysis_node2(sim_sym, list(assay = "counts",
                                            method = "Sum")),
        target_analysis_node2(sim_sym, list(assay = "cpm",
                                            method = "Sum")),
        target_analysis_node2(sim_sym, list(assay = "logcounts",
                                            method = "Mean")),
        target_analysis_node2(sim_sym, list(assay = "vstresiduals",
                                            method = "Mean")),
        target_analysis_node2(sim_sym, list(assay = "vstresiduals",
                                            method = "None")),
        target_analysis_node2(sim_sym, list(assay = "logcounts",
                                            method = "None"))
    )
}

target_analysis_node2 <- function(sim, params) {
    list(tar_target_raw(
        name = paste0(as.character(sim), "_", params$assay, "_", params$method),
        command = substitute(aggregate_assay(data = SIM,
                                             assay = ASSAY,
                                             method = METHOD),
                             list(SIM = sim,
                                  ASSAY = params$assay,
                                  METHOD = params$method))
    ))
}

# Return the full target list:
c(
    target_analysis_base("Kang"),
    target_analysis_base("LPS")
)
