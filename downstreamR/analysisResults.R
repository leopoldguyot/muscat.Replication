library(tidyverse)

metadata <- read.csv("dataOutput/modelResults/resultMetadata.csv")


suppressMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})
perf_metrics <- function(metadata, selectedData, selectedNC, selectedProps, padjLoc) {
    subset <- metadata %>%
        filter(data == selectedData, NC == selectedNC, props == selectedProps, model != "MMdream2")

    df_list <- lapply(seq_len(nrow(subset)), function(i) {
        file <- subset$path[i]
        print(file)
        df <- readRDS(file) %>%
            mutate(E = (sim_mean.A + sim_mean.B) / 2) %>%
            filter(E > 0.1) %>%
            mutate(sid = paste(gene, cluster_id, sep = "_"))

        method_name <- paste(subset$model[i], subset$aggregationType[i], subset$valueType[i], sep = "_")

        # Return a data frame with method column names
        data.frame(
            sid = df$sid,
            pval = df$p_val,
            padj.loc = df$p_adj.loc,
            padj.glb = df$p_adj.glb,
            is_de = df$is_de,
            method = method_name
        )
    })

    combined_df <- bind_rows(df_list)

    pval_df <- combined_df %>%
        select(sid, method, pval) %>%
        pivot_wider(names_from = method, values_from = pval) %>%
        column_to_rownames("sid")

    padj.loc_df <- combined_df %>%
        select(sid, method, padj.loc) %>%
        pivot_wider(names_from = method, values_from = padj.loc) %>%
        column_to_rownames("sid")
    padj.glb_df <- combined_df %>%
        select(sid, method, padj.glb) %>%
        pivot_wider(names_from = method, values_from = padj.glb) %>%
        column_to_rownames("sid")

    truth_df <- combined_df %>%
        select(sid, is_de) %>%
        distinct() %>%
        column_to_rownames("sid")

    cobdata <- COBRAData(
        padj = as.data.frame(if (padjLoc) padj.loc_df else padj.glb_df),
        truth = as.data.frame(truth_df)
    )
    perf <- calculate_performance(cobdata,
                                  binary_truth = "is_de",
                                  aspects = "fdrtpr",
                                  maxsplit = Inf,
                                  thrs = c(0.01, 0.05, 0.1, 0.2)
    ) %>%
        fdrtpr() %>%
        mutate(thr = as.numeric(sub("thr", "", thr)))
    return(perf)
}
make_prop_plot <- function(data) {
    ggplot(data, aes(x = FDR, y = TPR, color = method)) +
        geom_vline(
            xintercept = c(0.01, 0.05, 0.1),
            linetype = "dashed", color = "grey50", linewidth = 0.3
        ) +
        geom_point(size = 2.5, alpha = 0.8) +
        geom_line(size = 0.7) +
        scale_x_continuous(
            trans = "sqrt",
            limits = c(0, 1),
            breaks = c(0.01, 0.05, 0.2, 0.4, 0.8, 1),
            labels = scales::label_number()
        ) +
        scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, 0.2)
        ) +
        theme_minimal() +
        labs(
            x = "False Discovery Rate (sqrt scale)",
            y = "True Positive Rate",
            color = "Method"
        ) +
        facet_grid(rows = vars(padjType), cols = vars(prop)) +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            legend.position = "bottom"
        )

}
make_size_plot <- function(data) {
    ggplot(data, aes(x = FDR, y = TPR, color = NC)) +
        geom_vline(
            xintercept = c(0.01, 0.05, 0.1),
            linetype = "dashed", color = "grey50", linewidth = 0.3
        ) +
        geom_point(size = 2.5, alpha = 0.8) +
        geom_line(size = 0.7) +
        scale_x_continuous(
            trans = "sqrt",
            limits = c(0, 1),
            breaks = c(0.01, 0.05, 0.2, 0.4, 0.8, 1),
            labels = scales::label_number()
        ) +
        scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, 0.2)
        ) +
        theme_minimal() +
        labs(
            x = "False Discovery Rate (sqrt scale)",
            y = "True Positive Rate",
            color = "Method"
        ) + facet_wrap(~method) +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            legend.position = "bottom"
        )

}

## Process LPS results
perfsPropLPS <- lapply(unique(metadata$props), function(prop) {
    cobDataLoc <- perf_metrics(metadata, "simLPS", "NC400", prop, TRUE) %>%
        mutate(padjType = "local",
               prop = prop)
    cobDataGlb <- perf_metrics(metadata, "simLPS", "NC400", prop, FALSE) %>%
        mutate(padjType = "global",
               prop = prop)

    rbind(cobDataLoc, cobDataGlb)
})

perfsSizeLPS <- lapply(unique(metadata$NC), function(size) {
     perf_metrics(metadata, "simLPS", size, "de10", FALSE) %>%
        mutate(padjType = "global",
               NC = size)
})

## Process Kang results

perfsPropKang <- lapply(unique(metadata$props), function(prop) {
    cobDataLoc <- perf_metrics(metadata, "simKang", "NC400", prop, TRUE) %>%
        mutate(padjType = "local",
               prop = prop)
    cobDataGlb <- perf_metrics(metadata, "simKang", "NC400", prop, FALSE) %>%
        mutate(padjType = "global",
               prop = prop)

    rbind(cobDataLoc, cobDataGlb)
})

perfsSizeKang <- lapply(unique(metadata$NC), function(size) {
    perf_metrics(metadata, "simKang", size, "de10", FALSE) %>%
        mutate(padjType = "global",
               NC = size)
})

## Process testis results

perfsPropTestis <- lapply(unique(metadata$props), function(prop) {
    cobDataLoc <- perf_metrics(metadata, "simtestis", "NC400", prop, TRUE) %>%
        mutate(padjType = "local",
               prop = prop)
    cobDataGlb <- perf_metrics(metadata, "simtestis", "NC400", prop, FALSE) %>%
        mutate(padjType = "global",
               prop = prop)

    rbind(cobDataLoc, cobDataGlb)
})

perfsSizeTestis <- lapply(unique(metadata$NC), function(size) {
    perf_metrics(metadata, "simtestis", size, "de10", FALSE) %>%
        mutate(padjType = "global",
               NC = size)
})

## aggregate processed results

perfsPropLPSDF <- do.call(rbind, perfsPropLPS)
perfsSizeLPSDF <- do.call(rbind, perfsSizeLPS)
perfsPropKangDF <- do.call(rbind, perfsPropKang)
perfsSizeKangDF <- do.call(rbind, perfsSizeKang)
perfsPropTestisDF <- do.call(rbind, perfsPropTestis)
perfsSizeTestisDF <- do.call(rbind, perfsSizeTestis)


## plots LPS
plotPropLPS <- make_prop_plot(perfsPropLPSDF)

ggsave("figs/fdrtpr_prop_method_LPS.pdf", plotPropLPS, width = 9, height = 5)

plotSizeLPS <- make_size_plot(perfsSizeLPSDF)

ggsave("figs/fdrtpr_size_method_LPS.pdf", plotSizeLPS, width = 7, height = 5)

## plots Kang

plotPropKang <- make_prop_plot(perfsPropKangDF)

ggsave("figs/fdrtpr_prop_method_Kang.pdf", plotPropKang, width = 9, height = 5)

plotSizeKang <- make_size_plot(perfsSizeKangDF)


ggsave("figs/fdrtpr_size_method_Kang.pdf", plotSizeKang, width = 7, height = 5)

## plots testis

plotPropTestis <- make_prop_plot(perfsPropTestisDF)

ggsave("figs/fdrtpr_prop_method_Testis.pdf", plotPropTestis, width = 9, height = 5)

plotSizeTestis <- make_size_plot(perfsSizeTestisDF)

ggsave("figs/fdrtpr_size_method_Testis.pdf", plotSizeTestis, width = 7, height = 5)
