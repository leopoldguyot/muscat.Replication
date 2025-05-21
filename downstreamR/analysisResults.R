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
                                  thrs = seq(0.01, 0.5, by = 0.05)
    ) %>%
        fdrtpr() %>%
        mutate(thr = as.numeric(sub("thr", "", thr)))
    return(perf)
}

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

perfsPropLPSDF <- do.call(rbind, perfsPropLPS)
perfsSizeLPSDF <- do.call(rbind, perfsSizeLPS)
perfsPropKangDF <- do.call(rbind, perfsPropKang)
perfsSizeKangDF <- do.call(rbind, perfsSizeKang)

plot <- ggplot(perfsPropLPSDF, aes(x = FDR, y = TPR, color = method)) +
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
        title = "FDR vs TPR (square root x-axis)",
        x = "False Discovery Rate (sqrt scale)",
        y = "True Positive Rate",
        color = "Method"
    ) + facet_grid(rows = vars(padjType), cols = vars(prop)) + theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )

ggsave("figs/fdrtpr_prop_method_LPS.pdf", plot)


plot <- ggplot(perfsSizeLPSDF, aes(x = FDR, y = TPR, color = NC)) +
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
        title = "FDR vs TPR (square root x-axis)",
        x = "False Discovery Rate (sqrt scale)",
        y = "True Positive Rate",
        color = "Method"
    ) + facet_wrap(~method) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )

ggsave("figs/fdrtpr_size_method_LPS.pdf", plot)

plot <- ggplot(perfsPropKangDF, aes(x = FDR, y = TPR, color = method)) +
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
        title = "FDR vs TPR (square root x-axis)",
        x = "False Discovery Rate (sqrt scale)",
        y = "True Positive Rate",
        color = "Method"
    ) + facet_grid(rows = vars(padjType), cols = vars(prop)) + theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )

ggsave("figs/fdrtpr_prop_method_Kang.pdf", plot)


plot <- ggplot(perfsSizeKangDF, aes(x = FDR, y = TPR, color = NC)) +
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
        title = "FDR vs TPR (square root x-axis)",
        x = "False Discovery Rate (sqrt scale)",
        y = "True Positive Rate",
        color = "Method"
    ) + facet_wrap(~method) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )

ggsave("figs/fdrtpr_size_method_Kang.pdf", plot)
