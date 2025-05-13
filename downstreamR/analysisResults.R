library(tidyverse)

metadata <- read.csv("dataOutput/modelResults/resultMetadata.csv")


suppressMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})
curve_NC400_DP <- function(metadata) {
    subset <- metadata %>%
        filter(data == "simKang", NC == "NC400", props == "de10", model != "scDD")

    df_list <- lapply(seq_len(nrow(subset)), function(i) {
        file <- subset$path[i]
        df <- readRDS(file) %>%
            mutate(E = (sim_mean.A + sim_mean.B) / 2) %>%
            filter(E > 0.1) %>%
            mutate(sid = paste(gene, cluster_id, sep = "_"))

        method_name <- paste(subset$model[i], subset$aggregationType[i], subset$valueType[i], sep = "_")

        # Return a data frame with method column names
        data.frame(
            sid = df$sid,
            pval = df$p_val,
            padj = df$p_adj.loc,
            is_de = df$is_de,
            method = method_name
        )
    })

    combined_df <- bind_rows(df_list)

    pval_df <- combined_df %>%
        select(sid, method, pval) %>%
        pivot_wider(names_from = method, values_from = pval) %>%
        column_to_rownames("sid")

    padj_df <- combined_df %>%
        select(sid, method, padj) %>%
        pivot_wider(names_from = method, values_from = padj) %>%
        column_to_rownames("sid")

    truth_df <- combined_df %>%
        select(sid, is_de) %>%
        distinct() %>%
        column_to_rownames("sid")

    cobdata <- COBRAData(
        pval = as.data.frame(pval_df),
        padj = as.data.frame(padj_df),
        truth = as.data.frame(truth_df)
    )

    return(cobdata)
}

t <- curve_NC400_DP(metadata)
perf <- calculate_performance(t,
                              binary_truth = "is_de",
                              aspects = "fdrtpr",
                              maxsplit = Inf
                              ) %>%
    fdrtpr() %>%
    mutate(thr = as.numeric(sub("thr", "", thr)))

plot <- ggplot(perf, aes(x = FDR, y = TPR, color = method)) +
    geom_vline(xintercept = c(0.01, 0.05, 0.1),
               linetype = "dashed", color = "grey50", size = 0.3) +
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
    )

ggsave("figs/NC400_DP.pdf", plot)
