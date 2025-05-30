library(dplyr)
library(tidyr)
library(ggplot2)
library(targets)

tar_runtime <- tar_meta(fields = c("name", "seconds")) %>%
  filter(!is.na(seconds)) %>%
  separate(
    name,
    into = c("dataset", "NC", "prop", "counts", "aggregation", "model"),
    sep = "_",
    fill = "right",
    remove = FALSE  # Keep original name column if you want to inspect
  )

tar_runtime_filtered <- tar_runtime %>%
  filter(NC == "NC400",
         prop == "de10",
         dataset == "simKang")

tar_runtime_model <- tar_runtime_filtered[!is.na(tar_runtime_filtered$model) & tar_runtime_filtered$model != "MMdream2",]
tar_runtime_counts <- tar_runtime_filtered[is.na(tar_runtime_filtered$model) & !is.na(tar_runtime_filtered$aggregation),]
tar_runtime_total <- tar_runtime_model %>%
  left_join(
    tar_runtime_counts %>%
      select(dataset, NC, prop, counts, aggregation, seconds) %>%
      rename(counts_seconds = seconds),
    by = c("dataset", "NC", "prop", "counts", "aggregation")
  ) %>%
  mutate(total_seconds = counts_seconds + seconds) %>%
  unite(method, counts, aggregation, model, sep = "_", na.rm = TRUE)

# View final result
plotKang <- tar_runtime_total %>%
  select(dataset, NC, prop, method, counts_seconds, model_seconds = seconds, total_seconds) %>%
  ggplot(aes(x = reorder(method, total_seconds), y = total_seconds, fill = method)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +  # Horizontal bars for better readability
  labs(
    title = "Total Runtime by Method",
    x = "Method (Counts_Aggregation_Model)",
    y = "Total Time (seconds)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
ggsave("figs/plot_runtime_Kang.pdf", plotKang, width = 6, height = 4)
