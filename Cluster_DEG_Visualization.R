library(ggrepel)

# Markers from FindAllMarkers
markers = markers %>%
  filter(p_val_adj < 0.05) %>%
  mutate(label = ifelse(avg_log2FC < 0, "sigDown", "sigUp"))

# Labelled genes
topgene = markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>% 
  bind_rows(markers_fib %>%
      group_by(cluster) %>%
      slice_min(avg_log2FC, n = 5)) %>%
  ungroup()

# Colors
col_par = markers %>%
  group_by(cluster) %>%
  summarise(
    low = round(min(avg_log2FC) - 0.5),
    up  = round(max(avg_log2FC) + 0.5)
  )
scatter_col = c("#376795", "#F08A48")  # Color of Dots
ctys = unique(markers$cluster)
barcol = scales::hue_pal()(length(ctys))  # Color of Clusters


# Plotting
p1 = ggplot() +
  # Background
  geom_col(aes(x = cluster, y = low), col_par, fill = "#dcdcdc", alpha = 0.6) +
  geom_col(aes(x = cluster, y = up), col_par, fill = "#dcdcdc", alpha = 0.6) +
  # Dots
  geom_jitter(aes(x = cluster, y = avg_log2FC, color = label), markers_fib,
              width = 0.4, size = 1) +
  scale_color_manual(values = scatter_col) +
  # Cluster Bars
  geom_tile(aes(x = ctys, y = 0), height = 0.5, fill = barcol, show.legend = FALSE) +
  geom_text(aes(x = ctys, y = 0, label = ctys), size = 3, fontface = "bold") +
  # Labelled Genes
  geom_text_repel(aes(x = cluster, y = avg_log2FC, label = gene), 
                  data = topgene, size = 3, max.overlaps = Inf) +
  # Axis & Themes
  labs(x = "Cluster", y = "Average log2FoldChange", title = "Differential expression genes") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.98, 0.96),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 12)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

p1

ggsave(p1, filename = 'figures/DEG_Clusters.pdf', height = 2000, width = 3000, units = 'px')
