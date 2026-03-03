library(monocle3)
cds = as.cell_data_set(seu)
cds = cluster_cells(cds, reduction_method = 'UMAP')
cds = learn_graph(cds, use_partition = T)
expr = logcounts(cds)
naive_genes = intersect(c("SELL","CCR7","TCF7","LEF1","IL7R"), rownames(expr))  # indicate starter
naive_score = Matrix::colMeans(expr[naive_genes, , drop = FALSE])
thr = quantile(naive_score, 0.99, na.rm = TRUE)
root_cells = names(naive_score)[naive_score >= thr]
cds = order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)
plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)

