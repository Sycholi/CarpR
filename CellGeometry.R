library(cellGeometry)
mat = seu@assays$RNA@layers$counts  # JoinLayers()
rownames(mat) = rownames(seu)
meta = seu@meta.data
subcl = meta$subtype
mk = cellMarkers(mat,
                 subclass = subcl,
                 dual_mean = TRUE,
                 cores = 80)
signature_heatmap(mk)
spillover_heatmap(mk)
dput(unique(seu$subtype))

# TCGA-LIHC
load('TCGA-LIHC-DATA.RData')
count_matrix %<>% t()  # rownames: Genes/ colnames: Samples
count_matrix[1:5,1:5]
mk = updateMarkers(mk, bulkdata = count_matrix)
fit = deconvolute(mk, count_matrix, use_filter = FALSE)
fit

phenom %<>% dplyr::filter(submitter_id.samples %in% rownames(count_matrix))
df = data.frame(CAF = fit[["subclass"]][["percent"]][,2],
                Group = phenom$sample_type.samples)
df %<>% dplyr::filter(Group %in% c('Primary Tumor', 'Solid Tissue Normal'))

p1 = ggplot(data = df) +
  geom_boxplot(aes(y = CAF, colour = Group)) + 
  scale_color_manual(values = c("#F08A48","#376795")) 
ggsave(p1, filename = 'figures/pseudobulk_bar.pdf',
       height = 1500, width = 1500, units = 'px')
