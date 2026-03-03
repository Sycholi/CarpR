
deg = markers %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5)
geneList = split(deg$gene, deg$cluster)
all.genes = unique(deg$gene)
id.map = bitr(all.genes,
              fromType = "SYMBOL",
              toType   = "ENTREZID",
              OrgDb    = org.Hs.eg.db)
geneList.entrez = lapply(geneList, function(gs){id.map$ENTREZID[match(gs, id.map$SYMBOL)] %>% na.omit()})

# Enrichment Analysis: compareCluster + enrichGO
ego.compare = compareCluster(
  geneCluster   = geneList.entrez,
  fun           = "enrichKEGG",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2)
ego.compare@compareClusterResult %<>% filter(category != 'Human Diseases')

# Visualization
dp = dotplot(ego.compare,
             showCategory = 10, 
             x = "Cluster",         
             color = "p.adjust") +
  scale_color_gradient(low = "#376795", high = "#F08A48") +
  # coord_flip() + 
  theme_bw(base_size = 12) +
  theme(
    axis.text.y  = element_text(size = 10),    
    axis.text.x  = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
dp
ggsave(dp, filename = 'figures/Clusters_KEGG.pdf',
       height = 3000, width = 2000, units = 'px')
