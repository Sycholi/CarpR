library(clusterProfiler);library(enrichplot)
library(org.Hs.eg.db);library(org.Mm.eg.db)
markers = markers_from_FindMarkers
markers %<>% rownames_to_column('Gene')
markers %<>% arrange(-avg_log2FC)

gene_symbol = bitr(markers$Gene,
                   fromType = 'SYMBOL',
                   toType = 'ENTREZID',
                   OrgDb = 'org.Hs.eg.db')
markers = merge(x = gene_symbol, by.x = 'SYMBOL',
                y = markers, by.y = 'Gene', all.x = T, all.y = F, sort = F)
KEGG_result = enrichKEGG(gene = markers$ENTREZID,
                         organism = 'hsa')
KEGG_result = setReadable(KEGG_result, keyType = 'ENTREZID',OrgDb = org.Hs.eg.db)
KEGG_result@result %<>% filter(category != 'Human Diseases')
KEGG_result@result %<>% subset(subcategory %in% c('Immune system','Signal transduction','Energy metabolism'))
# KEGG_result@result %<>% subset(subcategory %in% c('Immune system','Signal transduction','Energy metabolism'))
p1 = ggplot(KEGG_result, showCategory = 10, 
            aes(zScore, fct_reorder(Description, zScore))) + 
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_gradient2(low = '#F08A48', mid = '#F08A48', high = '#376795') +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  ylab(NULL) ;p1
ggsave(p1, filename = paste0('figures/KEGG_enrich.pdf'),
       height = 1200, width = 2000, units = 'px', dpi = 300)
