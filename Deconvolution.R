# CellGeometry ----
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




# MCP counter ----
library(MCPcounter)
exprs = seu@assays$RNA@layers$counts
pd = seu@meta.data
Estimates = MCPcounter.estimate(exprs, featuresType = "HUGO_symbols")

df = data.frame(t(Estimates)) %>%
  mutate(group = pd$Group) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(cols = colnames(.)[2:11],
               names_to = "cell.type",
               values_to = 'value')
df %<>% filter(cell.type == 'NK.cells')

p1 = ggplot(df,aes(cell.type,value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 80,vjust = 0.5,size = 14,face = "italic", colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 14,colour = 'black'))+
  scale_fill_manual(values = c('#1B3955','#7DB391')) +
  stat_compare_means(aes(group = group),label = "p.format",size = 3,method = "kruskal.test")
ggsave(p1,
       filename = 'MCP.pdf',
       height = 1500,
       width = 1000,
       units = 'px')




# CIBERSORT ----
library(CIBERSORT)
sig_matrix = AggregateExpression(seu, group.by = 'Group')
# sig_matrix = system.file("extdata", "LM22.txt", package = "CIBERSORT")
mixture = as.matrix(exprs)
res = cibersort(sig_matrix = sig_matrix, mixture_file = mixture, perm = 1000, QN = T)
saveRDS(res, 'res.rds')

res %<>% as.data.frame()

res$NK = res$`NK cells activated` + res$`NK cells resting`
res = res[,c(26,1:25)]
res1 = data.frame(res[,1:22]) %>%
  mutate(group = pd$Group) %>%
  rownames_to_column("Sample")%>%
  pivot_longer(cols = colnames(.)[2:23],
               names_to = "cell.type",
               values_to = 'value')

ggplot(res1,aes(cell.type,value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 14,face = "italic",colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 14,colour = 'black'))+
  scale_fill_nejm()+
  stat_compare_means(aes(group = group),label = "p.format",size=3,method = "kruskal.test")
