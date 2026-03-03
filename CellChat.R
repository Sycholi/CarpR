seu %<>% JoinLayers()
library(CellChat)
seu.NC = seu %>% subset(Group == 'SOFT')
seu.TR = seu %>% subset(Group == 'STIFF')

cellchat.NC = createCellChat(object = seu.NC, group.by = "subtype", assay = 'RNA')
cellchat.TR = createCellChat(object = seu.TR, group.by = "subtype", assay = 'RNA')
# saveRDS(cellchat.NC, 'cellchat.NC.rds')
# saveRDS(cellchat.TR, 'cellchat.TR.rds')

CellChatDB = CellChatDB.human
CellChatDB.use = subsetDB(CellChatDB, search = c("Secreted Signaling","Cell-Cell Contact"), key = "annotation")
cellchat.NC@DB = CellChatDB.use
cellchat.TR@DB = CellChatDB.use

source('~/r_repo/repository/cellchat_process.r')

cellchat.NC %<>% cellchat_process()
cellchat.TR %<>% cellchat_process()

object.list = list(TR = cellchat.TR, NC = cellchat.NC)
cellchat = mergeCellChat(object.list, add.names = names(object.list))
cellchat
saveRDS(cellchat, file = 'cellchat.rds')

gg1 = compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 = compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, sources.use = 2)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", sources.use = 2)

gg1 = netVisual_heatmap(cellchat)
gg2 = netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

num.link = sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax = c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg = list()
for (i in 1:length(object.list)) {
  gg[[i]] = netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg) + xlim(c(0,15)) + ylim(c(0,20))

netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD4_C01_CTLA4_Treg", signaling = c('CXCL','CCL'))
gg2 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD4_C01_CTLA4_Treg", signaling.exclude = c("MIF"))
patchwork::wrap_plots(plots = list(gg1,gg2))

p1 = gg1+gg2
ggsave(p1, filename = 'tmp.pdf', height = 1500, width = 3000, units = 'px', dpi = 300)

cellchat = computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat = netEmbedding(cellchat, type = "functional")
cellchat = netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

rankSimilarity(cellchat, type = "functional")

gg1 = rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = 2, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 = rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = 2, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union = union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 10, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 10, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

netVisual_bubble(cellchat, 
                 sources.use = 2, 
                 targets.use = 12:13,
                 comparison = c(2, 1), 
                 angle.x = 45)

gg1 = netVisual_bubble(cellchat, sources.use = 2, targets.use = 8:13,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in TR", angle.x = 45, remove.isolate = T)
gg2 = netVisual_bubble(cellchat, sources.use = 2, targets.use = 8:13,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in NC", angle.x = 45, remove.isolate = T)
gg1 + gg2

gg1 = netVisual_aggregate(cellchat.NC, signaling = 'CXCL', sources.use = 2, targets.use = 1:20)
gg2 = netVisual_aggregate(cellchat.TR, signaling = 'CXCL', sources.use = 2, targets.use = 1:20)

pairLR.CXCL = extractEnrichedLR(cellchat, signaling = 'CXCL', geneLR.return = FALSE)
LR.show <- pairLR.CXCL[5,] # show one ligand-receptor pair
# Hierarchy plot
netVisual_individual(cellchat.TR, signaling = 'CXCL',  pairLR.use = LR.show, vertex.receiver = c(2,3,12,13,15,16))
#> [[1]]
# Circle plot
netVisual_individual(cellchat.TR, signaling = 'CXCL', pairLR.use = LR.show, layout = "chord")
netVisual_individual(cellchat.NC, signaling = 'CXCL', pairLR.use = LR.show, layout = "chord")

pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 2, targets.use = 8:13, lab.cex = 0.5, title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
}

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NC", "TR")) # set factor level
pp = plotGeneExpression(cellchat, 
                        signaling = c("CXCL"), 
                        split.by = "datasets", 
                        colors.ggplot = T, 
                        type = "violin",
                        color.use = c("#4575B4","#A50026"));pp

VlnPlot(seu, 
        c('Cxcl16','Cxcr6'),
        stack = T,
        split.by = 'Group',
        cols = c("#4575B4","#A50026"))
