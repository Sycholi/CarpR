## TCR readout ----
suppressPackageStartupMessages(library(scRepertoire))

tcr_sample = list.files("original_data/TCRs/")

contig.list = lapply(tcr_sample, function(sid){
  loadContigs(read.csv(file.path("original_data/TCRs", sid, "filtered_contig_annotations.csv"),
             stringsAsFactors = FALSE), format = "10X")
})
names(contig.list) = tcr_sample
combined.tcr = combineTCR(input.data = contig.list, samples = tcr_sample)
seu = combineExpression(input.data = combined.tcr, sc.data = seu, cloneCall = "strict", chain = "both")
table(is.na(seu$CTstrict))
DimPlot(seu, group.by = 'cloneSize')
