# Alex Abbas
# 170103
# load project, save as eset data

library(ExpressionPlot) # internal Genentech package for RNAseq project management

ngs1030 = ep.find.project("NGS1030")
eset = ep.ExpressionSet(proj = ngs1030, feature.type = "gene", stat = "rpkm", attach.annot = T)
save(eset,file="~/ipf/ipfcore.git/data/eset.RData")
