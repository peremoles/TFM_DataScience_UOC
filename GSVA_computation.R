library(GSVA)
library(BiocParallel)

spe <- readRDS("spe.rds")
genesets <- readRDS("genesets.rds")

gsvaPar <- gsvaParam(spe, genesets, assay = "logcounts", kcdf="none", minSize = 2)
system.time(res <- gsva(gsvaPar,  BPPARAM=MulticoreParam(workers=4)))

saveRDS(res, file="gsva_res.rds")
