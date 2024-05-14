library(nnSVG)
library(BiocParallel)

gsva_res <- readRDS("gsva_res.rds")
nnsvg_res <- nnSVG(gsva_res, assay_name = "es", BPPARAM=MulticoreParam(workers=4))
saveRDS(nnsvg_res, file="nnsvg_res.rds")
