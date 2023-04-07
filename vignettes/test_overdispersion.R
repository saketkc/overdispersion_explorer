suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(glmGamPoi)
  library(tidyverse)
  library(ggridges)
})

theme_set(theme_classic(base_size = 14))

# load datasets (both same technology)

# I expect this to not be as heterogeneous (theta should be high)
chrom.v3.pbmc <- readRDS("~/github/scRNA_NB_comparison/data/rds_filtered/PBMC__ChromiumV3.rds")

# I expect this to be more heterogenous (theta should be low)
chrom.v3.hek <- readRDS("~/github/scRNA_NB_comparison/data/rds_filtered/HEK__ChromiumV3.rds")

umi.pbmc <- GetAssayData(object = chrom.v3.pbmc, assay = "RNA", slot = "counts")
umi.hek <- GetAssayData(object = chrom.v3.hek, assay = "RNA", slot = "counts")


vst.out.pbmc <- sctransform::vst(umi = umi.pbmc, n_genes = nrow(umi.pbmc), n_cells = ncol(umi.pbmc), vst.flavor = "v2", return_gene_attr = T, return_cell_attr = T, return_corrected_umi = F)
vst.out.hek <- sctransform::vst(umi = umi.hek, n_genes = nrow(umi.hek), n_cells = ncol(umi.hek), vst.flavor = "v2", return_gene_attr = T, return_cell_attr = T, return_corrected_umi = F)


model_parameters.pbmc <- vst.out.pbmc$model_pars
model_parameters.hek <- vst.out.hek$model_pars



common.genes <- intersect(rownames(model_parameters.hek), rownames(model_parameters.pbmc))

df.default <- data.frame(hek=model_parameters.hek[common.genes, "theta"],
                         pbmc=model_parameters.pbmc[common.genes, "theta"],
                         gene = common.genes)

# replace large theta values by 10000
df.default[df.default$hek>1000, "hek"] <- 1000
df.default[df.default$pbmc>1000, "pbmc"] <- 1000

df.long <- df.default %>% pivot_longer(!gene, names_to = "type", values_to = "overdispersion")

# plot one one on one comaprison (these are usually not very informative)
ggplot(df.default, aes(hek, pbmc)) + geom_point() + geom_abline(color="red")

# plot summary
ggplot(df.long, aes(overdispersion, type, fill=type)) + stat_density_ridges(quantile_lines = T, quantiles = 2)  + scale_x_log10() + scale_fill_brewer(type = "qual", palette = "Set2")
