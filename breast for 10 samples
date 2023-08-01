library(dplyr)
library(Seurat)
library(patchwork)
library(reticulate)
library(SeuratWrappers)
library(sceasy)
library(ggplot2)
library(harmony)
library(Rcpp)
library(RCurl)
library(cowplot)
options(Seurat.object.assay.version = 'v5')
# install glmGamPoi
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("glmGamPoi")
# install sctransform from Github
install.packages("sctransform")
library(sctransform)

setwd("/disk2/user/merelm/breastsinglecell/GSE235326_RAW/")

read_allsinglecelldata <- lapply(Sys.glob("*filtered*", dirmark = FALSE)[seq(1,10)], Read10X_h5)
#in order to read all files which has "filtered" word in it.
head(read_allsinglecelldata)
seurat_object <- lapply(read_allsinglecelldata, CreateSeuratObject)
names_file <-Sys.glob("*filtered*", dirmark = FALSE)[seq(1,10)]
for (i in 1:length(seurat_object)) {
  seurat_object[[i]]$orig.ident <- names_file[[i]]
}
big_merge <- merge(seurat_object[[1]], y = seurat_object[seq(2,length(seurat_object))], merge.data = TRUE)
#normalized_data <- NormalizeData(big_merge)
big_merge <- SCTransform(big_merge)

#top2000 <- head(VariableFeatures(normalized_data), 2000)
#normalized_data <- normalized_data[top2000]

#integ_features <- SelectIntegrationFeatures(object.list = norm_seurat_list, nfeatures = 3000) 


big_merge<- RunPCA(big_merge,npcs = 30, verbose = FALSE)

#digestion methods were not taken into account.
big_merge <- RunHarmony(big_merge,
                                  group.by.vars = c("orig.ident"),
                                  reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

big_merge<- RunUMAP(big_merge, reduction = "harmony", assay = "SCT", dims = 1:30, verbose = TRUE)


big_merge <- RunUMAP(big_merge, reduction = "harmony",  dims = 1:30,  n.neighbors = "nneigh", verbose = TRUE )
big_merge <- FindNeighbors(big_merge, reduction = "pca", dims = 1:30, nneigh = 50, verbose = TRUE)
big_merge <- FindClusters(big_merge, reduction = "pca", dims = 1:30, resolution = 1, algorithm = 1, verbose = TRUE)
DimPlot(big_merge)
FeaturePlot(big_merge, "PTPRC")


#head(normalized_data)
#normalized_data <- IntegrateLayers(
 # object = normalized_data, method = scVIIntegration,
  #new.reduction = "integrated.scvi",
  #conda_env = "/disk2/user/merelm/.local/share/r-miniconda/envs/scvi", verbose = FALSE
#)


#scVIIntegration
