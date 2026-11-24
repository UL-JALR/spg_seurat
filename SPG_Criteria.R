
#Code by JALR

###Load potentially relevant packages
library(Seurat, quietly = TRUE)
library(Rcpp, quietly = TRUE)
library(harmony, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(patchwork, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(BPCells, quietly = TRUE)
library(presto, quietly = TRUE)
#library(glmGamPoi, quietly = TRUE)
library(hdf5r, quietly = TRUE)
library(plotly, quietly = TRUE)
library(plot3D, quietly = TRUE)
library(colorspace, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(gdata, quietly = TRUE)
library(pheatmap, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(DoubletFinder, quietly = TRUE)
library(future, quietly = TRUE)


#Setting seed
set.seed(777)

###Setting script settings
plan("multicore", workers = 6)


###Specify usage
message("Usage: Rscript SPG_Criteria.R <fixed data directory (ambient RNA corrected, see README)> <output directory>")

###Deal arguments
args = commandArgs(T)

datadir = args[1]
outdir = args[2]


tmpList = list.files(datadir)  #Create a list containing all tissue directories

message("====== Reading in filtered feature bc matrices ======")
sampListALL = list()
for(i in tmpList) {
  sampListALL[[i]] = Read10X(file.path(datadir, i))
}
names(sampListALL) = sub("\\.fixed_feature_bc_matrix", "", names(sampListALL)) #Rename for ease


message("====== Quick data structure check ======")

message("--- Data dimensions ---")
dims = sapply(sampListALL, dim) 
rownames(dims) = c("Number of genes", "Number of cells")
dims

message("--- Data structure ---")
views = sapply(sampListALL, function(x) x[1:5, 1:5]) 
views

message("--- # of genes of interest check (HOX) ---")
sums = sapply(sampListALL, function(z) rowSums(z[grep(rownames(z),pattern="^HOX"),]))
nrow(sums)


message("====== Seurat object creation ======")
SeuALL = lapply(names(sampListALL), function(j) {
  counts = sampListALL[[j]]
  seurats = CreateSeuratObject(counts = counts, project = j)
  seurats$sample = j
  seurats
})
names(SeuALL) = names(sampListALL)


message("====== Adding metadata ======")
#Flag individual objects by experimental individual and tissue (Specified in KRT)
for(i in seq_along(SeuALL)) {

  if (grepl("S1Neck|S1LowBack", names(SeuALL[i]))) {
    SeuALL[[i]]@meta.data$individual = "i1"
    SeuALL[[i]]@meta.data$protocol = "scRNAseq"
    SeuALL[[i]]@meta.data$devstage = "adult"
  }

  if (grepl("S1Thigh|S1Bicep|S1Uropygium", names(SeuALL[i]))) {
    SeuALL[[i]]@meta.data$individual = "i2"
    SeuALL[[i]]@meta.data$protocol = "snRNAseq"
    SeuALL[[i]]@meta.data$devstage = "adult"

  }

  if (grepl("S2", names(SeuALL[i]))) {
    SeuALL[[i]]@meta.data$individual = "i3"
    SeuALL[[i]]@meta.data$protocol = "snRNAseq"
    SeuALL[[i]]@meta.data$devstage = "adult"
  }

  if (grepl("S3", names(SeuALL[i]))) {
    SeuALL[[i]]@meta.data$individual = "i4"
    SeuALL[[i]]@meta.data$protocol = "snRNAseq"
    SeuALL[[i]]@meta.data$devstage = "adult"
  }
  
  if (grepl("EmbUpper", names(SeuALL[i]))) {
    SeuALL[[i]]@meta.data$individual = "e1"
    SeuALL[[i]]@meta.data$protocol = "scRNAseq"
    SeuALL[[i]]@meta.data$devstage = "embryo"
  }
  
  if (grepl("EmbLower", names(SeuALL[i]))) {
    SeuALL[[i]]@meta.data$individual = "e2"
    SeuALL[[i]]@meta.data$protocol = "scRNAseq"
    SeuALL[[i]]@meta.data$devstage = "embryo"
  }
  
  if (grepl("EmbFore", names(SeuALL[i]))) {
    SeuALL[[i]]@meta.data$individual = "e3"
    SeuALL[[i]]@meta.data$protocol = "scRNAseq"
    SeuALL[[i]]@meta.data$devstage = "embryo"
  }
  
  if (grepl("EmbHind", names(SeuALL[i]))) {
    SeuALL[[i]]@meta.data$individual = "e4"
    SeuALL[[i]]@meta.data$protocol = "scRNAseq"
    SeuALL[[i]]@meta.data$devstage = "embryo"
  }

  SeuALL[[i]]@meta.data$tissue = sub("^..", "", SeuALL[[i]]@meta.data$sample)

}


#Define mitochondrial genes and determine mitochondrial percentage in each sample
mito_symbols = c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3",
                 "ND3","ND4L","ND4","ND5","CYTB","ND6")

for(i in seq_along(SeuALL)) {
  
  mito_genes = intersect(mito_symbols, rownames(SeuALL[[i]]))
  SeuALL[[i]]@meta.data$percent.mt = PercentageFeatureSet(SeuALL[[i]], features = mito_genes)
  
}

message("====== Merging into a single object and saving it ======")
#Seurat merge
merged = merge(SeuALL[[1]], SeuALL[-1], add.cell.ids = names(SeuALL))

#Saving object
seuratDir = file.path(outdir, "SeuratOBJs")
if (!dir.exists(seuratDir)) {
  dir.create(seuratDir, recursive = TRUE)
}

saveRDS(merged, file.path(seuratDir, "preQC_merged.rds"))

message("====== Merged Seurat object successfully created and saved ======")


message("====== Creating criteria selection figures ======")
#Feature violin plots
v = VlnPlot(merged, 
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
            group.by = "individual",
            ncol = 3)

#Feature scatterplots and save them (group.by individual was very unclear for these)
scatters = vector("list", length(unique(merged$individual)))
names(scatters) = unique(merged$individual)


for(i in unique(merged$individual)){
  
  indiv = subset(merged, subset = individual == i)
  s1 = FeatureScatter(indiv, 
                      feature1 = "nCount_RNA", 
                      feature2 = "nFeature_RNA",
                      pt.size = 0.2) +
    ggtitle(i)
  s2 = FeatureScatter(indiv, 
                      feature1 = "nCount_RNA", 
                      feature2 = "percent.mt",
                      pt.size = 0.2) +
    ggtitle(i)
  
  s3 = FeatureScatter(indiv, 
                      feature1 = "nFeature_RNA", 
                      feature2 = "percent.mt",
                      pt.size = 0.2) +
    ggtitle(i)
  
  
  scatters[[i]] = s1 + s2 + s3

}


#Save plots
criteriaDir = file.path(outdir, "CriteriaPlots")
if (!dir.exists(criteriaDir)) {
  dir.create(criteriaDir, recursive = TRUE)
}


for(i in seq_along(1:length(scatters))) {
  
  ggsave(file.path(criteriaDir, paste0(scatters[[i]]@labels$title, "_scatter.pdf")), scatters[[i]], width = 10, height = 7)
  
}

ggsave(file.path(criteriaDir, "Violin.pdf"), v, width = 10, height = 7)



###Finished
message("<><><><><><><><><><> All done <><><><><><><><><><>")
