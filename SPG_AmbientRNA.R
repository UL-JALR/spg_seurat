
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
library(glmGamPoi, quietly = TRUE)
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
library(openxlsx, quietly = TRUE)
library(topGO, quietly = TRUE)
library(ggrepel, quietly = TRUE)
library(SoupX, quietly = TRUE)
library(DropletUtils, quietly = TRUE)


#Setting seed
set.seed(777)

###Specify usage
message("Usage: Rscript SPG_AmbientRNA.R <raw matrix directory> <filtered matrix directory>")
#Sys.sleep(15) #Pause to allow user to read warning

###Deal arguments
args = commandArgs(T)

rawdir = args[1]
fildir = args[2]


filList = list.files(fildir)  #Create a list containing all tissue directories


message("====== Reading in filtered feature bc matrices ======")
sampListALL = list()
for(i in filList) {
  sampListALL[[i]] = Read10X(file.path(fildir, i))
}
names(sampListALL) = sub("\\.filtered_feature_bc_matrix", "", names(sampListALL)) #Rename for ease

message("====== Seurat object creation ======")
SeuALL = lapply(names(sampListALL), function(j) {
  counts = sampListALL[[j]]
  seurats = CreateSeuratObject(counts = counts, project = j)
  seurats$sample = j
  seurats
})
names(SeuALL) = names(sampListALL)


message("====== Run all standard processing steps (NOT integration) per sample (we want ambient RNA correction per sample) ======")
clusALL = list()
for(i in names(SeuALL)) {
  
  clusALL[[i]] = SeuALL[[i]] %>%
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:15) %>%
    FindClusters(resolution = 0.5, random.seed = 777)
  
}

#Explicitly check number of clusters created
sapply(clusALL, \(x) length(unique(Idents(x))))


message("====== Running SoupX adjustment, loading both filtered and raw, requirement for SoupX ======")
rawList = list.files(rawdir)  #Create a list containing all tissue directories
rawListALL = list()
contamination = data.frame()
for(i in rawList) {
  
  message("Currently processing ", i)
  
  rawListALL[[i]] = load10X(file.path(rawdir, i))
  rawListALL[[i]] = setClusters(rawListALL[[i]], Idents(clusALL[[i]]))
  rawListALL[[i]] = autoEstCont(rawListALL[[i]], doPlot = TRUE)
  
  contamination = rbind(contamination, data.frame(sample = i, cont_fraction = rawListALL[[i]]$fit$rhoEst))
  
}


#Actually producing the adjusted data
fixDir = file.path(getwd(), "SoupMatrices")
if (!dir.exists(fixDir)) {
  dir.create(fixDir, recursive = TRUE)
}

contaminatedDir = file.path(getwd(), "ContaminationFracSamples")
if (!dir.exists(contaminatedDir)) {
  dir.create(contaminatedDir, recursive = TRUE)
}

correctedALL = list()

for(i in rawList) {
  
  message("Currently adjusting counts for ", i)
  correctedALL[[i]] = adjustCounts(rawListALL[[i]])
  
}

for(i in names(correctedALL)) {
  
  write10xCounts(
    path = file.path(fixDir, paste0(i, ".fixed_feature_bc_matrix")),
    x = correctedALL[[i]],
    overwrite = TRUE
  )
  
}


write.xlsx(contamination, file.path(contaminatedDir, "contaminated_fraction_per_sample.xlsx"), overwrite = TRUE)
