
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
library(openxlsx, quietly = TRUE)
library(topGO, quietly = TRUE)

#Setting seed
set.seed(777)

###Setting script settings -- Turn on if not on MacOS device
options(future.globals.maxSize = 40 * 1024^3)
plan("multisession", workers = 12)


###Specify usage
message("Usage: Rscript SPG_ProcessingAnnotation.R <output directory (SAME AS IN SPG_Criteria.R)> <SPG_Parameters.yaml file (see README)>")

###Deal arguments
args = commandArgs(T)

outdir = args[1]
parameters = args[2]

###Loading parameters
handlers = list("bool#no" = function(x){if ( x %in% c("false", "FALSE") ) FALSE else x}, "bool#yes" = function(x){if ( x %in% c("true", "TRUE") ) TRUE else x})
parms = yaml::yaml.load_file(parameters, handlers = handlers)

#Create output directory for informative documents
annotateDir = file.path(outdir, "Annotation")
if (!dir.exists(annotateDir)) {
  dir.create(annotateDir, recursive = TRUE)
}

###Create directory for plots made in this document
#This is the final plotting directory, which contains the plots from the subsequent SPG_Plotting.R script
#For easier comparison I already make the directory here so pre-integration plots can be stored in the same directory as the post-integration plots
plottingDir = file.path(outdir, "FinalPlots")
if (!dir.exists(plottingDir)) {
  dir.create(plottingDir, recursive = TRUE)
}

###Load preQC merged object
merged = readRDS(file.path(outdir, "SeuratOBJs/preQC_merged.rds"))


###Filter on criteria specified in parameters file
message("====== Filtering specified criteria ======")
indiFilter = (
  
  #Individual 1
  (merged$individual == "i1" &
     merged$nFeature_RNA > parms$filter$individual1$MINnFeature_RNA &
     merged$nFeature_RNA < parms$filter$individual1$MAXnFeature_RNA &
     merged$percent.mt < parms$filter$individual1$MAXpercent.mito
  ) |
  
  #Individual 2
  (merged$individual == "i2" &
     merged$nFeature_RNA > parms$filter$individual2$MINnFeature_RNA &
     merged$nFeature_RNA < parms$filter$individual2$MAXnFeature_RNA &
     merged$percent.mt < parms$filter$individual2$MAXpercent.mito
  ) |
    
  #Individual 3
  (merged$individual == "i3" &
     merged$nFeature_RNA > parms$filter$individual3$MINnFeature_RNA &
     merged$nFeature_RNA < parms$filter$individual3$MAXnFeature_RNA &
     merged$percent.mt < parms$filter$individual3$MAXpercent.mito
  ) |
    
  #Individual 4
  (merged$individual == "i4" &
     merged$nFeature_RNA > parms$filter$individual4$MINnFeature_RNA &
     merged$nFeature_RNA < parms$filter$individual4$MAXnFeature_RNA &
     merged$percent.mt < parms$filter$individual4$MAXpercent.mito
  ) |
  
  #Embryo 1
  (merged$individual == "e1" &
      merged$nFeature_RNA > parms$filter$embryo1$MINnFeature_RNA &
      merged$nFeature_RNA < parms$filter$embryo1$MAXnFeature_RNA &
      merged$percent.mt < parms$filter$embryo1$MAXpercent.mito
  ) |
    
    #Embryo 2
    (merged$individual == "e2" &
       merged$nFeature_RNA > parms$filter$embryo2$MINnFeature_RNA &
       merged$nFeature_RNA < parms$filter$embryo2$MAXnFeature_RNA &
       merged$percent.mt < parms$filter$embryo2$MAXpercent.mito
    ) |
    
    #Embryo 3
    (merged$individual == "e3" &
       merged$nFeature_RNA > parms$filter$embryo3$MINnFeature_RNA &
       merged$nFeature_RNA < parms$filter$embryo3$MAXnFeature_RNA &
       merged$percent.mt < parms$filter$embryo3$MAXpercent.mito
    ) |
    
    #Embryo 4
    (merged$individual == "e4" &
       merged$nFeature_RNA > parms$filter$embryo4$MINnFeature_RNA &
       merged$nFeature_RNA < parms$filter$embryo4$MAXnFeature_RNA &
       merged$percent.mt < parms$filter$embryo4$MAXpercent.mito
    )

)

#Filtering
filtered = subset(merged, cells = colnames(merged)[indiFilter])

#Storing what cells were removed by filtering
filteredCells = colnames(merged)[!indiFilter]
fcDF = data.frame(filteredCells)

write.xlsx(fcDF, file = file.path(annotateDir, "filtered.cells.xlsx"), sheetName = "Cells filtered out")


message("====== Doublet detection and removal ======")
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes


#Merged object is nicely filtered and contains all samples, need to detect per sample
doubs = SplitObject(filtered, split.by = "sample")


#Preprocess data, required for doublet detection
doubletDir = file.path(outdir, "DoubletPlots")
if (!dir.exists(doubletDir)) {
  dir.create(doubletDir, recursive = TRUE)
}

doubs = lapply(doubs, function(x) {
  
  x %>%
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE) %>%
    ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% 
    RunPCA() %>%
    RunUMAP(dims = 1:15, reduction = "pca", n.neighbors = 30, 
          n.components = 2L, umap.method = "uwot") #Settings copied from Bo
  
})

for(i in seq_along(doubs)) {
  
  #Set pN
  pN = 0.25
  
  #Get pK -- Adjusted very slightly during troubleshooting
  dims = 1:20
  sweep.res.list = paramSweep(doubs[[i]], PCs = dims, sct = FALSE)
  sweep.stats = summarizeSweep(sweep.res.list)
  bcmvn = find.pK(as.data.frame(sweep.stats))
  bcmetric = bcmvn[["BCmetric"]]
  pkvals = bcmvn[["pK"]]
  pK = as.numeric(as.character(pkvals[which.max(bcmetric)])) 
  
  #Get nExp
  rate = 7.6 * 10^-6 * ncol(doubs[[i]]) + 5.27 * 10^-4
  nExp_poi = round(as.numeric(rate) * ncol(doubs[[i]])) 
  
  #Identify doublets
  doubs[[i]] = doubletFinder(
    doubs[[i]], 
    PCs = dims, 
    pN = pN, 
    pK = pK, 
    nExp = nExp_poi, 
    reuse.pANN = NULL, #This on FALSE will result in an error, needs to be NULL
    sct = FALSE)
  colnames(doubs[[i]]@meta.data)[grep('pANN', colnames(doubs[[i]]@meta.data))] = "pANN"
  colnames(doubs[[i]]@meta.data)[grep('DF.classifications', colnames(doubs[[i]]@meta.data))] = "classifications"
  
  d = DimPlot(
    doubs[[i]], 
    reduction = "pca", 
    group.by = "classifications", 
    cols = c("Singlet" = "black", "Doublet" = "red")) + 
    ggtitle(names(doubs)[[i]])
  
  u = DimPlot(
    doubs[[i]],
    group.by = "classifications", 
    cols = c("Singlet" = "black", "Doublet" = "red")) +	
    ggtitle(names(doubs)[[i]])

  ggsave(file.path(doubletDir, paste0(names(doubs)[[i]], ".doubletPlot.pdf")), d, width = 10, height = 7)
  ggsave(file.path(doubletDir, paste0(names(doubs)[[i]], ".doubletUMAP.pdf")), u, width = 10, height = 7)
}


#Remove doublets from each sample
sampDoublets = unlist(lapply(doubs, function(x) {
  colnames(x)[x$classifications == "Doublet"]
}))

length(sampDoublets) #Check number of doublets
filteredSinglets = subset(filtered, cells = setdiff(colnames(filtered), sampDoublets))

#Storing what cells were removed
doubIdent = filtered$orig.ident[colnames(filtered) %in% sampDoublets]
doubletDF = data.frame(sampDoublets, doubIdent)

write.xlsx(doubletDF, file = file.path(doubletDir, "removed.doublets.xlsx"), sheetName = "Removed doublet IDs")  #Store doublet identities



message("====== Separating data by developmental stage to avoid incorrect parameter estimation ======")
devdata = SplitObject(filteredSinglets, split.by = "devstage")

#Defining function to do all subsequent steps per developmental stage dataset
processing = function(obj, devstage) {
  
  prefix = paste0(devstage, "_")

  message("====== Data normalization, cellcycle scoring, scaling, and dim reduction ======")
  processed = obj %>%
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) %>%
    ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% 
    RunPCA() 

  ####################################### Section added to create UMAPs with pre-integration to assess integration effect
  preint = processed %>%
    RunUMAP(dims = 1:15, reduction = "pca", n.neighbors = 30, 
            n.components = 2L, umap.method = "uwot") #Settings copied from Bo

  for (i in 1) {
  
    #Make plots per element to assess whether integration worked
    sample = DimPlot(preint, group.by = "sample", alpha = 0.3, raster = FALSE) 
    tissue = DimPlot(preint, group.by = "tissue", alpha = 0.3, raster = FALSE)
    protocol = DimPlot(preint, group.by = "protocol", alpha = 0.3, raster = FALSE)
    individual = DimPlot(preint, group.by = "individual", alpha = 0.3, raster = FALSE)
    cycle = DimPlot(preint, group.by = "Phase", alpha = 0.3, raster = FALSE) #Added this to assess whether cell cycle stage was properly regressed out also, not strictly integration assessment

    ggsave(file.path(plottingDir, paste0(prefix, "pre-integration.sample_UMAP.pdf")), sample, width = 24, height = 16)
    ggsave(file.path(plottingDir, paste0(prefix, "pre-integration.tissue_UMAP.pdf")), tissue, width = 24, height = 16)
    ggsave(file.path(plottingDir, paste0(prefix, "pre-integration.protocol_UMAP.pdf")), protocol, width = 24, height = 16)
    ggsave(file.path(plottingDir, paste0(prefix, "pre-integration.individual_UMAP.pdf")), individual, width = 24, height = 16)
    ggsave(file.path(plottingDir, paste0(prefix,  "pre-integration.cell_cycle_stage_UMAP.pdf")), cycle, width = 24, height = 16)

  }
#######################################

  message("====== Performing data integration ======")
  integrated = RunHarmony(processed, "orig.ident", plot_convergence = FALSE)


  message("====== Identifying clusters ======")
  integratedElbow = ElbowPlot(integrated, reduction = "harmony")
  criteriaDir = file.path(outdir, "CriteriaPlots")
  ggsave(file.path(criteriaDir, paste0(prefix, "postHarmony_Elbowplot.pdf")), integratedElbow, width = 10, height = 7)

  clustered = integrated %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>% #Not sure whether to base dims on harmony or pca ElbowPlot
    FindClusters(resolution = 0.5, random.seed = 777) #Adjusted back to default


  message("====== Running UMAP ======")
  UM = clustered %>%
    RunUMAP(dims = 1:15, reduction = "harmony", n.neighbors = 30, 
            n.components = 2L, umap.method = "uwot") #Unsure, copied settings from Bo, but kept dims the same as in FindNeighbours based on Cecile

  ##Checking cluster characteristics

  #Number of cells per cluster
  length(levels(Idents(UM))) #Check number of clusters explicitly
  cluschar = UM@meta.data %>%
    dplyr::count(seurat_clusters, sample)

  clustCells = ggplot(cluschar) +
    geom_col(aes(x = seurat_clusters, y = n, fill = sample)) + 
    scale_fill_manual(values = c(
      "#925E9F", "#76C8DC", "#FDAF91", "#FF7777", "#AD002A",
      "#42C1BB", "#0F8074", "#0A7C2E", "#B8D24D", "#FAB158", #Colours there for when 5 samples become 15
      "#792244", "#AD556B", "#B09F91", "#DEB8A1", "#4C4E4E" 
    )) +
    labs(
      x = "Cluster",
      y = "Number of cells"
    ) +
    theme_classic() + 
    RotatedAxis()

  ggsave(file.path(annotateDir, paste0(prefix, "Sample_Cluster_cellDistribution.pdf")), clustCells, width = 10, height = 7)


  #Which cells in which cluster
  clusterStat = data.frame(UM$seurat_clusters, UM$tissue, UM$sample)
  colnames(clusterStat) = c("Cluster", "Tissue", "Sample")
  clusterStat = clusterStat %>%
    rownames_to_column("Cell ID")

  write.xlsx(clusterStat, file = file.path(annotateDir, paste0(prefix, "Cells.cluster.list.xlsx")), sheetName = "Cell + cluster combination")

  message("====== Finding marker genes ======")
  markers = UM %>%   #First need to do JoinLayers; FindAllMarkers does not work otherwise
    JoinLayers() %>%
    FindAllMarkers(logfc.threshold = 0.5, min.pct = 0.4, only.pos = T, test.use = 'MAST')


  #Number of DE genes per cluster
  clusDE = markers %>%
    dplyr::count(cluster, gene)

  clusDEstat = clusDE %>%
    group_by(cluster) %>%
    summarise(nDE = length(cluster))
  colnames(clusDEstat) = c("Cluster", "Number of DE genes")

  DEclust = ggplot(clusDE) +
    geom_col(aes(x = cluster, y = n, colour = cluster)) +
    labs(
      x = "Cluster",
      y = "Number of DE genes"
    ) +
    theme_classic() +
    theme(legend.position = "none")

  ggsave(file.path(annotateDir, paste0(prefix, "DEgene.count.per.cluster.pdf")), DEclust, width = 10, height = 7)
  write.xlsx(clusDEstat, file = file.path(annotateDir, paste0(prefix, "DEgene.count.per.cluster.xlsx")), sheetName = "DE gene counts")



  message("====== Saving cluster marker object and final Seurat object ======")
  seuratDir = file.path(outdir, "SeuratOBJs")
  saveRDS(markers, file.path(seuratDir, paste0(prefix, "markers.rds")))
  saveRDS(UM, file.path(seuratDir, paste0(prefix, "UM.rds")))


  message("====== Creating canonical marker dotplots and DE gene list ======")
  #Specifying canonical markers for different cell types
  clusters = levels(markers$cluster)

  myofibroblast = c("MYF5","NRXN1", "FGFR4", "MYO18A", "MYLK", "MYL9") 
  myocyte = c("MYH11", "MYLK", "MYL4") 
  keratinocyte = c("WNT6", "WNT4", "KRT14", "KRT13", "KRT5", "PRLR") 
  fibroblast = c("APOA1", "COL1A1", "CXCL14", "COL1A2", "PDGFRA") 
  mesenchymal = c("BARX2", "FABP5", "COL17A1")
  neurocyte = c("SMC2", "TOP2A") 
  macrophage = c("LCP1", "RGS1") 
  melanocyte = c("MLANA", "PMEL", "EDNRB2") 
  immunecluster = c("TCF7", "CD24", "CD48", "MSX2") 
  endothelial = c("APOLD1", "VWF", "TMSB4X", "GNG11", "RHOA") 
  dendritic = c("C1QA", "AT2") 
  adipocyte = c("ELOVL4", "ELOVL7", "ADIPOQ")
  pericyte = c("RGS4", "ACTA2", "RGS5")
  tcell = c("TCF7", "GJA1")


  cellgenes = list(myofibroblast, myocyte, keratinocyte, fibroblast, mesenchymal, neurocyte, 
                       macrophage, melanocyte, immunecluster, endothelial, dendritic, adipocyte,
                       pericyte, tcell)
  names(cellgenes) = c("myofibroblast", "myocyte", "keratinocyte", "fibroblast", "mesenchymal", "neurocyte", 
                     "macrophage", "melanocyte", "immunecluster", "endothelial", "dendritic", "adipocyte",
                     "pericyte", "tcell")

  #Creating dotplots per type
  for (type in names(cellgenes)) {
    
    d = DotPlot(UM, features = cellgenes[[type]]) + 
      scale_color_gradient(
        low = "black", 
        high = "red") +
      labs(
        title = type,
        y = "Cluster",
        x = "Gene"
      ) +
      guides(color = guide_colorbar(title = "Z-scaled expression (log2FC)"))

    ggsave(file.path(annotateDir, paste0(prefix, type, ".dotPlot.pdf")), d, width = 10, height = 7)

  }


  #Saving DE genes per clusters for future annotation
  topDE = markers %>%
    arrange(cluster, desc(avg_log2FC)) #Sort so highest log2FC is at the top

  genelist = split(topDE$gene, topDE$cluster) #Create dataframe 
  maxLength = max(lengths(genelist))

  genedf = as.data.frame(
    lapply(genelist, function (x) {
      c(x, rep("", maxLength - length(x))) #Avoids issues with different lengths
    }),
    check.names = FALSE
  )

  write.xlsx(list(
    "Marker genes - log2FC + p" = markers,
    "Top DE genes (log2FC desc)" = genedf), 
    file = file.path(annotateDir, paste0(prefix, "DEgeneList.xlsx")))


  #Getting average expression per cluster for preliminary overview (prior to successful annotation)
  annAvg = as.data.frame(AverageExpression(UM, group.by = "seurat_clusters")$RNA)
  colnames(annAvg) = sub("g", "Cluster ", colnames(annAvg))
  annAvg = annAvg %>%
    rownames_to_column("Gene ID")

  write.xlsx(annAvg, file = file.path(annotateDir, paste0(prefix, "AllGenes.AvgExp.cluster.xlsx")), sheetName = "Average expression per cluster")


  #Getting percentage of cells expressing the gene in that cluster
  UMstat = JoinLayers(UM[["RNA"]])
  statmat = GetAssayData(UMstat, layer = "data")
  cbc = split(colnames(UM), Idents(UM))

  fracCells = sapply(cbc, function (x) {
    
    (Matrix::rowSums(statmat[, x] > 0) / length(x)) * 100
    
  })
  colnames(fracCells) = paste("Cluster", colnames(fracCells), sep = " ")
  pctExp = as.data.frame(fracCells) %>%
    rownames_to_column("Gene ID")

  write.xlsx(pctExp, file = file.path(annotateDir, paste0(prefix, "AllGenes.PctCellsExpr.cluster.xlsx")), sheetName = "% of cells expressing gene")



  #Storing metadata
  xlMeta = UM@meta.data %>%
    rownames_to_column("Sample-Cell_barcode")

  write.xlsx(xlMeta, file = file.path(annotateDir, paste0(prefix, "metadata.xlsx")), sheetName = "Cell Seurat metadata")


  ###Finished processing
  message("<><><><><><><><><><> Processing done <><><><><><><><><><>")



if(isTRUE(parms$printer$list)) {
  
  message("::: Processing finished, this script will now print top 25 genes based on average log2FC and adjusted p value for easy ChatGPT annotation, as well as cluster GO terms based on these top genes :::")
  message("<><><> Currently printing for ", devstage, " <><><>")
  Sys.sleep(10) #Time to read message
  
  
  clusters = levels(markers$cluster)
  annotation = read.delim(parms$topGO_annotation) %>% #Genes and their GO association for topGO
    mutate(annoID = if_else(
      Gene.name != "",
      Gene.name,
      Gene.stable.ID
    )) %>%
    filter(GO.term.accession != "") %>%
    group_by(annoID) %>%
    summarise(GO = list(unique(GO.term.accession)),
              .groups = "drop") %>%
    tibble::deframe()
  anno_universe = rownames(UM) #Gene universe for topGO
  GO_LIST = list() #Empty GO vector for final printer
  P_LIST = list() #Empty vector for top 25 genes based on P value
  L_LIST = list() #Empty vector for top 25 genes based on log2FC
  
  #Mitochondrial genes to remove for annotation
  mito_symbols = c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3",
                 "ND3","ND4L","ND4","ND5","CYTB","ND6")

  #Top 25 log2FC genes (ONLY NON-ENSGALGS BECAUSE THOSE ARE NOT INFORMATIVE)
  logtop = markers %>%
    filter(
     !grepl("^RP[LS]", gene), #Filter out housekeeping genes
     !grepl("^ENSGAL", gene), #Filter out unknown genes
     !gene %in% mito_symbols) %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    top_n(n = 25, wt = avg_log2FC)
  
  #Top 25 pvaladj genes - can be more since many share 0 
  pvaltop = markers %>%
    filter( 
     !grepl("^RP[LS]", gene), #Filter out housekeeping genes
     !gene %in%	mito_symbols) %>%
    group_by(cluster) %>%
    top_n(n = -25, wt = p_val_adj) #- is to get smallest values
  
  #Get list per cluster
  for(i in seq_along(clusters)) {
    
    cluster = clusters[i]
    cat("Currently working on cluster", cluster) #Progress marker
    Lgenes = logtop$gene[as.character(logtop$cluster) == cluster] #Fold change top genes
    Pgenes = pvaltop$gene[as.character(pvaltop$cluster) == cluster] #P value top genes
    L_LIST[[i]] = Lgenes
    P_LIST[[i]] = Pgenes
    
    combgenes = union(Lgenes, Pgenes) #For GO term annotation
    
    #GO term annotation
    anno_list = factor(as.integer(anno_universe %in% combgenes))
    names(anno_list) = anno_universe
    
    clusGO = new("topGOdata",ontology = "BP",allGenes = anno_list,annot = annFUN.gene2GO,gene2GO = annotation)
    anno_result = runTest(clusGO, algorithm = "weight01", statistic = "fisher")
    
    GOtable = GenTable(
      clusGO, weightFisher = anno_result, orderBy = "weightFisher", topNodes = 25
    )
    GO_LIST[[i]] = GOtable
    
    

  }
  
  for(i in seq_along(clusters)) {
    
    cluster = clusters[i]
    cat(" ==============", "\n", "Cluster", cluster, "\n", "==============","\n","top DE genes by avg.log2FC are:", L_LIST[[i]], "\n \n",  
        "this cluster is associated with the following GO terms (top 25 based on Fisher weight):", GO_LIST[[i]]$Term, "\n \n")
    
  }
  
}

}


#Applying the function to the different developmental stages
for(devstage in names(devdata)) {
  processing(devdata[[devstage]], devstage)
}




message("..... Complete .....")