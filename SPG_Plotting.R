
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
library(reshape2, quietly = TRUE)
library(topGO, quietly = TRUE)
library(ggrepel, quietly = TRUE)
library(future, quietly = TRUE)
library(vegan, quietly = TRUE)
library(openxlsx, quietly = TRUE)
library(lme4, quietly = TRUE)
library(lmerTest, quietly = TRUE)

#Setting seed
set.seed(777)



###Setting script settings -- Turn on if not on MacOS device
#plan("multicore", workers = 12)


###Specify usage
message("Usage: Rscript SPG_Plotting.R <output directory (SAME AS IN SPG_Criteria.R)> <SPG_Parameters.yaml file (see README)>")
message("WARNING: Note that the number of annotated clusters in SPG_Parameters.yaml must match the actual number of clusters, see output/Annotation/DEgeneList.xlsx")
#Sys.sleep(15) #Pause to allow user to read warning if running live

###Deal arguments
#args = commandArgs(T)

outdir = "../output" #args[1]
parameters = "SPG_Parameters.yaml" #args[2]

###Loading parameters
handlers = list("bool#no" = function(x){if ( x %in% c("false", "FALSE") ) FALSE else x}, "bool#yes" = function(x){if ( x %in% c("true", "TRUE") ) TRUE else x})
parms = yaml::yaml.load_file(parameters, handlers = handlers)

###Create directories for plots made in this document (main directory should already exist from SPG_ProcessingAnnotation.R)
plottingDir = file.path(outdir, "FinalPlots")
if (!dir.exists(plottingDir)) {
  dir.create(plottingDir, recursive = TRUE)}

hoxDir = file.path(plottingDir, "HOX gene expression")
if (!dir.exists(hoxDir)) {
  dir.create(hoxDir, recursive = TRUE)}

###Load list of spatial patterning genes
spgfile = read.xlsx(parms$SPGlist, sheet = "SPG_Master")
spgs = toupper(spgfile$Catalogue_gene_symbol) #Not all genes are capitalized but chicken genes are, so need to specify toupper

###Defining a function for heatmap scaling
safeHeat = function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  minimum = min(x, na.rm = TRUE)
  maximum = max(x, na.rm = TRUE)
  
  if(maximum == minimum) {
    return(rep(0, length(x)))
  }
  (x - minimum) / (maximum - minimum)
}

###Defining full devstage specific plotting function to be applied to both embryo and adult data
#plotfunc = function(devstage) {
  
  #Define prefix
  prefix = paste0(devstage, "_")
  
  
  #Load devstage files
  UM = readRDS(file.path(outdir, "SeuratOBJs", paste0(prefix, "UM.rds")))
  markers = readRDS(file.path(outdir, "SeuratOBJs", paste0(prefix, "markers.rds")))
  
  
  #Setting cluster-cell type annotations
  message("====== Setting cluster annotations for ", devstage, " from SPG_Parameters.yaml ======")
  anno = parms[[paste0(prefix, "annotation")]]
  annLabels = sapply(names(anno), function(x) { #Create naming vector
    clusterN = sub("^clus", "", x)
    paste(clusterN, anno[[x]]$name)})
  
  annLookup = setNames(annLabels, sub("^clus", "", names(anno))) #Create lookup vector
  
  print(parms[[paste0(prefix, "annotation")]]) #Check because changes to parameter file were not showing up for all plots
  
  
  
  ###Plotting from this point onwards
  
  
  
  ########### HEATMAP
  message("====== Creating Z-scaled HOX heatmap for ", devstage," ======")
  
  #Determining average expression
  UM[["celltype"]] = Idents(UM)
  avg = AverageExpression(UM, group.by = "celltype")
  
  #Creating object to make plotting easier
  TrueAvg = avg$RNA
  colnames(TrueAvg) = annLabels
  
  #Set heatmap colours
  color = colorRampPalette(c("#00468B", "white", "#AD000A"))(50) #Gradient
  
  #Filtering for HOX genes and then clustering them by biological paralog groups
  a = TrueAvg[grep("^HOX", rownames(TrueAvg)), ]
  a = a %>%
    as.data.frame() %>%
    mutate(
      hoxN = as.integer(str_extract(rownames(a), "\\d+$")),
      hoxL = str_extract(rownames(a), "(?<=^HOX)[ABCD]")
    ) %>%
    arrange(hoxN, hoxL) %>%
    dplyr::select(-hoxN, -hoxL) %>% #dplyr added because I got errors, potentially because BioMart was loaded
    as.matrix()
  
  #Z-scaling by row (not done typically I feel but Bo indicated that he attempted this but did not recall how exactly. Scaling using the built in "scale = "row" returned values between -4 and 4, and we want them to remain between 0 and 1. I do not know if this is more accurate, but Bo indicated that I should do this, and afterwards my figure looked much more like his. ChatGPT suggested the below code/function to do so, and indicated that the following should be included in the paper if this visualization is used: "Values represent row-wise z-scores rescaled to the 0–1 range for visualization.")
  az = t(scale(t(a)))
  a1 = t(apply(az, 1, safeHeat))
  
  #Creating and saving heatmap
  heat = pheatmap(a1, 
                  color = color, 
                  cluster_rows = FALSE, 
                  cluster_cols = FALSE,
                  silent = TRUE)
  
  ggsave(file.path(plottingDir, paste0(prefix, "HOX_Heatmap.pdf")), heat, width = 6, height = 8)
  ###########
  
  
  
  ########### CORRELATION MATRIX
  message("====== Creating correlation matrix for ", devstage, " ======")
  corAvg = AverageExpression(UM)$RNA
  cormatrix = cor(as.matrix(corAvg))
  corOrd = colnames(corAvg)
  cormatrix = cormatrix[corOrd, corOrd]
  rownames(cormatrix) = sub("g", "Cluster ", rownames(cormatrix))
  colnames(cormatrix) = sub("g", "Cluster ", colnames(cormatrix))
  
  corDF = reshape2::melt(cormatrix) #Swapping to long format so I can set shape size and add value using ggplot since pheatmap can't
  colnames(corDF) = c("Cluster1", "Cluster2", "Correlation")
  
  corHeat = ggplot(corDF, aes(x = Cluster1, y = Cluster2)) +
    geom_point(size = 12, color = "grey90", shape = 0) +
    geom_point(aes(size = Correlation, fill = Correlation), shape = 22) +
    
    scale_fill_gradient2(low = "#00468B", mid = "white", high = "#AD002A", midpoint = 0.5) +
    scale_size(range = c(0,12)) +
    geom_text(aes(label = round(Correlation, 3)), size = 3, fontface = "bold") +
    labs(
      x = NULL,
      y = NULL,
      fill = NULL
    ) +
    guides(
      size = "none"
    ) +
    theme_classic() +
    RotatedAxis()
  
  ggsave(file.path(plottingDir, paste0(prefix, "Cluster_Correlations.pdf")), corHeat, width = 14, height = 12)
  ###########
  
  
  
  ########### UMAPs
  message("====== Creating UMAPs for ", devstage," ======")
  
  ##Quick check UMAPs
  tissueCols = c("#792244", "#5C298F", "#0F8074", "#FAB158", "#6C6DA4")
  paralogCols = c("#990000", "#C21807", "#E53917", "#F46D1A", 
                  "#D9EF8B", "#A6D96A", "#66BD63", "#1A9850", 
                  "#74ADD1", "#4575B4", "#2C7BB6", "#313695", 
                  "#1F2A7C")
  sampleCols = c(
    "#925E9F", "#76C8DC", "#FDAF91", "#FF7777", "#AD002A",
    "#42C1BB", "#0F8074", "#0A7C2E", "#B8D24D", "#FAB158", #Extra colours ready for expansion from 5 to 15 samples
    "#792244", "#AD556B", "#B09F91", "#DEB8A1", "#4C4E4E"
  )
  
  dimSamp = DimPlot(UM, reduction = "umap", group.by = "sample", raster = FALSE)
  dimTiss = DimPlot(UM, reduction = "umap", group.by = "tissue", cols = tissueCols, raster = FALSE)
  
  ggsave(file.path(plottingDir, paste0(prefix, "sample_UMAP.pdf")), dimSamp, width = 14, height = 10)
  ggsave(file.path(plottingDir, paste0(prefix, "tissue_UMAP.pdf")), dimTiss, width = 14, height = 10)
  
  ##Full UMAPs
  #Coordinate extraction for labels
  UMAP.coor = UM@reductions$umap@cell.embeddings
  UMAP.coor = as.data.frame(UMAP.coor)
  
  #Cell annotation, numbering and annotation separately to label with numbers and color with anno
  celltype = Idents(UM)
  UMAP.coor = cbind(UMAP.coor, celltype)
  annotation = annLookup[as.character(UMAP.coor$celltype)] 
  UMAP.coor = cbind(UMAP.coor, annotation)
  
  #Fix annotation so it will sort the legend like 0 1 2 3 and not 0 1 10 11
  labels = unique(UMAP.coor$annotation)
  getOrder = function(x) as.integer(sub("^\\s*(\\d+).*$","\\1", x)) #Captures number only
  FinalOrder = labels[order(getOrder(labels), labels)] #Sorts by the number
  UMAP.coor$annotation = factor(UMAP.coor$annotation, levels = FinalOrder)
  
  #Defining colours
  ClusColours = c(
    "#00468B", "#B8D24D", "#FD0000", "#5B4232",
    "#76C8DC", "#792244", "#FAB158", "#0A7C2E", "#5C298F",
    "#DEB8A1", "#0099B4", "#E67E74", "#5377A7", "#42C1BB",
    "#AD002A", "#EDE447", "#6C6DA4", "#28AA6C", "#FDAF91",
    "#3B81AB", "#AE8691", "#76D1B1", "#925E9F", "#FF7777",
    "#4C4E4E", "#759EDD", "#CE9573", "#0F8074", "#AD556B",
    "#B09F91"
  )
  
  #Defining labels
  label.coor = UMAP.coor %>% #Finding middle of UMAP clusters to put label
    group_by(celltype) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))
  
  #Creating and saving plot
  Ufull = ggplot(
    data = UMAP.coor, 
    aes(
      x = umap_1, 
      y = umap_2, 
      label = celltype)) +
    geom_point(
      aes(
        color = annotation), 
      size = 1, 
      alpha = 0.5) +
    geom_text(
      data = label.coor, 
              color = 'black', 
              size = 5, 
              fontface = "bold") +
    scale_color_manual(
      values = ClusColours) +
    labs(
      color = "Celltype") +
    theme_classic()
  
  ggsave(file.path(plottingDir, paste0(prefix, "annotated_UMAP.pdf")), Ufull, width = 14, height = 10)
  
  
  
  ###Also creating UMAPs per sample, several ways to be informative
  UM[["celltype"]] = unname(annLookup[as.character(UM$seurat_clusters)])
  for (i in unique(UM$sample)) {
    
    #Only sample cells
    subbed = subset(UM, subset = sample == i)
    empty = DimPlot(subbed, group.by = "sample", cols = "#6C6DA4", alpha = 0.3, raster = FALSE) 
    
    #Only sample cells separated into annotated clusters
    grouped = DimPlot(subbed, group.by = "celltype", label = TRUE, alpha = 0.3, raster = FALSE) 
    
    #This sample compared to all
    refcells = colnames(UM)[UM$sample == i] #Sample cells in the full UMAP
    reference = DimPlot(UM, cells.highlight = refcells, cols = "grey80", cols.highlight = "#AD002A", alpha = 0.3, raster = FALSE) 
    
    #Plot full UMAP, highlight all cells from the same tissue, colour based on samples -> compare similar UMAP patterns across samples
    tissue = sub("^..", "", i) #Identify tissue like this so I don't need to create a separate metadata column for sample + tissue
    tiscells = split(colnames(UM)[UM$tissue == tissue], UM$sample[UM$tissue == tissue]) #Get tissue specific cells but still all samples
    relative = DimPlot(UM, cells.highlight = tiscells, cols = "grey80", cols.highlight = c("#AD002A", "#6C6DA4", "#0A7C2E"), alpha = 0.3, raster = FALSE) 
    
    combed = (empty + grouped) / (reference + relative)
    
    ggsave(file.path(plottingDir, paste0(prefix, i, ".sample_UMAP.pdf")), combed, width = 24, height = 16)
  }
  
  ###Also creating UMAPs specifically to check integration quality -- for i in 1 just because i do not know if i use these object names already
  for (i in 1) {
    
    #Make plots per element to assess whether integration worked
    sample = DimPlot(UM, group.by = "sample", alpha = 0.3, raster = FALSE) 
    tissue = DimPlot(UM, group.by = "tissue", alpha = 0.3, raster = FALSE)
    protocol = DimPlot(UM, group.by = "protocol", alpha = 0.3, raster = FALSE)
    individual = DimPlot(UM, group.by = "individual", alpha = 0.3, raster = FALSE)
    
    #With custom labels
    celltype = DimPlot(UM, group.by = "celltype", alpha = 0.3, raster = FALSE, label = FALSE) + NoLegend()
    celltype = LabelClusters(plot = celltype, id = "celltype", size = 6, fontface = "bold")
    
    cycle = DimPlot(UM, group.by = "Phase", alpha = 0.3, raster = FALSE) #Added this to assess whether cell cycle stage was properly regressed out also, not strictly integration assessment
    
    ggsave(file.path(plottingDir, paste0(prefix, "integration.sample_UMAP.pdf")), sample, width = 24, height = 16)
    ggsave(file.path(plottingDir, paste0(prefix, "integration.tissue_UMAP.pdf")), tissue, width = 24, height = 16)
    ggsave(file.path(plottingDir, paste0(prefix, "integration.protocol_UMAP.pdf")), protocol, width = 24, height = 16)
    ggsave(file.path(plottingDir, paste0(prefix, "integration.individual_UMAP.pdf")), individual, width = 24, height = 16)
    ggsave(file.path(plottingDir, paste0(prefix, "integration.celltype_UMAP.pdf")), celltype, width = 24, height = 16)
    ggsave(file.path(plottingDir, paste0(prefix, "cell_cycle_stage_UMAP.pdf")), cycle, width = 24, height = 16)
    
  }
  ###########
  
  
  
  ########### HOX FEATUREPLOTS
  #Define HOX genes for use in this section, already isolated when making heatmap
  HOXgenes = rownames(a)
  
  #Creating featureplots per individual HOX gene
  message("====== Creating HOX gene featureplots for ", devstage, " ======")
  
  for (i in HOXgenes) {
    
    fp = FeaturePlot(UM, features = i, 
                     split.by = "tissue",
                     cols = c("lightgrey", "#AD002A"),
                     keep.scale = "feature",
                     pt.size = 1.5,
                     order = TRUE)
    formatted = fp + plot_layout(ncol = 3, nrow = 2, guides = "collect") & theme(legend.position = "right")
    
    ggsave(file.path(hoxDir, paste0(prefix, i, ".expression.pdf")), formatted, width = 20, height = 10)
  }
  ###########
  
  
  
  ########### DOTPLOTS
  message("====== Creating dot plots comparing similar tissues for ", devstage, " ======")
  
  #Define tissue names for devstages
  if(devstage == "adult") {
    tissuesTB = c("Bicep", "Thigh")
    tissuesUL = c("Neck", "LowBack")
  } else if(devstage == "embryo") {
    tissuesTB = c("EmbFore", "EmbHind")
    tissuesUL = c("EmbUpper", "EmbLower")
  }
  
  #Subset, filter out cases where a tissue shows up very few times in a cluster to avoid noise in the dotplots
  removenoise = UM@meta.data %>%
    dplyr::count(tissue, celltype) %>%
    dplyr::filter(n >= 10)
  removenoise$group = paste(removenoise$tissue, removenoise$celltype, sep = "_")
  UM$dot_group = paste(UM$tissue, UM$celltype, sep = "_")
  
  UMnoise = subset(UM, subset = dot_group %in% removenoise$group)
  
  BicepThigh = subset(UMnoise, subset = tissue %in% tissuesTB)
  UpperLower = subset(UMnoise, subset = tissue %in% tissuesUL)
  
  #Group tissue and cluster for future identification
  BicepThigh$cluster_tissue = paste(as.character(Idents(BicepThigh)), BicepThigh$tissue, sep = "__")
  UpperLower$cluster_tissue = paste(as.character(Idents(UpperLower)), UpperLower$tissue, sep = "__")
  
  #Define levels
  all.clusters = levels(Idents(UMnoise))
  levelsTB = as.vector(outer(all.clusters, tissuesTB, paste, sep = "__"))
  levelsUL = as.vector(outer(all.clusters, tissuesUL, paste, sep = "__"))
  BicepThigh$cluster_tissue = factor(BicepThigh$cluster_tissue, levels = levelsTB)
  UpperLower$cluster_tissue = factor(UpperLower$cluster_tissue, levels = levelsUL)
  
  
  #Get base DotPlot (does not include clusters without expression)
  pTB = DotPlot(BicepThigh, features = HOXgenes, group.by = "cluster_tissue") + 
    scale_x_discrete(limits = rev) + 
    scale_y_discrete(drop = FALSE) +
    RotatedAxis()
  
  dpTB = pTB@data
  
  dpTB$id = factor(dpTB$id, levels = levelsTB)
  dpTB$features.plot = factor(dpTB$features.plot, levels = HOXgenes)
  
  #Fill DotPlot data with empty data for the clusters not present in this tissue
  dpTB2 = complete(dpTB,
                   id = levelsTB,
                   features.plot = HOXgenes,
                   fill = list(pct.exp = 0, avg.exp = NA_real_, avg.exp.scaled = NA_real_))
  
  dpTB2 = dpTB2 %>% #Status flag to visually filter out 0% expression
    mutate(
      status = case_when(
        pct.exp > 0 ~ "Expressed",
        is.na(avg.exp.scaled) ~ "Cluster/gene not present in tissue",
        TRUE ~ "Not expressed"
      )
    )
  
  #Tag clusters as celltypes instead of just sample + cluster number
  dpTB2 = dpTB2 %>%
    mutate(
      tissue = sub("^\\d+__([A-Za-z]{5}).*$", "\\1", id),
      cluster = sub("^(\\d+).*", "\\1", id),
      celltype = annLookup[as.character(cluster)],
      tag = paste(celltype, tissue, sep = "_")
    )
  
  #Repeat above steps for skin datasets
  pUL = DotPlot(UpperLower, features = HOXgenes, group.by = "cluster_tissue") + 
    scale_x_discrete(limits = rev) + 
    scale_y_discrete(drop = FALSE) +
    RotatedAxis()
  
  dpUL = pUL@data
  
  dpUL$id = factor(dpUL$id, levels = levelsUL)
  dpUL$features.plot = factor(dpUL$features.plot, levels = HOXgenes)
  
  dpUL2 = complete(dpUL,
                   id = levelsUL,
                   features.plot = HOXgenes,
                   fill = list(pct.exp = 0, avg.exp = NA_real_, avg.exp.scaled = NA_real_))
  
  dpUL2 = dpUL2 %>% #Status flag to visually filter out 0% expression
    mutate(
      status = case_when(
        pct.exp > 0 ~ "Expressed",
        is.na(avg.exp.scaled) ~ "Cluster/gene not present in tissue",
        TRUE ~ "Not expressed"
      )
    )
  
  #Tag clusters as celltypes instead of just sample + cluster number
  dpUL2 = dpUL2 %>%
    mutate(
      tissue = sub("^\\d+__([A-Za-z]{5}).*$", "\\1", id),
      cluster = sub("^(\\d+).*", "\\1", id),
      celltype = annLookup[as.character(cluster)],
      tag = paste(celltype, tissue, sep = "_")
    )
  
  #Create function to order HOX genes by paralog groups
  orderHOX = function(genes, letter_order = c("A", "B", "C", "D")) {
    
    g = unique(as.character(genes))
    letter = sub("^HOX([ABCD]).*$", "\\1", g)
    num = suppressWarnings(as.integer(sub("^HOX[ABCD](\\d+).*$", "\\1", g)))
    letter = factor(letter, levels = letter_order)
    g[order(num, letter, g)]
    
  }
  
  #Apply function (these will be the same (at least currently), not sure if this could change)
  OrderedHOXTB = orderHOX(dpTB2$features.plot)
  OrderedHOXUL = orderHOX(dpUL2$features.plot)
  
  #Apply ordering to HOX gene column
  dpTB2$features.plot = factor(dpTB2$features.plot, levels = OrderedHOXTB)
  dpUL2$features.plot = factor(dpUL2$features.plot, levels = OrderedHOXUL)
  
  
  #Create function to order by cluster number in the tag (this is entirely unnecessary and will not be done with the final data, it just looks nicer when the grouping of clusters currently has no meaning)
  orderTags = function(x, tissues) {
    
    labs = unique(as.character(x))
    num = as.integer(sub("^\\s*(\\d+).*","\\1", labs))
    tissue = sub(".*_", "", labs)
    tissueOrder = match(tissue, tissues)
    labs[order(num, tissueOrder, labs)]
    
  }
  
  #Apply meaningless nicer order to x axis labels
  TBtags = orderTags(dpTB2$tag, tissuesTB)
  dpTB2$tag = factor(dpTB2$tag, levels = TBtags)
  
  ULtags = orderTags(dpUL2$tag, tissuesUL)
  dpUL2$tag = factor(dpUL2$tag, levels = ULtags)
  
  
  #Create DotPlots
  duoVerticals = seq(2.5, length(unique(dpUL2$tag)), by = 2)
  #Bicep vs Thigh
  tb = ggplot(dpTB2,
              aes(
                x = tag,
                y = features.plot)
  ) +
    geom_point(
      aes(
        size = pct.exp,
        color = avg.exp.scaled,
        shape = status),
      stroke = 1
    ) +
    scale_size_continuous(breaks = c(0, 25, 50, 75, 100), #Set 0-100 size scale for Michael
                          limits = c(0,100),
                          range = c(1.5, 9),
                          name = "% of cells expressing") +
    scale_color_gradient(low = "lightgrey",
                         high = "red",
                         name = "Z-scaled expression (log2FC)") +
    scale_shape_manual(values =
                         c('Expressed' = 16,
                           'Not expressed' = 4,
                           'Cluster/gene not present in tissue' = NA),
                       breaks =
                         c('Not expressed'),
                       name = "") + #Michael does not want to include expression shape
    scale_y_discrete(
      limits = rev(OrderedHOXTB),
      drop = FALSE
    ) +
    theme_classic() +
    theme(
      panel.background = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      legend.key = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = NA, color = NA)
    ) +
    labs(
      x = "Cluster, cell type, tissue",
      y = "Gene"
    ) +
    guides(
      size = guide_legend(order = 1),
      color = guide_colourbar(order = 2),
      shape = guide_legend(order = 3)
    ) +
    geom_vline(
      xintercept = duoVerticals,
      linetype = "solid"
    ) +
    RotatedAxis()
  
  #Upper vs Lower
  ul = ggplot(dpUL2,
              aes(
                x = tag, 
                y = features.plot)
  ) +
    geom_point(
      aes(
        size = pct.exp, 
        color = avg.exp.scaled, 
        shape = status),
      stroke = 1
    ) +
    scale_size_continuous(breaks = c(0, 25, 50, 75, 100), #Set 0-100 size scale for Michael
                          limits = c(0,100),
                          range = c(1.5, 9),
                          name = "% of cells expressing") +
    scale_color_gradient(low = "lightgrey", 
                         high = "red",
                         name = "Z-scaled expression (log2FC)") +
    scale_shape_manual(values =
                         c('Expressed' = 16,
                           'Not expressed' = 4,
                           'Cluster/gene not present in tissue' = NA),
                       breaks =
                         c('Not expressed'),
                       name = "") + #Michael does not want to include expression shape
    scale_y_discrete(
      limits = rev(OrderedHOXUL), 
      drop = FALSE
    ) +
    theme_classic() +
    theme(
      panel.background = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      legend.key = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = NA, color = NA)
    ) +
    labs(
      x = "Cluster, cell type, tissue",
      y = "Gene"
    ) +
    guides(
      size = guide_legend(order = 1),
      color = guide_colourbar(order = 2),
      shape = guide_legend(order = 3)
    ) +
    geom_vline(
      xintercept = duoVerticals,
      linetype = "solid"
    ) +
    RotatedAxis()
  
  #Saving dotplots
  ggsave(file.path(plottingDir, paste0(prefix, "HOX_dotplot_ForeHindMuscle.pdf")), tb, width = 14, height = 10)
  ggsave(file.path(plottingDir, paste0(prefix, "HOX_dotplot_UpperLowerSkin.pdf")), ul, width = 20, height = 14)
  ###########
  
  
  
  ########### BAR PLOTS
  message("====== Creating barplots capturing variability in HOX paralog expression across tissues for ", devstage, " ======")
  barEXP = AverageExpression(UM, group.by = "tissue")
  barDF = barEXP$RNA
  barHOX = as.data.frame(barDF[grep("^HOX", rownames(barDF)), ])
  
  finalbar = barHOX %>%
    mutate(paralog = sub("^HOX[ABCD]([0-9]+)$", "\\1", rownames(barHOX)),
           paralog = paste0("HOX", paralog)) %>%
    group_by(paralog) %>%
    summarize(across(where(is.numeric), sum), .groups = "drop") %>%
    pivot_longer(cols = -paralog, names_to = "tissue", values_to = "expression") %>%
    mutate(hoxN = as.integer(str_extract(paralog, "\\d+$"))) %>%
    arrange(hoxN) %>%
    mutate(paralog = factor(paralog, levels = unique(paralog))) %>%
    dplyr::select(-hoxN)
  if(devstage == "adult") {
  finalbar$tissue = factor(finalbar$tissue, levels = c("Neck", "Bicep", "LowBack", "Thigh", "Uropygium"))
  } else if(devstage == "embryo") {
    finalbar$tissue = factor(finalbar$tissue, levels = c("EmbUpper", "EmbFore", "EmbLower", "EmbHind"))
  }
  
  relaTiss = finalbar %>%
    group_by(tissue) %>%
    mutate(relative = expression / sum(expression)) %>%
    ungroup()
  
  relaPara = finalbar %>%
    group_by(paralog) %>%
    mutate(relative = expression / sum(expression)) %>%
    ungroup()
  
  #One for all HOX and coloured by tissue
  barplot = ggplot(finalbar) +
    geom_col(aes(x = paralog, y = expression, fill = tissue)) +
    scale_fill_manual(values = tissueCols) +
    theme_classic() +
    labs(
      x = "HOX paralog",
      y = "Summed paralog expression")
  
  #Same as above but relative as a fraction of 0-1
  relaplot = ggplot(relaPara) +
    geom_col(aes(x = paralog, y = relative, fill = tissue)) +
    scale_fill_manual(values = tissueCols) +
    theme_classic() +
    labs(
      x = "HOX paralog",
      y = "Relative tissue contribution")
  
  #Same as 2 above but dodged
  dodgeplot = ggplot(finalbar) +
    geom_col(aes(x = paralog, y = expression, fill = tissue), position = "dodge") +
    scale_fill_manual(values = tissueCols) +
    theme_classic() +
    labs(
      x = "HOX paralog",
      y = "Paralog expression")
  
  #One for all tissues and coloured by HOX
  barplot2 = ggplot(finalbar) +
    geom_col(aes(x = tissue, y = expression, fill = paralog)) +
    scale_fill_manual(values = paralogCols) +
    theme_classic() +
    labs(
      x = "Tissue",
      y = "Total expression")
  
  #Same as above but relative as a fraction of 0-1
  relaplot2 = ggplot(relaTiss) +
    geom_col(aes(x = relative, y = tissue, fill = paralog)) +
    scale_fill_manual(values = paralogCols) +
    theme_classic() +
    scale_y_discrete(limits = rev) +
    labs(
      x = "Relative paralog expression",
      y = "Tissue")
  
  #Same as 2 above but dodged
  dodgeplot2 = ggplot(relaTiss) +
    geom_col(aes(x = tissue, y = relative, fill = paralog), position = "dodge") +
    scale_fill_manual(values = paralogCols) +
    theme_classic() +
    labs(
      x = "Tissue",
      y = "Paralog expression")
  
  #Saving barplots
  ggsave(file.path(plottingDir, paste0(prefix, "profilesHOXparalogs_GENE.pdf")), barplot, width = 8, height = 6)
  ggsave(file.path(plottingDir, paste0(prefix, "profilesHOXparalogs_TISSUE.pdf")), barplot2, width = 8, height = 6)
  ggsave(file.path(plottingDir, paste0(prefix, "relativeHOXparalogs_GENE.pdf")), relaplot, width = 8, height = 6)
  ggsave(file.path(plottingDir, paste0(prefix, "relativeHOXparalogs_TISSUE.pdf")), relaplot2, width = 7, height = 4)
  ggsave(file.path(plottingDir, paste0(prefix, "dodgedHOXparalogs_GENE.pdf")), dodgeplot, width = 8, height = 6)
  ggsave(file.path(plottingDir, paste0(prefix, "dodgedHOXparalogs_TISSUE.pdf")), dodgeplot2, width = 8, height = 6)
  ###########
  
  
  
  ########### HOX PROFILES
  message("====== Creating HOX profiles for ", devstage, " ======")
  profiles = barHOX %>%
    rownames_to_column(var = "gene") %>%
    pivot_longer(cols = -gene, names_to = "tissue", values_to = "expression")
  
  OrderedProfiles = orderHOX(profiles$gene)
  profiles$gene = factor(profiles$gene, levels = OrderedProfiles)
  if(devstage == "adult") {
  profiles$tissue = factor(profiles$tissue, levels = c("Bicep", "Neck", "Uropygium", "Thigh", "LowBack"))
  } else if(devstage == "embryo") {
    profiles$tissue = factor(profiles$tissue, levels = c("EmbUpper", "EmbFore", "EmbLower", "EmbHind"))
  }
  
  prof = ggplot(data = profiles) +
    geom_col(aes(x = gene, y = expression), fill = "#925E9F") +
    facet_wrap(~tissue) + 
    theme_minimal() +
    RotatedAxis()
  
  ggsave(file.path(plottingDir, paste0(prefix, "HOX.profiles.pdf")), prof, width = 18, height = 10)
  ###########
  
  
  
  ########### PER TISSUE HEATMAPS
  message("====== Creating per tissue heatmaps for ", devstage, " ======")
  #Scaled by row (this becomes by gene)
  t1 = apply(barHOX, 1, safeHeat)
  hoxN = as.integer(str_extract(colnames(t1), "\\d+$"))
  hoxL = str_extract(colnames(t1), "(?<=^HOX)[ABCD]")
  hoxO = order(hoxN, hoxL)
  if(devstage == "adult") {
  rowOrder = c("Neck", "Bicep", "LowBack", "Thigh", "Uropygium") 
  } else if(devstage == "embryo") {
    rowOrder = c("EmbUpper", "EmbFore", "EmbLower", "EmbHind") 
  }
  t2 = as.data.frame(t1)[, hoxO]
  t2 = t2[rowOrder, ]
  
  geneHeat = pheatmap(t2, color = color, 
                      cluster_rows = FALSE, 
                      cluster_cols = FALSE,
                      silent = TRUE) 
  
  
  #Scaled by column (this becomes by tissue))
  t3 = apply(barHOX, 2, safeHeat)
  t4 = t(t3)
  t4 = as.data.frame(t4)[, hoxO]
  t4 = t4[rowOrder, ]
  
  tisHeat = pheatmap(t4, color = color, 
                     cluster_rows = FALSE, 
                     cluster_cols = FALSE,
                     silent = TRUE)
  
  #Saving heatmaps
  ggsave(file.path(plottingDir, paste0(prefix, "byTissue_Heatmap.pdf")), tisHeat, width = 8, height = 4)
  ggsave(file.path(plottingDir, paste0(prefix, "byGene_Heatmap.pdf")), geneHeat, width = 8, height = 4)
  ###########



  ########### POSTERIORITY SCORING
  message("====== Creating posteriority scoring for ", devstage, " ======")
  
  ##Set comparisons
  if(devstage == "adult") {
    comparison = list(
      Limb = c("Bicep", "Thigh"),
      Skin = c("Neck", "LowBack"))
  } else if(devstage == "embryo") {
    comparison = list(
      Limb = c("EmbFore", "EmbHind"),
      Skin = c("EmbUpper", "EmbLower"))
  }
  
  ##Calculate posteriority score for each tissue (WHOLE TISSUE LEVEL COMPARISON)
  
  #This weighted stuff is only relevant for HOX 
  #Weighted HOX gene posteriority scores to assess whether HOX genes towards 5' are more present in a tissue
  pbulk = AggregateExpression(object = UM, assays = "RNA", group.by = "sample", return.seurat = FALSE)$RNA
  saveRDS(pbulk, file = file.path(outdir, "SeuratOBJs", paste0(prefix, "sample_aggregated.rds")))
  
  hoxpbulk = pbulk[HOXgenes, , drop = FALSE]
  paralog = as.numeric(sub("^HOX[A-D]", "", rownames(hoxpbulk)))
  tispost = colSums(hoxpbulk * paralog) / colSums(hoxpbulk)
  
  #Store relevant info
  sampmeta = UM[[]] %>% dplyr::select(sample, tissue, individual) %>% distinct()
  postdf = data.frame( #Save for future  interdevstage comparisons
    sample = colnames(hoxpbulk),
    devstage = devstage,
    posteriority = tispost,
    totalHOX = colSums(hoxpbulk),
    detectedHOX = colSums(hoxpbulk > 0)
  )
  postdf = left_join(postdf, sampmeta, by = "sample")
  saveRDS(postdf, file = file.path(outdir, "SeuratOBJs", paste0(prefix, "posteriority.rds")))
  
  #Define comparisons
  samplepairs = list()
  for(c in names(comparison)) {
    
    anterior = sampmeta %>% filter(tissue == comparison[[c]][1])
    posterior = sampmeta %>% filter(tissue == comparison[[c]][2])
    
    if(devstage == "adult") {
     pairs = inner_join(anterior, posterior, by = "individual", suffix = c(".ant", ".post"))
    } else if(devstage == "embryo") {
      pairs = cross_join(anterior, posterior, suffix = c(".ant", ".post"))
    }
   
    pairs$comparison = c
    samplepairs[[c]] = pairs
     
  }
  samplepairs = bind_rows(samplepairs)
  
  #Calculate posteriority deltas
  score = setNames(postdf$posteriority, postdf$sample)
  deltapost = samplepairs
  deltapost$p.ant = score[deltapost$sample.ant]
  deltapost$p.post = score[deltapost$sample.post]
  deltapost$delta = deltapost$p.post - deltapost$p.ant
  deltapost$devstage = devstage
  saveRDS(deltapost, file = file.path(outdir, "SeuratOBJs", paste0(prefix, "deltapost.rds"))) #Save for future  interdevstage comparisons
  
  #Plotting posteriority differences
  monodeltaplot = ggplot(deltapost, aes(x = comparison, y = delta)) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey30"
    ) +
    geom_point(
      aes(color = comparison),
      position = position_jitter(width = 0.1, height = 0)
    ) +
    labs(
      title = paste0("Weighted posteriority comparison for ", devstage),
      y = "Posteriority score delta",
      x = "Tissue type"
    ) +
    theme_bw()
  
  ggsave(file.path(plottingDir, paste0(prefix, "within_devstage_postdelta.pdf")), monodeltaplot, width = 6, height = 6)
  
  
  
  
  ##Log2fc comparison
  pbulkcpm = sweep(pbulk, MARGIN = 2, STATS = colSums(pbulk), FUN = "/") * 1e6 #Normalize to cpm
  saveRDS(pbulkcpm, file = file.path(outdir, "SeuratOBJs", paste0(prefix, "sample_aggregated_cpm.rds"))) #Future comparisons
  geneFC = list()
  
  for(i in seq_len(nrow(samplepairs))) {
    
    ant = samplepairs$sample.ant[i]
    post = samplepairs$sample.post[i]
    
    geneFC[[i]] = data.frame(
      gene = rownames(pbulkcpm),
      sample.ant = ant,
      sample.post = post,
      comparison = samplepairs$comparison[i],
      CPM.ant = pbulkcpm[, ant],
      CPM.post = pbulkcpm[, post],
      log2FC = log2((pbulkcpm[, post] + 1)/(pbulkcpm[, ant] + 1)),
      row.names = NULL
    )
    
  }
  geneFC = bind_rows(geneFC)
  geneFC$devstage = devstage
  saveRDS(geneFC, file = file.path(outdir, "SeuratOBJs", paste0(prefix, "geneFC.rds"))) #Save for future  interdevstage comparisons
  
  #Specify genes of interest
  spgFC = geneFC %>% filter(gene %in% spgs)
  hoxFC = geneFC %>% filter(gene %in% HOXgenes)
  spgnohoxFC = spgFC %>% filter(!gene %in% HOXgenes)
  
  
  #Log2FC heatmap for HOX genes
  hoxFC$gene = factor(hoxFC$gene, levels = rev(HOXgenes))
  hoxFC$contrast = paste(hoxFC$sample.ant, hoxFC$sample.post, sep = " → ")
  hoxFC = hoxFC %>%
    left_join(
      deltapost %>% dplyr::select(sample.ant, sample.post), by = c("sample.ant", "sample.post"))
  
  l2fcheat = ggplot(hoxFC, aes(x = contrast, y = gene, fill = log2FC)) +
    geom_tile(
      colour = "white", linewidth = 0.2) +
    facet_grid(
      . ~ comparison,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_fill_gradient2(
      low = 
    )
  
  ggsave(file.path(plottingDir, paste0(prefix, "within_devstage_postdelta.pdf")), monodeltaplot, width = 6, height = 6)
  




  ###########


  
#}



#!#!#!#!#!#!#!#!#!#!#
#TEMPORARY CODE TO RUN LOCALLY
prefix = paste0("embryo", "_")
seurat = file.path(outdir, "SeuratOBJs", paste0(prefix, "UM.rds"))
markers = file.path(outdir, "SeuratOBJs", paste0(prefix, "markers.rds"))
devstage = "embryo"
#!#!#!#!#!#!#!#!#!#!#

###Run function for both the embryo and adult data
# for(i in c("embryo", "adult")) {
#   
#   prefix = paste0(i, "_")
#   seurat = file.path(outdir, "SeuratOBJs", paste0(prefix, "UM.rds"))
#   markers = file.path(outdir, "SeuratOBJs", paste0(prefix, "markers.rds"))
#   
#   if(file.exists(seurat) && file.exists(markers)) {
#     plotfunc(i)
#   }
#   
# }




# #Below is only code from the previous HOX only pipeline, keeping in for now for future reference
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #Also making a dotplot that includes all five tissues, to avoid having 150 colums only including cell types present in all tissues
# # Five = UM
# # Five$cluster_tissue = paste(as.character(Idents(Five)), Five$tissue, sep = "__")
# # dp5 = DotPlot(Five, features = HOXgenes, group.by = "cluster_tissue")
# # data5 = dp5@data %>%
# #   mutate(
# #     tissue = sub("^\\d+__([A-Za-z]{2}).*$", "\\1", id),
# #     cluster = sub("^(\\d+).*", "\\1", id),
# #     celltype = annLookup[as.character(cluster)],
# #     tag = paste(celltype, tissue, sep = "_"),
# #     status = case_when(
# #       pct.exp > 0 ~ "Expressed",
# #       is.na(avg.exp.scaled) ~ "Cluster/gene not present in tissue",
# #       TRUE ~ "Not expressed")
# #   ) %>%
# #   group_by(cluster) %>%
# #   filter(n_distinct(tissue) == n_distinct(Five$tissue)) %>%
# #   ungroup() 
# # 
# # 
# # #tagsALL = orderTags(data5$tag)
# # #data5$tag = factor(data5$tag, levels = tagsALL)
# # 
# # alltagorder = c("Bi","Th","Ne","Lo","Ur")
# # ALLorder = data5 %>%
# #   distinct(tag) %>%
# #   mutate(
# #     ct = str_extract(tag, "^\\d+") |> as.numeric(),
# #     ts = sub(".*_","", tag)
# #   ) %>%
# #   arrange(ct, match(ts, alltagorder)) %>%
# #   pull(tag)
# # 
# # data5$tag = factor(data5$tag, levels = ALLorder)
# # 
# # 
# # Allverticals = seq(5.5, length(unique(data5$tag)), by = 5)
# # 
# # dotfive = ggplot(data5,
# #        aes(
# #          x = tag, 
# #          y = features.plot)
# # ) +
# #   geom_point(
# #     aes(
# #       size = pct.exp, 
# #       color = avg.exp.scaled, 
# #       shape = status),
# #     stroke = 1
# #   ) +
# #   scale_size_continuous(breaks = c(0, 25, 50, 75, 100), #Set 0-100 size scale for Michael
# #                         limits = c(0,100),
# #                         range = c(1.5, 9),
# #                         name = "% of cells expressing") +
# #   scale_color_gradient(low = "lightgrey", 
# #                        high = "red",
# #                        name = "Z-scaled expression (log2FC)") +
# #   scale_shape_manual(values =
# #                        c('Expressed' = 16,
# #                          'Not expressed' = 4,
# #                          'Cluster/gene not present in tissue' = NA),
# #                      breaks =
# #                        c('Not expressed'),
# #                      name = "") + #Michael does not want to include expression shape
# #   scale_y_discrete(
# #     limits = rev, 
# #     drop = FALSE
# #   ) +
# #   theme_classic() +
# #   theme(
# #     panel.background = element_rect(
# #       color = "black",
# #       fill = NA,
# #       linewidth = 1
# #     ),
# #     legend.key = element_rect(fill = NA, color = NA),
# #     legend.background = element_rect(fill = NA, color = NA)
# #   ) +
# #   labs(
# #     x = "Cluster, cell type, tissue",
# #     y = "Gene"
# #   ) +
# #   guides(
# #     size = guide_legend(order = 1),
# #     color = guide_colourbar(order = 2),
# #     shape = guide_legend(order = 3)
# #   ) +
# #   RotatedAxis() +
# #   geom_vline(
# #     xintercept = Allverticals, 
# #     linetype = "solid")
# # 
# # ggsave(file.path(plottingDir, "HOX_dotplot_AllTissues.pdf"), dotfive, width = 14, height = 10)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###Also making heatmaps per tissue since this is included in master file (page 25 also)
# 
# 
# 
# #NOT VERY NICE, KEEP FOR ADDITIONAL SAMPLES, INCLUDE OR EXCLUDE ACCORDINGLY
# #Making barplots capturing expression ratio between anterior and posterior, both ways for now, could choose either or set lims
# #rate = AverageExpression(UM, group.by = "tissue", layer = "data")$RNA
# #ratioFront = as.data.frame(rate[grep("^HOX", rownames(rate)), ]) %>%
# #  mutate(
# #    muscleRatioFront = log2(Bicep / Thigh),
# #    skinRatioFront = log2(Neck / LowBack)
# #  ) %>%
# #  dplyr::select(-c(1:5)) %>%
# #  rownames_to_column(var = "gene") %>%
# #  pivot_longer(
# #    cols = -gene,
# #    names_to = "group",
# #    values_to = "ratio"
# #  ) %>%
# #  mutate(type = sub("RatioFront", "", group))
# 
# #ratioBack = as.data.frame(rate[grep("^HOX", rownames(rate)), ]) %>%
# #  mutate(
# #    muscleRatioTail = log2(Thigh / Bicep),
# #    skinRatioTail = log2(LowBack / Neck)
# #  ) %>%
# #  dplyr::select(-c(1:5)) %>%
# #  rownames_to_column(var = "gene") %>%
# #  pivot_longer(
# #    cols = -gene,
# #    names_to = "group",
# #    values_to = "ratio"
# #  ) %>%
# #  mutate(type = sub("RatioTail", "", group))
# #
# #ratioFlvls = orderHOX(ratioFront$gene)
# #ratioFront$gene = factor(ratioFront$gene, levels = ratioFlvls)
# #
# #ratioBlvls = orderHOX(ratioBack$gene)
# #ratioBack$gene = factor(ratioBack$gene, levels = ratioBlvls)
# #
# #
# #ggplot(ratioFront) +
# #  geom_col(aes(x = gene, y = ratio, fill = type), position = "dodge") +
# #  labs(
# #    title = "Ratio front:back"
# #  )
# #ggplot(ratioBack) +
# #  geom_col(aes(x = gene, y = ratio, fill = type), position = "dodge") +
# #  labs(
# #    title = "Ratio back:front")
# 
# 
# #################################################
# #Pseudobulk and volcano plots
# #################################################
# 
# ###Perform pseudobulks
# pseudoDir = file.path(plottingDir, "Pseudobulk_integration_check")
# if (!dir.exists(pseudoDir)) {
#   dir.create(pseudoDir, recursive = TRUE)
# }
# 
# UMct = UM #Created only to isolate pseudobulk object from general object to avoid changing contents for other code
# 
# ##Perform pseudobulk - INTEGRATION CHECK - CELLTYPES SHOULD BE SIMILAR BETWEEN SAMPLES CONTAINING THOSE CELLS
# plan(sequential)
# pseudo = AggregateExpression(UMct, assays = "RNA", return.seurat = TRUE, group.by = c("sample", "celltype"))
# plan(multicore, workers = 12)
# 
# #Add informative metadata
# Idents(pseudo) = pseudo$celltype
# pseudo$sample = sub("_(.*)$", "", Cells(pseudo))
# pseudo$celltype = sub("^[^_]+_", "", Cells(pseudo))
# meta_map = unique(UMct@meta.data[, c("sample", "tissue", "individual", "protocol")])
# rownames(meta_map) = meta_map$sample
# pseudo$tissue = meta_map[pseudo$sample, "tissue"]
# pseudo$individual = meta_map[pseudo$sample, "individual"]
# pseudo$protocol = meta_map[pseudo$sample, "protocol"]
# 
# tail(Cells(pseudo)) #check
# 
# 
# for(i in unique(pseudo$celltype)) {
#   
# #Determine which samples contain very few cells for a given cell type to eliminate noise
#   keep_samples = UMct@meta.data %>%
#     dplyr::count(sample, celltype) %>%
#     filter(celltype == i, n >= 100) %>% #Minimum # of cells of a celltype to reduce noise
#     pull(sample)
# 
#   if(length(keep_samples) < 2) {
#     message("Skipping", i, ": fewer than two samples contain more than 100 of this cell type")
#     next
#   }
#   
# #Get celltype subset with only samples containing at least 25 cells of that cell type 
#   sub = subset(pseudo, subset = celltype == i & sample %in% keep_samples)
#     
#   mat = GetAssayData(sub, layer = "data")
#   cor_mat = cor(mat)
# 
#   h = pheatmap(cor_mat,
#            clustering_distance_rows = "correlation",
#            clustering_distance_cols = "correlation",
#            color = colorRampPalette(c("#00468B", "white", "#AD000A"))(100),
#            breaks = seq(0, 1, length.out = 101),
#            angle_col = 45)
#   
#   ggsave(file.path(pseudoDir, paste0("pseudobulk_correlation_heatmap_", make.names(i), ".pdf")), h, width = 8, height = 6)
#   
# }
# 
# 
# ####################################
# 
# 
# ##Perform pseudobulk - DE CHECK BETWEEN TISSUES 
# tisDir = file.path(plottingDir, "Pseudobulk_DESeq_by_tissue")
# if (!dir.exists(tisDir)) {
#   dir.create(tisDir, recursive = TRUE)
# }
# 
# plan(sequential)
# tispseudo = AggregateExpression(UMct, assays = "RNA", return.seurat = TRUE, group.by = c("sample", "tissue"))
# plan(multicore, workers = 12)
# 
# #Add informative metadata
# tispseudo$sample = sub("_(.*)$", "", Cells(tispseudo))
# tispseudo$celltype = sub("^[^_]+_", "", Cells(tispseudo))
# 
# Idents(tispseudo) = tispseudo$celltype
# 
# meta_map2 = unique(UMct@meta.data[, c("sample", "tissue", "individual", "protocol")])
# rownames(meta_map2) = meta_map2$sample
# tispseudo$tissue = meta_map2[tispseudo$sample, "tissue"]
# tispseudo$individual = meta_map2[tispseudo$sample, "individual"]
# tispseudo$protocol = meta_map2[tispseudo$sample, "protocol"]
# 
# tiscounts = GetAssayData(tispseudo, assay = "RNA", layer = "counts")
# 
# tismeta = tispseudo@meta.data[, c("sample", "tissue", "individual", "protocol")]
# rownames(tismeta) = colnames(tiscounts)
# 
# #Running DESeq on it 
# tisdds = DESeqDataSetFromMatrix(
#   countData = round(tiscounts),
#   colData = tismeta,
#   design = ~ tissue
# )
# tisdds = DESeq(tisdds)
# 
# skintisres = results(
#   tisdds,
#   contrast = c("tissue", "Neck", "LowBack")
# )
# 
# muscletisres = results(
#   tisdds,
#   contrast = c("tissue", "Bicep", "Thigh")
# )
# 
# #Defining volcano plot function
# makevolcano = function(results, title) {
#   
#   df = as.data.frame(results) %>%
#     rownames_to_column("gene") %>%
#     mutate(
#       neglog10 = -log10(padj),
#       sig = case_when(
#         padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
#         padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
#         TRUE ~ "Non-significant change"
#       )
#     )
#   ggplot(df, aes(x = log2FoldChange, y = neglog10)) +
#     geom_point(aes(colour = sig, alpha = sig), size = 1.5) +
#     geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#     geom_text_repel(data = df %>% filter(grepl("^HOX", gene), padj < 0.05),
#                     aes(label = gene),
#                     size = 4,
#                     fontface = "bold",
#                     max.overlaps = 30) +
#     scale_color_manual(values = c("Upregulated" = "#AD002A", "Downregulated" = "#3B81AB", "Non-significant change" = "black")) +
#     scale_alpha_manual(values = c("Upregulated" = 0.8, "Downregulated" = 0.8, "Non-significant change" = 0.2), guide = "none") +
#     labs(title = title,
#          x = "log2FC",
#          y = "-log10(padj)",
#          color = "DE status") +
#     theme_minimal()
#   
# }
# 
# #Making volcano plots
# skintis = makevolcano(skintisres, "Upper neck vs lower back skin")
# muscletis = makevolcano(muscletisres, "Bicep vs thigh muscle")
# 
# ggsave(file.path(tisDir, "skin_volcanoplot.pdf"), skintis, width = 8, height = 6)
# ggsave(file.path(tisDir, "muscle_volcanoplot.pdf"), muscletis, width = 8, height = 6)
# 
# #Also adding a combined plot for comparison of the tissue types
# combinedvolc = (skintis + muscletis) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
# ggsave(file.path(tisDir, "combined_volcanoplot.pdf"), combinedvolc, width = 8, height = 6)
# 
# 
# ###############################
# 
# 
# ###Also making these sample DE volcano plots between tissues with a minimum of 10 cells in a cell type filter present to see whether contamination affects these a lot
# 
# #####For skin
# keepskin = UMct@meta.data %>%
#   filter(tissue %in% c("Neck", "LowBack")) %>%
#   dplyr::count(tissue, celltype) %>%
#   pivot_wider(names_from = tissue, values_from = n, values_fill = 0) %>%
#   filter(Neck >= 10 & LowBack >= 10) %>%
#   pull(celltype)
# 
# UMcleanskin = subset(UMct, subset = tissue %in% c("Neck", "LowBack") & celltype %in% keepskin)
# 
# plan(sequential)
# cleanskin = AggregateExpression(UMcleanskin, assays = "RNA", return.seurat = TRUE, group.by = c("sample", "tissue"))
# plan(multicore, workers = 12)
# 
# #Add informative metadata
# cleanskin$sample = sub("_(.*)$", "", Cells(cleanskin))
# 
# meta_map9 = unique(UMct@meta.data[, c("sample", "tissue", "individual", "protocol")])
# rownames(meta_map9) = meta_map9$sample
# cleanskin$tissue = meta_map9[cleanskin$sample, "tissue"]
# cleanskin$individual = meta_map9[cleanskin$sample, "individual"]
# cleanskin$protocol = meta_map9[cleanskin$sample, "protocol"]
# 
# skincounts = GetAssayData(cleanskin, assay = "RNA", layer = "counts")
# 
# skinmeta = cleanskin@meta.data[, c("sample", "tissue", "individual", "protocol")]
# rownames(skinmeta) = colnames(skincounts)
# 
# #Running DESeq on it 
# skindds = DESeqDataSetFromMatrix(
#   countData = round(skincounts),
#   colData = skinmeta,
#   design = ~ tissue
# )
# skindds = DESeq(skindds)
# 
# cleanskinres = results(
#   skindds,
#   contrast = c("tissue", "Neck", "LowBack")
# )
# 
# #Make volcano
# cleanskinvolc = makevolcano(cleanskinres, "Filtered upper neck vs lower back skin")
# ggsave(file.path(tisDir, "filtered_skin_volcanoplot.pdf"), cleanskinvolc, width = 8, height = 6)
# 
# 
# #####For muscle
# keepmuscle = UMct@meta.data %>%
#   filter(tissue %in% c("Bicep", "Thigh")) %>%
#   dplyr::count(tissue, celltype) %>%
#   pivot_wider(names_from = tissue, values_from = n, values_fill = 0) %>%
#   filter(Bicep >= 10 & Thigh >= 10) %>%
#   pull(celltype)
# 
# UMcleanmuscle = subset(UMct, subset = tissue %in% c("Bicep", "Thigh") & celltype %in% keepmuscle)
# 
# plan(sequential)
# cleanmuscle = AggregateExpression(UMcleanmuscle, assays = "RNA", return.seurat = TRUE, group.by = c("sample", "tissue"))
# plan(multicore, workers = 12)
# 
# #Add informative metadata
# cleanmuscle$sample = sub("_(.*)$", "", Cells(cleanmuscle))
# 
# meta_map8 = unique(UMct@meta.data[, c("sample", "tissue", "individual", "protocol")])
# rownames(meta_map8) = meta_map8$sample
# cleanmuscle$tissue = meta_map8[cleanmuscle$sample, "tissue"]
# cleanmuscle$individual = meta_map8[cleanmuscle$sample, "individual"]
# cleanmuscle$protocol = meta_map8[cleanmuscle$sample, "protocol"]
# 
# musclecounts = GetAssayData(cleanmuscle, assay = "RNA", layer = "counts")
# 
# musclemeta = cleanmuscle@meta.data[, c("sample", "tissue", "individual", "protocol")]
# rownames(musclemeta) = colnames(musclecounts)
# 
# #Running DESeq on it 
# muscledds = DESeqDataSetFromMatrix(
#   countData = round(musclecounts),
#   colData = musclemeta,
#   design = ~ tissue
# )
# muscledds = DESeq(muscledds)
# 
# cleanmuscleres = results(
#   muscledds,
#   contrast = c("tissue", "Bicep", "Thigh")
# )
# 
# #Make volcano
# cleanmusclevolc = makevolcano(cleanmuscleres, "Filtered bicep vs thigh")
# ggsave(file.path(tisDir, "filtered_muscle_volcanoplot.pdf"), cleanmusclevolc, width = 8, height = 6)
# 
# #Make combined volcano plot for filtered data
# combinedcleanvolc = (cleanskinvolc + cleanmusclevolc) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
# ggsave(file.path(tisDir, "combined_filtered_volcanoplot.pdf"), combinedcleanvolc, width = 8, height = 6)
# 
# 
# 
# #Also making DE volcanos where tissues are compared but individual cell types are removed as negative controls
# #Should allow us to see if DE between tissues changes significantly as a result of removing any single cell type
# 
# #Defining a function to run 1 by 1 leave out DESeqs
# 
# leaveout = function(obj, tissues, contrast, outdir, prefix) {
#   
#   dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
#   
#   celltypes = sort(unique(obj$celltype))
#   removeSet = c("NONE_all_cells", celltypes)
#   reslist = list()
#   summarytable = data.frame()
#   
#   for(i in removeSet) {
#     
#     message("Running leave one cell type out DESeq for ", i)
#     
#     if(i == "NONE_all_cells") {
#       subobj = subset(obj, subset = tissue %in% tissues)
#     } else {
#       subobj = subset(obj, subset = tissue %in% tissues & celltype != i)
#     }
#     
#     #Pseudobulking
#     plan(sequential)
#     pseudo = AggregateExpression(subobj, assays = "RNA", return.seurat = TRUE, group.by = c("sample", "tissue"))
#     plan(multicore, workers = 12)
#     
#     #Setting info
#     pseudo$sample = sub("_(.*)$", "", Cells(pseudo))
#     metamap = unique(obj@meta.data[, c("sample", "tissue", "individual", "protocol")])
#     rownames(metamap) = metamap$sample
#     pseudo$tissue = metamap[pseudo$sample, "tissue"]
#     pseudo$individual = metamap[pseudo$sample, "individual"]
#     pseudo$protocol = metamap[pseudo$sample, "protocol"]
#     
#     counts = GetAssayData(pseudo, assay = "RNA", layer = "counts")
#     
#     meta = pseudo@meta.data[, c("sample", "tissue", "individual", "protocol")]
#     rownames(meta) = colnames(counts)
#     meta = droplevels(meta)
#     meta$tissue = factor(meta$tissue, levels = tissues)
#     
#     #Skipping cell types only present in one tissue
#     if(length(unique(meta$tissue)) < 2) {
#       message("Skipping ", i, ": only one tissue left")
#       next
#     }
#     
#     #Run DESeq
#     dds = DESeqDataSetFromMatrix(countData = round(counts), colData = meta, design = ~ tissue)
#     
#     #Default dds = DESeq(dds) was giving issues locally, make it use a different fittype if necessary
#     dds = tryCatch({
#       message("Trying fitType parametric")
#       DESeq(dds, fitType = "parametric")
#       }, error = function(e1) {
#         message("Parametric failed: ", e1$message)
#         message("Trying fitType local")
#         tryCatch({
#           DESeq(dds, fitType = "local")
#         }, error = function(e2) {
#           message("Local failed: ", e2$message)
#           message("Trying fitType mean")
#           tryCatch({
#             DESeq(dds, fitType = "mean")
#           }, error = function(e3) {
#             message("Mean also failed: ", e3$message)
#             message("Skipping ", i)
#             return(NULL)
#           })
#         })
#       })
#     if(is.null(dds)) next
#     
#     res = results(dds, contrast = contrast)
#     resdf = as.data.frame(res) %>%
#       rownames_to_column("gene") %>%
#       mutate(removed_celltype = i)
#     
#     #Ensure naming
#     safect = make.names(i)
#     reslist[[safect]] = resdf
#     
#     #Make and save plots
#     p = makevolcano(res, paste0(prefix, " | removed: ", i))
#     ggsave(file.path(outdir, paste0(prefix, "_removed_", safect, "_volcano.pdf")), p, width = 7, height = 5)
#     
#     #Store numerical output
#     hox = grepl("^HOX", resdf$gene)
#     summarytable = rbind(summarytable, data.frame(
#       removed_celltype = i,
#       n_DE_all = sum(resdf$padj < 0.05, na.rm = TRUE),
#       n_DE_absLog2FC1 = sum(resdf$padj < 0.05 & abs(resdf$log2FoldChange) > 1, na.rm = TRUE),
#       n_HOX_DE = sum(hox & resdf$padj < 0.05, na.rm = TRUE),
#       n_HOX_DE_absLog2FC1 = sum(hox & resdf$padj < 0.05 & abs(resdf$log2FoldChange) > 1, na.rm = TRUE)))
#     
#   }
#   
#   write.xlsx(reslist, file.path(outdir, paste0(prefix, "leave_one_celltype_out_all_results.xlsx")), overwrite = TRUE)
#   write.xlsx(summarytable, file.path(outdir, paste0(prefix, "leave_one_celltype_out_summary.xlsx")), overwrite = TRUE)
#   
# }
# 
# #Using the function to get the leave one out volcanos
# leaveout(obj = UMct, 
#          tissues = c("Bicep", "Thigh"), 
#          contrast = c("tissue", "Bicep", "Thigh"), 
#          outdir = file.path(tisDir, "Muscle_leave_one_out"),
#          prefix = "Bicep_vs_Thigh")
# 
# leaveout(obj = UMct,
#          tissues = c("Neck", "LowBack"),
#          contrast = c("tissue", "Neck", "LowBack"),
#          outdir = file.path(tisDir, "Skin_leave_one_out"),
#          prefix = "Neck_vs_LowBack")
# 
# 
# 
# 
# ###############################
# 
# 
# 
# 
# 
# 
# 
# 
# ##Perform pseudobulk - DE CHECK BETWEEN TISSUES PER CELLTYPE
# cttisDir = file.path(plottingDir, "Pseudobulk_DESeq_by_tissue_per_celltype")
# if (!dir.exists(cttisDir)) {
#   dir.create(cttisDir, recursive = TRUE)
# }
# 
# plan(sequential)
# ctpseudo = AggregateExpression(UMct, assays = "RNA", return.seurat = TRUE, group.by = c("sample", "tissue", "celltype"))
# plan(multicore, workers = 12)
# 
# #Add informative metadata
# parts = do.call(rbind, strsplit(Cells(ctpseudo), "_"))
# 
# ctpseudo$sample = parts[,1]
# ctpseudo$tissue = parts[,2]
# ctpseudo$celltype = parts[,3]
# 
# Idents(ctpseudo) = ctpseudo$celltype
# 
# meta_map2 = unique(UMct@meta.data[, c("sample", "individual", "protocol")])
# rownames(meta_map2) = meta_map2$sample
# ctpseudo$individual = meta_map2[ctpseudo$sample, "individual"]
# ctpseudo$protocol = meta_map2[ctpseudo$sample, "protocol"]
# 
# ctcounts = GetAssayData(ctpseudo, assay = "RNA", layer = "counts")
# 
# ctmeta = ctpseudo@meta.data[, c("sample", "tissue", "celltype", "individual", "protocol")]
# rownames(ctmeta) = colnames(ctcounts)
# 
# 
# #Perform DESeq and make plots per cell type
# #For skin
# skinreslist = list()
# skinplotlist = list()
# skinsignif = data.frame()
# 
# for (i in unique(ctmeta$celltype)) {
#   
#   keep = ctmeta$celltype == i &
#     ctmeta$tissue %in% c("Neck", "LowBack")
#   message("RUNNING DESEQ FOR ", i)
#   
#   count = ctcounts[, keep]
#   meta = ctmeta[keep,]
#   
#   print(table(meta$tissue))
#   print(meta[, c("sample", "tissue", "celltype", "individual")])
#   
#   #Not all cell types are present in both tissue, so filter those cell types out
#   if(length(unique(meta$tissue)) < 2) {
#     message("Skipping ", i, " because only one of the tissues contains this cell type")
#     next
#   }
# 
#   #We also cannot compare when a cell type is only present in 1 individual since we can't account for effects of an individual
#   if(length(unique(meta$individual)) < 2) {
#     message("Skipping ", i, " because fewer than 2 individuals contain this cell type")
#     next
#   }
#   
#   #Even when the cell types are present in both tissue, they may not be present in every sample, so filter those out also
#   pairtab = table(meta$individual, meta$tissue)
#   
#   if(!all(c("Neck", "LowBack") %in% colnames(pairtab)) ||
#      any(pairtab[, c("Neck", "LowBack")] == 0)) {
#     message("Skipping ", i, " because not all samples of both tissues contain the celltype")
#     next
#   }
#   
#   #Adding this because I am getting errors
#   meta = droplevels(meta)
#   meta$tissue = factor(meta$tissue, levels = c("Neck", "LowBack"))
#   meta$individual = factor(meta$individual)
# 
#   #Continue 
#   dds = DESeqDataSetFromMatrix(
#     countData = round(count),
#     colData = meta,
#     design = ~ individual + tissue
#   )
#   
#   dds = DESeq(dds)
#   
#   skinreslist[[i]] = results(dds, contrast = c("tissue", "Neck", "LowBack"))
#   skinplotlist[[i]] = makevolcano(skinreslist[[i]], paste0(i, ": upper neck vs lower back"))
#   
#   safename = gsub("/", "_", i) #Cell types with a / in the name mess up ggsave
#   ggsave(file.path(cttisDir, paste0("skin_", safename, "_volcanoplot.pdf")), skinplotlist[[i]], width = 7, height = 5)
#   
#   
#   #Also make a df for a barchart reflecting whether any HOX genes are significantly differentially expressed in any specific cell types
#   resdf = as.data.frame(skinreslist[[i]])
#   hox = grepl("^HOX", rownames(resdf))
#   
#   skinsignif = rbind(
#     skinsignif,
#     data.frame(
#       celltype = i,
#       comparison = "Neck_vs_LowBack",
#       n_hox_sig = sum(
#         hox &
#           resdf$padj < 0.05,
#         na.rm = TRUE
#       ),
#       n_hox_logfc_sig = sum(
#         hox &
#           resdf$padj < 0.05 &
#           abs(resdf$log2FoldChange) > 1,
#         na.rm = TRUE
#       ),
#       n_hox_up = sum(
#         hox &
#           resdf$padj < 0.05 &
#           resdf$log2FoldChange > 1,
#         na.rm = TRUE
#       ),
#       n_hox_down = sum(
#         hox &
#           resdf$padj < 0.05 &
#           resdf$log2FoldChange < -1,
#         na.rm = TRUE
#       )
#     )
#     
#   )
#   
# }
# 
# #Make the plot and save it
# skinhoxcount = ggplot(skinsignif, aes(x = reorder(celltype, n_hox_logfc_sig), y = n_hox_logfc_sig)) +
#   geom_col() +
#   coord_flip() +
#   theme_minimal()
# ggsave(file.path(cttisDir, "skin_SignificantHOX_DE_count_per_celltype.pdf"), skinhoxcount, width = 10, height = 8)
# 
# 
# #For muscle
# musclereslist = list()
# muscleplotlist = list()
# musclesignif = data.frame()
# 
# for (i in unique(ctmeta$celltype)) {
#   
#   keep = ctmeta$celltype == i &
#     ctmeta$tissue %in% c("Bicep", "Thigh")
#   message("RUNNING DESEQ FOR ", i)
#   
#   count = ctcounts[, keep]
#   meta = ctmeta[keep,]
#   
#   print(table(meta$tissue))
#   print(meta[, c("sample", "tissue", "celltype", "individual")])
#   
#   #Not all cell types are present in both tissue, so filter those cell types out
#   if(length(unique(meta$tissue)) < 2) {
#     message("Skipping", i, "because only one of the tissues contains this cell type")
#     next
#   }
# 
#   #We also cannot compare when a cell type is only present in 1 individual since we can't account for effects of an individual
#   if(length(unique(meta$individual)) < 2) {
#     message("Skipping ", i, " because fewer than 2 individuals contain this cell type")
#     next
#   }
#   
#   #Even when the cell types are present in both tissue, they may not be present in every sample, so filter those out also
#   pairtab = table(meta$individual, meta$tissue)
#   
#   if(!all(c("Bicep", "Thigh") %in% colnames(pairtab)) ||
#      any(pairtab[, c("Bicep", "Thigh")] == 0)) {
#     message("Skipping", i, "because not all samples of both tissues contain the celltype")
#     next
#   }
# 
#   #Adding this because I am getting errors 
#   meta = droplevels(meta)
#   meta$tissue = factor(meta$tissue, levels = c("Bicep", "Thigh"))
#   meta$individual = factor(meta$individual)
# 
#   #Continue   
#   dds = DESeqDataSetFromMatrix(
#     countData = round(count),
#     colData = meta,
#     design = ~ individual + tissue
#   )
#   
#   dds = DESeq(dds)
#   
#   musclereslist[[i]] = results(dds, contrast = c("tissue", "Bicep", "Thigh"))
#   muscleplotlist[[i]] = makevolcano(musclereslist[[i]], paste0(i, ": bicep vs thigh"))
#   
#   safename = gsub("/", "_", i) #Cell types with a / in the name mess up ggsave
#   
#   ggsave(file.path(cttisDir, paste0("muscle_", safename, "_volcanoplot.pdf")), muscleplotlist[[i]], width = 7, height = 5)
#   
#   
#   #Also make a df for a barchart reflecting whether any HOX genes are significantly differentially expressed in any specific cell types
#   resdf = as.data.frame(musclereslist[[i]])
#   hox = grepl("^HOX", rownames(resdf))
#   
#   musclesignif = rbind(
#     musclesignif,
#     data.frame(
#       celltype = i,
#       comparison = "Bicep_vs_Thigh",
#       n_hox_sig = sum(
#         hox &
#           resdf$padj < 0.05,
#         na.rm = TRUE
#       ),
#       n_hox_logfc_sig = sum(
#         hox &
#           resdf$padj < 0.05 &
#           abs(resdf$log2FoldChange) > 1,
#         na.rm = TRUE
#       ),
#       n_hox_up = sum(
#         hox &
#           resdf$padj < 0.05 &
#           resdf$log2FoldChange > 1,
#         na.rm = TRUE
#       ),
#       n_hox_down = sum(
#         hox &
#           resdf$padj < 0.05 &
#           resdf$log2FoldChange < -1,
#         na.rm = TRUE
#       )
#     )
#     
#   )
#   
# }
# 
# musclehoxcount = ggplot(musclesignif, aes(x = reorder(celltype, n_hox_logfc_sig), y = n_hox_logfc_sig)) +
#   geom_col() +
#   coord_flip() +
#   theme_minimal()
# 
# ggsave(file.path(cttisDir, "muscle_SignificantHOX_DE_count_per_celltype.pdf"), musclehoxcount, width = 10, height = 8)
# 
# 
# 
# ############################################################################
# 
# ##Also making filtered versions of the per cell type volcano plots requiring at least 10 cells of a cell type per tissue
# 
# #Perform DESeq and make plots per cell type
# #For skin
# filskinreslist = list()
# filskinplotlist = list()
# filskinsignif = data.frame()
# 
# for (i in unique(ctmeta$celltype)) {
#   
#   keep = ctmeta$celltype == i &
#     ctmeta$tissue %in% c("Neck", "LowBack")
#   message("RUNNING DESEQ FOR ", i)
#   
#   count = ctcounts[, keep]
#   meta = ctmeta[keep,]
#   
#   print(table(meta$tissue))
#   print(meta[, c("sample", "tissue", "celltype", "individual")])
#   
#   cellcounts = UMct@meta.data %>%
#     filter(tissue %in% c("Neck", "LowBack"),
#            celltype == i) %>%
#     dplyr::count(sample, tissue)
#   
#   
#   #Not all cell types are present in both tissue, so filter those cell types out
#   if(length(unique(meta$tissue)) < 2) {
#     message("Skipping ", i, " because only one of the tissues contains this cell type")
#     next
#   }
#   
#   #For this version we only want cell types for which there are 10 or more cells of that cell type present in both tissues
#   if(any(cellcounts$n < 10)) {
#     message("Skipping ", i, ": one or more pseudobulks have fewer than 10 cells of this cell type")
#     next
#   }
#   
#   #We also cannot compare when a cell type is only present in 1 individual since we can't account for effects of an individual
#   if(length(unique(meta$individual)) < 2) {
#     message("Skipping ", i, " because fewer than 2 individuals contain this cell type")
#     next
#   }
#   
#   #Even when the cell types are present in both tissue, they may not be present in every sample, so filter those out also
#   pairtab = table(meta$individual, meta$tissue)
#   
#   if(!all(c("Neck", "LowBack") %in% colnames(pairtab)) ||
#      any(pairtab[, c("Neck", "LowBack")] == 0)) {
#     message("Skipping ", i, " because not all samples of both tissues contain the celltype")
#     next
#   }
#   
#   #Adding this because I am getting errors
#   meta = droplevels(meta)
#   meta$tissue = factor(meta$tissue, levels = c("Neck", "LowBack"))
#   meta$individual = factor(meta$individual)
#   
#   #Continue 
#   dds = DESeqDataSetFromMatrix(
#     countData = round(count),
#     colData = meta,
#     design = ~ individual + tissue
#   )
#   
#   dds = DESeq(dds)
#   
#   filskinreslist[[i]] = results(dds, contrast = c("tissue", "Neck", "LowBack"))
#   filskinplotlist[[i]] = makevolcano(filskinreslist[[i]], paste0(i, ": filtered upper neck vs lower back"))
#   
#   safename = gsub("/", "_", i) #Cell types with a / in the name mess up ggsave
#   ggsave(file.path(cttisDir, paste0("filtered_skin_", safename, "_volcanoplot.pdf")), filskinplotlist[[i]], width = 7, height = 5)
#   
#   
#   #Also make a df for a barchart reflecting whether any HOX genes are significantly differentially expressed in any specific cell types
#   resdf = as.data.frame(filskinreslist[[i]])
#   hox = grepl("^HOX", rownames(resdf))
#   
#   filskinsignif = rbind(
#     filskinsignif,
#     data.frame(
#       celltype = i,
#       comparison = "Neck_vs_LowBack",
#       n_hox_sig = sum(
#         hox &
#           resdf$padj < 0.05,
#         na.rm = TRUE
#       ),
#       n_hox_logfc_sig = sum(
#         hox &
#           resdf$padj < 0.05 &
#           abs(resdf$log2FoldChange) > 1,
#         na.rm = TRUE
#       ),
#       n_hox_up = sum(
#         hox &
#           resdf$padj < 0.05 &
#           resdf$log2FoldChange > 1,
#         na.rm = TRUE
#       ),
#       n_hox_down = sum(
#         hox &
#           resdf$padj < 0.05 &
#           resdf$log2FoldChange < -1,
#         na.rm = TRUE
#       )
#     )
#     
#   )
#   
# }
# 
# #Make the plot and save it
# filskinhoxcount = ggplot(filskinsignif, aes(x = reorder(celltype, n_hox_logfc_sig), y = n_hox_logfc_sig)) +
#   geom_col() +
#   coord_flip() +
#   theme_minimal()
# ggsave(file.path(cttisDir, "filtered_skin_SignificantHOX_DE_count_per_celltype.pdf"), filskinhoxcount, width = 10, height = 8)
# 
# 
# #For muscle
# filmusclereslist = list()
# filmuscleplotlist = list()
# filmusclesignif = data.frame()
# 
# for (i in unique(ctmeta$celltype)) {
#   
#   keep = ctmeta$celltype == i &
#     ctmeta$tissue %in% c("Bicep", "Thigh")
#   message("RUNNING DESEQ FOR ", i)
#   
#   count = ctcounts[, keep]
#   meta = ctmeta[keep,]
#   
#   print(table(meta$tissue))
#   print(meta[, c("sample", "tissue", "celltype", "individual")])
#   
#   cellcounts = UMct@meta.data %>%
#     filter(tissue %in% c("Bicep", "Thigh"),
#            celltype == i) %>%
#     dplyr::count(sample, tissue)
#   
#   
#   #Not all cell types are present in both tissue, so filter those cell types out
#   if(length(unique(meta$tissue)) < 2) {
#     message("Skipping ", i, " because only one of the tissues contains this cell type")
#     next
#   }
#   
#   #For this version we only want cell types for which there are 10 or more cells of that cell type present in both tissues
#   if(any(cellcounts$n < 10)) {
#     message("Skipping ", i, ": one or more pseudobulks have fewer than 10 cells of this cell type")
#     next
#   }
#   
#   #We also cannot compare when a cell type is only present in 1 individual since we can't account for effects of an individual
#   if(length(unique(meta$individual)) < 2) {
#     message("Skipping ", i, " because fewer than 2 individuals contain this cell type")
#     next
#   }
#   
#   #Even when the cell types are present in both tissue, they may not be present in every sample, so filter those out also
#   pairtab = table(meta$individual, meta$tissue)
#   
#   if(!all(c("Bicep", "Thigh") %in% colnames(pairtab)) ||
#      any(pairtab[, c("Bicep", "Thigh")] == 0)) {
#     message("Skipping", i, "because not all samples of both tissues contain the celltype")
#     next
#   }
#   
#   dds = DESeqDataSetFromMatrix(
#     countData = round(count),
#     colData = meta,
#     design = ~ individual + tissue
#   )
#   
#   dds = DESeq(dds)
#   
#   filmusclereslist[[i]] = results(dds, contrast = c("tissue", "Bicep", "Thigh"))
#   filmuscleplotlist[[i]] = makevolcano(filmusclereslist[[i]], paste0(i, ": bicep vs thigh"))
#   
#   safename = gsub("/", "_", i) #Cell types with a / in the name mess up ggsave
#   
#   ggsave(file.path(cttisDir, paste0("filtered_muscle_", safename, "_volcanoplot.pdf")), filmuscleplotlist[[i]], width = 7, height = 5)
#   
#   
#   #Also make a df for a barchart reflecting whether any HOX genes are significantly differentially expressed in any specific cell types
#   resdf = as.data.frame(filmusclereslist[[i]])
#   hox = grepl("^HOX", rownames(resdf))
#   
#   filmusclesignif = rbind(
#     filmusclesignif,
#     data.frame(
#       celltype = i,
#       comparison = "Bicep_vs_Thigh",
#       n_hox_sig = sum(
#         hox &
#           resdf$padj < 0.05,
#         na.rm = TRUE
#       ),
#       n_hox_logfc_sig = sum(
#         hox &
#           resdf$padj < 0.05 &
#           abs(resdf$log2FoldChange) > 1,
#         na.rm = TRUE
#       ),
#       n_hox_up = sum(
#         hox &
#           resdf$padj < 0.05 &
#           resdf$log2FoldChange > 1,
#         na.rm = TRUE
#       ),
#       n_hox_down = sum(
#         hox &
#           resdf$padj < 0.05 &
#           resdf$log2FoldChange < -1,
#         na.rm = TRUE
#       )
#     )
#     
#   )
#   
# }
# 
# filmusclehoxcount = ggplot(filmusclesignif, aes(x = reorder(celltype, n_hox_logfc_sig), y = n_hox_logfc_sig)) +
#   geom_col() +
#   coord_flip() +
#   theme_minimal()
# 
# ggsave(file.path(cttisDir, "filtered_muscle_SignificantHOX_DE_count_per_celltype.pdf"), filmusclehoxcount, width = 10, height = 8)
# 
# 
# 
# ############################################################################
# 
# 
# 
# 
# 
# 
# #########Composition and composition ratio plots -- to see whether tissue level DE might just be the result of different cell type abundances
# 
# message("====== Making composition plots ======")
# 
# #Get proportions
# comp = UMct@meta.data %>% 
#   dplyr::count(sample, tissue, celltype) %>%
#   group_by(sample) %>%
#   mutate(prop = n/sum(n)) %>%
#   ungroup()
# 
# #Order by cluster number
# comp$celltype = factor(
#   comp$celltype, 
#   levels = comp %>%
#     distinct(celltype) %>%
#     arrange(as.numeric(sub(" .*", "", celltype))) %>%
#     pull(celltype)
# )
# 
# #Plot
# compplot = ggplot(comp, aes(x = sample, y = prop, fill = celltype)) +
#   geom_col(position = "fill") +
#   facet_grid(~ tissue, scales = "free_x", space = "free_x") +
#   scale_y_continuous(labels = scales::percent_format()) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "right"
#   ) +
#   labs(
#     x = "Sample",
#     y = "Cell type proportion",
#     fill = "Cell type"
#   ) +
#   scale_fill_manual(values = ClusColours)
# 
# #Save
# ggsave(file.path(plottingDir, "sample_celltype_composition.pdf"), compplot, width = 14, height = 8)
# 
# 
# ######Ratios to check cell type relative abundances explicitly
# #Filters
# musclekeep = UMct@meta.data %>%
#   filter(tissue %in% c("Bicep", "Thigh")) %>%
#   dplyr::count(tissue, celltype) %>%
#   pivot_wider(
#     names_from = tissue,
#     values_from = n,
#     values_fill = 0
#   ) %>%
#   filter(Bicep >= 10 & Thigh >= 10) %>%
#   pull(celltype)
# 
# skinkeep = UMct@meta.data %>%
#   filter(tissue %in% c("Neck", "LowBack")) %>%
#   dplyr::count(tissue, celltype) %>%
#   pivot_wider(
#     names_from = tissue,
#     values_from = n,
#     values_fill = 0
#   ) %>%
#   filter(Neck >= 10 & LowBack >= 10) %>%
#   pull(celltype)
# 
# #Ratios
# comp_ratioM = UMct@meta.data %>%
#   filter(celltype %in% musclekeep) %>%
#   dplyr::count(sample, tissue, celltype) %>%
#   group_by(sample) %>%
#   mutate(prop = n / sum(n)) %>%
#   ungroup() %>%
#   group_by(tissue, celltype) %>%
#   summarise(mean_prop = mean(prop), .groups = "drop") %>%
#   filter(tissue %in% c("Bicep", "Thigh")) %>%
#   pivot_wider(names_from = tissue, values_from = mean_prop, values_fill = 0) %>%
#   mutate(ratio_Bicep_vs_Thigh = (Bicep + 1e-6) / (Thigh + 1e-6))
# 
# comp_ratioS = UMct@meta.data %>%
#   filter(celltype %in% skinkeep) %>%
#   dplyr::count(sample, tissue, celltype) %>%
#   group_by(sample) %>%
#   mutate(prop = n / sum(n)) %>%
#   ungroup() %>%
#   group_by(tissue, celltype) %>%
#   summarise(mean_prop = mean(prop), .groups = "drop") %>%
#   filter(tissue %in% c("Neck", "LowBack")) %>%
#   pivot_wider(names_from = tissue, values_from = mean_prop, values_fill = 0) %>%
#   mutate(ratio_Neck_vs_LowBack = (Neck + 1e-6) / (LowBack + 1e-6))
# 
# 
# #Plots
# muscleratio = ggplot(comp_ratioM, aes(x = reorder(celltype, ratio_Bicep_vs_Thigh),
#                        y = ratio_Bicep_vs_Thigh)) +
#   geom_col() +
#   coord_flip() +
#   theme_bw() +
#   labs(x = "Cell type", y = "Mean proportion ratio: Bicep / Thigh")
# 
# skinratio = ggplot(comp_ratioS, aes(x = reorder(celltype, ratio_Neck_vs_LowBack),
#                         y = ratio_Neck_vs_LowBack)) +
#   geom_col() +
#   coord_flip() +
#   theme_bw() +
#   labs(x = "Cell type", y = "Mean proportion ratio: Neck / LowBack")
# 
# #Save ratio DFs and plots
# comp_ratioM = comp_ratioM %>%
#   arrange(desc(ratio_Bicep_vs_Thigh))
# 
# comp_ratioS = comp_ratioS %>%
#   arrange(desc(ratio_Neck_vs_LowBack))
# 
# write.xlsx(comp_ratioM, file.path(plottingDir, "min10_celltype_ratios_MUSCLE.xlsx"), overwrite = TRUE)
# write.xlsx(comp_ratioS, file.path(plottingDir, "min10_celltype_ratios_SKIN.xlsx"), overwrite = TRUE)
# ggsave(file.path(plottingDir, "sample_celltype_ratios_MUSCLE.pdf"), muscleratio, width = 10, height = 8)
# ggsave(file.path(plottingDir, "sample_celltype_ratios_SKIN.pdf"), skinratio, width = 10, height = 8)
# 
# 
# ##### Table
# #Also adding this table to inspect the data and check for possible contamination
# celltable = as.data.frame.matrix(table(UM$sample, UM$celltype)) %>%
#   rownames_to_column("Sample")
# 
# write.xlsx(celltable, file.path(plottingDir, "celltype_counts_per_sample.xlsx"), overwrite = TRUE)
# 
# 
# 
# 
# 
# ########################## EXPLICIT PCAs #####################################
# 
# #Creating dirs
# pcaDir = file.path(plottingDir, "PCA")
# if (!dir.exists(pcaDir)) {
#   dir.create(pcaDir, recursive = TRUE)
# }
# 
# hoxpcamuscleDir = file.path(pcaDir, "celltype_HOX_PCAs_muscle")
# if (!dir.exists(hoxpcamuscleDir)) {
#   dir.create(hoxpcamuscleDir, recursive = TRUE)
# }
# 
# hoxpcaskinDir = file.path(pcaDir, "celltype_HOX_PCAs_skin")
# if (!dir.exists(hoxpcaskinDir)) {
#   dir.create(hoxpcaskinDir, recursive = TRUE)
# }
# 
# hvpcamuscleDir = file.path(pcaDir, "celltype_HVG_PCAs_muscle")
# if (!dir.exists(hvpcamuscleDir)) {
#   dir.create(hvpcamuscleDir, recursive = TRUE)
# }
# 
# hvpcaskinDir = file.path(pcaDir, "celltype_HVG_PCAs_skin")
# if (!dir.exists(hvpcaskinDir)) {
#   dir.create(hvpcaskinDir, recursive = TRUE)
# }
# 
# acpcamuscleDir = file.path(pcaDir, "Allcells_HVG_PCAs_muscle")
# if (!dir.exists(acpcamuscleDir)) {
#   dir.create(acpcamuscleDir, recursive = TRUE)
# }
# 
# acpcaskinDir = file.path(pcaDir, "Allcells_HVG_PCAs_skin")
# if (!dir.exists(acpcaskinDir)) {
#   dir.create(acpcaskinDir, recursive = TRUE)
# }
# 
# #Creating HVG PCA subsets per tissue type
# plan(sequential)
# UMpcaM = subset(UMct, subset = tissue %in% c("Bicep", "Thigh")) %>%
#   FindVariableFeatures() %>%
#   ScaleData()
# UMpcaM = RunPCA(UMpcaM, features = VariableFeatures(UMpcaM))
# 
# UMpcaS = subset(UMct, subset = tissue %in% c("Neck", "LowBack")) %>%
#   FindVariableFeatures() %>%
#   ScaleData()
# UMpcaS = RunPCA(UMpcaS, features = VariableFeatures(UMpcaS))
# plan(multicore, workers = 12)
# 
# ##########All cells - muscle
# 
# #Getting data
# pc.stdev.percentage = UMpcaM[["pca"]]@stdev / sum(UMpcaM[["pca"]]@stdev) * 100 #Find standard deviation for each PC
# cumu.percentage = cumsum(pc.stdev.percentage) #Find cumulative percentages for PCs
# co1 = which(cumu.percentage > 90 & pc.stdev.percentage < 5)[1] #Find PCs representing cumulative percent >90% and less than 5% associated with the single PC
# co1 #List PC
# co2 = sort(which((pc.stdev.percentage[1:length(pc.stdev.percentage) - 1] - pc.stdev.percentage[2:length(pc.stdev.percentage)]) > 0.01), 
#             decreasing = T)[1] #+ 1 #Find last PC where change in percent variation is more than 1% and add and extra PC to it
# 
# pcs = min(co1, co2) #Find the minimum PC from the 2 methods used above
# pcs < max(pcs) #List P
# pc.dims = 1:pcs #We will use these dimension later on
# 
# plot_df = data.frame(
#                       pc.stdev.percentage = pc.stdev.percentage, #Put PC values into dataframe for plotting
#                       cumu = cumu.percentage, 
#                       rank = 1:length(pc.stdev.percentage))
# #Cumulative plot 
# opt_PC_plot_muscle = ggplot(plot_df, aes(cumu.percentage, pc.stdev.percentage, label = rank, color = rank > pcs)) + 
#   geom_text() + geom_vline(xintercept = 90, color = "grey") + 
#   geom_hline(yintercept = min(pc.stdev.percentage[pc.stdev.percentage > 5]), color = "grey") + theme_bw() #Visualize PCs to use in elbow plot
# 
# ggsave(file.path(acpcamuscleDir, paste0("allcells_hvg_muscle_PCA_cumulative_plot.pdf")), opt_PC_plot_muscle, width = 10, height = 8)
# 
# #Optimal PC plot 
# elbowplotmuscle = ElbowPlot(UMpcaM, ndims = 50) + 
#   ggtitle(paste("Significant PCs", UMpcaM@meta.data[1,1], sep = " - ")) + 
#   geom_vline(xintercept = max(unlist(pc.dims)), linetype="dashed", color = "red", size=1) + 
#   theme(plot.title = element_text(size = 12))
# 
# ggsave(file.path(acpcamuscleDir, "allcells_hvg_muscle_PCA_elbowplot.pdf"), elbowplotmuscle, width = 10, height = 8)
# 
# #Create all PCA plots for the PCs that cumulatively explain >90% of variance and are less than 5% associated with the single PC
# for(i in seq(1, pcs, by = 2)) {
#   
#   p = DimPlot(UMpcaM, reduction = "pca", group.by = "tissue", raster = FALSE, dims = c(i, i + 1))
#   
#   ggsave(file.path(acpcamuscleDir, paste0("allcells_hvg_muscle_PCA_PC", i, "_PC", i + 1, ".pdf")), p, width = 10, height = 8)
#   
# }
# 
# #Get genes positively and negatively contributing to the PCs so we can say something about genes driving variance between tissues
# muscleload = Loadings(UMpcaM, reduction = "pca")
# toppcs_muscle = lapply(seq_len(ncol(muscleload)), function(pc) {
#   posgenes = names(sort(muscleload[, pc], decreasing = TRUE))[1:25]
#   neggenes = names(sort(muscleload[, pc], decreasing = FALSE))[1:25]
#   data.frame(
#     PC = paste0("PC", pc),
#     t(posgenes),
#     t(neggenes),
#     check.names = FALSE
#   )
# })
# 
# toppcs_muscle = bind_rows(toppcs_muscle)
# colnames(toppcs_muscle) = c("PC", paste0("Positive_", 1:25), paste0("Negative_", 1:25))
# 
# write.xlsx(toppcs_muscle, file.path(acpcamuscleDir, "allcells_hvg_muscle_PC_driver_genes.xlsx"), overwrite = TRUE)
# 
# 
# 
# ############All cells - skin
# #Getting data
# pc.stdev.percentage2 = UMpcaS[["pca"]]@stdev / sum(UMpcaS[["pca"]]@stdev) * 100 #Find standard deviation for each PC
# cumu.percentage2 = cumsum(pc.stdev.percentage2) #Find cumulative percentages for PCs
# co3 = which(cumu.percentage2 > 90 & pc.stdev.percentage2 < 5)[1] #Find PCs representing cumulative percent >90% and less than 5% associated with the single PC
# co3 #List PC
# co4 = sort(which((pc.stdev.percentage2[1:length(pc.stdev.percentage2) - 1] - pc.stdev.percentage2[2:length(pc.stdev.percentage2)]) > 0.01), 
#            decreasing = T)[1] #+ 1 #Find last PC where change in percent variation is more than 1% and add and extra PC to it
# 
# pcs2 = min(co3, co4) #Find the minimum PC from the 2 methods used above
# pcs2 < max(pcs2) #List P
# pc.dims2 = 1:pcs2 #We will use these dimension later on
# 
# plot_df2 = data.frame(
#   pc.stdev.percentage = pc.stdev.percentage2, #Put PC values into dataframe for plotting
#   cumu = cumu.percentage2, 
#   rank = 1:length(pc.stdev.percentage2))
# #Cumulative plot 
# opt_PC_plot_skin = ggplot(plot_df2, aes(cumu.percentage, pc.stdev.percentage, label = rank, color = rank > pcs)) + 
#   geom_text() + geom_vline(xintercept = 90, color = "grey") + 
#   geom_hline(yintercept = min(pc.stdev.percentage2[pc.stdev.percentage2 > 5]), color = "grey") + theme_bw() #Visualize PCs to use in elbow plot
# 
# ggsave(file.path(acpcaskinDir, paste0("allcells_hvg_skin_PCA_cumulative_plot.pdf")), opt_PC_plot_skin, width = 10, height = 8)
# 
# #Optimal PC plot 
# elbowplotskin = ElbowPlot(UMpcaS, ndims = 50) + 
#   ggtitle(paste("Significant PCs", UMpcaS@meta.data[1,1], sep = " - ")) + 
#   geom_vline(xintercept = max(unlist(pc.dims2)), linetype="dashed", color = "red", size=1) + 
#   theme(plot.title = element_text(size = 12))
# 
# ggsave(file.path(acpcaskinDir, paste0("allcells_hvg_skin_PCA_elbowplot.pdf")), elbowplotskin, width = 10, height = 8)
# 
# #Create all PCA plots for the PCs that cumulatively explain >90% of variance and are less than 5% associated with the single PC
# for(i in seq(1, ncol(Embeddings(UMpcaS, "pca")), by = 2)) {
#   
#   p = DimPlot(UMpcaS, reduction = "pca", group.by = "tissue", raster = FALSE, dims = c(i, i + 1))
#   
#   ggsave(file.path(acpcaskinDir, paste0("allcells_hvg_skin_PCA_PC", i, "_PC", i + 1, ".pdf")), p, width = 10, height = 8)
#   
# }
# 
# #Get genes positively and negatively contributing to the PCs so we can say something about genes driving variance between tissues
# skinload = Loadings(UMpcaS, reduction = "pca")
# toppcs_skin = lapply(seq_len(ncol(skinload)), function(pc) {
#   posgenes = names(sort(skinload[, pc], decreasing = TRUE))[1:25]
#   neggenes = names(sort(skinload[, pc], decreasing = FALSE))[1:25]
#   data.frame(
#     PC = paste0("PC", pc),
#     t(posgenes),
#     t(neggenes),
#     check.names = FALSE
#   )
# })
# 
# toppcs_skin = bind_rows(toppcs_skin)
# colnames(toppcs_skin) = c("PC", paste0("Positive_", 1:25), paste0("Negative_", 1:25))
# 
# write.xlsx(toppcs_skin, file.path(acpcaskinDir, "allcells_hvg_skin_PC_driver_genes.xlsx"), overwrite = TRUE)
# 
# 
# 
# 
# 
# #################################################
# #STATISTICAL TESTS
# #################################################
# 
# 
# ################################ PERMANOVA #################################
# #This section is predominantly used for PERMANOVA analysis, but for the highly variable gene loops PCA plots and files are also created to avoid pseudobulking the same type of data twice by doing it both
# #in the PCA section above and here.
# 
# 
# ########### PERMANOVA FOR HOX GENES SPECIFICALLY ##############
# plan(sequential)
# profpseudo = AggregateExpression(UMct, assays = "RNA", return.seurat = TRUE, group.by = c("sample", "tissue", "celltype"))
# plan(multicore, workers = 12)
# 
# #Add informative metadata
# profparts = do.call(rbind, strsplit(Cells(profpseudo), "_"))
# profpseudo$sample = profparts[,1]
# profpseudo$tissue = profparts[,2]
# profpseudo$celltype = profparts[,3]
# Idents(profpseudo) = profpseudo$celltype
# meta_map2 = unique(UMct@meta.data[, c("sample", "individual", "protocol")])
# rownames(meta_map2) = meta_map2$sample
# profpseudo$individual = meta_map2[profpseudo$sample, "individual"]
# profpseudo$protocol = meta_map2[profpseudo$sample, "protocol"]
# profcounts = GetAssayData(profpseudo, assay = "RNA", layer = "counts")
# profmeta = profpseudo@meta.data[, c("sample", "tissue", "celltype", "individual", "protocol")]
# rownames(profmeta) = colnames(profcounts)
# 
# #For skin ################################
# permaskin = data.frame()
# permalistskin = list()
# for(i in unique(profmeta$celltype)) {
#   
#   keep = profmeta$celltype == i &
#     profmeta$tissue %in% c("Neck", "LowBack")
#   message("RUNNING PERMANOVA FOR ", i)
#   
#   count = profcounts[, keep]
#   meta = profmeta[keep,]
#   
#   pairtab = table(meta$individual, meta$tissue)
#   
#   if(!all(c("Neck", "LowBack") %in% colnames(pairtab))) {
#     message("Skipping ", i, ": missing Neck or LowBack")
#     next
#   }
#   
#   if(length(unique(meta$individual)) < 2 ||
#      any(pairtab[ , c("Neck", "LowBack"), drop = FALSE] == 0)) {
#     message("Skipping ", i, ": not enough paired individuals")
#     next
#   }
#   
#   meta = droplevels(meta)
#   meta$tissue = factor(meta$tissue, levels = c("Neck", "LowBack"))
#   meta$individual = factor(meta$individual)
#   
#   ###
#   #This section normalizes the count data to avoid highly expressing HOX genes dominating the statistics calculation
#   #Better because adonis2 works by calculating distances
#   dds = DESeqDataSetFromMatrix(
#     countData = round(count),
#     colData = meta,
#     design = ~ tissue
#   )
#   vsd = varianceStabilizingTransformation(dds, blind = TRUE) #VST - Variance stabilizing transformation
#   hoxkeep = grepl("^HOX", rownames(assay(vsd)))
#   matrix = t(assay(vsd)[hoxkeep, ])
#   message(i, ":", ncol(matrix), " HOX genes")
#   ###
#   
#   if(ncol(matrix) < 2 || nrow(matrix) < 3) {
#     message("Skipping: ", i,  " too few HOX genes or samples")
#     next
#   }
# 
#   #Include PCA plots 
#   pca = prcomp(matrix, center = TRUE, scale. = FALSE)
#   scores = as.data.frame(pca$x)
#   scores$tissue = meta$tissue
#   scores$individual = meta$individual
#   percentVar = 100 * pca$sdev^2 / sum(pca$sdev^2)
#   maxpc = ncol(pca$x)
#   
#   for(j in seq(1, maxpc - 1, by = 2)){
#     
#     p = ggplot(scores, aes_string(x = paste0("PC", j),
#                                   y = paste0("PC", j + 1),
#                                   colour = "tissue",
#                                   shape = "individual")) +
#       geom_point(size = 3) +
#       theme_bw() +
#       labs(
#         x = paste0("PC", j, " (", round(percentVar[j], 1), "%)"),
#         y = paste0("PC", j + 1, " (", round(percentVar[j + 1], 1), "%)"),
#         title = paste0("PC", j, " vs PC", j + 1),
#         subtitle = i
#       )
#     
#     ggsave(file.path(hoxpcaskinDir, 
#                      paste0("PCA_hox_skin_PC", j, "_PC", j + 1, "_", make.names(i),".pdf")), 
#            p, width = 10, height = 8)
#   }
#   
#   #Store PC genes
#   loadings = pca$rotation
#   pcgenetable = lapply(seq_len(ncol(loadings)), function(pc) {
#     
#     pos = sort(loadings[, pc], decreasing = TRUE)[1:25]
#     neg = sort(loadings[, pc], decreasing = FALSE)[1:25]
#     data.frame(
#       PC = paste0("PC_", pc),
#       Positive_gene = names(pos),
#       Positive_loading = as.numeric(pos),
#       Negative_gene = names(neg),
#       Negative_loading = as.numeric(neg)
#     )
#   }) %>%
#     dplyr::bind_rows()
#   permalistskin[[i]] = pcgenetable
#   
#   #PERMANOVA
#   ad = adonis2(matrix ~ tissue,
#           data = meta,
#           strata = meta$individual,
#           method = "euclidean",
#           permutations = 9999)
#   
#   permaskin = rbind(permaskin, data.frame(celltype = i, R2 = ad$R2[1], F = ad$F[1], pval = ad$`Pr(>F)`[1]))
# 
# }
# 
# permaskin = permaskin %>%
#   arrange(desc(R2))
# 
# write.xlsx(permaskin, file.path(plottingDir, "PERMANOVA_HOXGENES_SKIN.xlsx"), overwrite = TRUE)
# names(permalistskin) = make.unique(
#   substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", names(permalistskin)), 1, 31)
# )
# write.xlsx(permalistskin, file.path(hoxpcaskinDir, "celltype_HOX_skin_PCA_genes_per_PC.xlsx"), overwrite = TRUE)  
# 
# 
# #For muscle
# permamusc = data.frame()
# permalistmuscle = list()
# for(i in unique(profmeta$celltype)) {
#   
#   keep = profmeta$celltype == i &
#     profmeta$tissue %in% c("Bicep", "Thigh")
#   message("RUNNING PERMANOVA FOR ", i)
#   
#   count = profcounts[, keep]
#   meta = profmeta[keep,]
#   
#   pairtab = table(meta$individual, meta$tissue)
#   
#   if(!all(c("Bicep", "Thigh") %in% colnames(pairtab))) {
#     message("Skipping ", i, ": missing Bicep or Thigh")
#     next
#   }
#   
#   if(length(unique(meta$individual)) < 2 ||
#      any(pairtab[ , c("Bicep", "Thigh"), drop = FALSE] == 0)) {
#     message("Skipping ", i, ": not enough paired individuals")
#     next
#   }
#   
#   meta = droplevels(meta)
#   meta$tissue = factor(meta$tissue, levels = c("Bicep", "Thigh"))
#   meta$individual = factor(meta$individual)
#   
#   ###
#   #This section normalizes the count data to avoid highly expressing HOX genes dominating the statistics calculation
#   #Better because adonis2 works by calculating distances
#   dds = DESeqDataSetFromMatrix(
#     countData = round(count),
#     colData = meta,
#     design = ~ tissue
#   )
#   vsd = varianceStabilizingTransformation(dds, blind = TRUE) #VST - Variance stabilizing transformation
#   hoxkeep = grepl("^HOX", rownames(assay(vsd)))
#   matrix = t(assay(vsd)[hoxkeep, ])
#   message(i, ":", ncol(matrix), " HOX genes")
#   ###
#   
#   if(ncol(matrix) < 2 || nrow(matrix) < 3) {
#     message("Skipping: ", i,  " too few HOX genes or samples")
#     next
#   }
# 
#   #Include PCA plots 
#   pca = prcomp(matrix, center = TRUE, scale. = FALSE)
#   scores = as.data.frame(pca$x)
#   scores$tissue = meta$tissue
#   scores$individual = meta$individual
#   percentVar = 100 * pca$sdev^2 / sum(pca$sdev^2)
#   maxpc = ncol(pca$x)
#   
#   for(j in seq(1, maxpc - 1, by = 2)){
#     
#     p = ggplot(scores, aes_string(x = paste0("PC", j),
#                                   y = paste0("PC", j + 1),
#                                   colour = "tissue",
#                                   shape = "individual")) +
#       geom_point(size = 3) +
#       theme_bw() +
#       labs(
#         x = paste0("PC", j, " (", round(percentVar[j], 1), "%)"),
#         y = paste0("PC", j + 1, " (", round(percentVar[j + 1], 1), "%)"),
#         title = paste0("PC", j, " vs PC", j + 1),
#         subtitle = i
#       )
#     
#     ggsave(file.path(hoxpcamuscleDir, 
#                      paste0("PCA_hox_muscle_PC", j, "_PC", j + 1, "_", make.names(i),".pdf")), 
#            p, width = 10, height = 8)
#   }
#   
#   #Store PC genes
#   loadings = pca$rotation
#   pcgenetable = lapply(seq_len(ncol(loadings)), function(pc) {
#     
#     pos = sort(loadings[, pc], decreasing = TRUE)[1:25]
#     neg = sort(loadings[, pc], decreasing = FALSE)[1:25]
#     data.frame(
#       PC = paste0("PC_", pc),
#       Positive_gene = names(pos),
#       Positive_loading = as.numeric(pos),
#       Negative_gene = names(neg),
#       Negative_loading = as.numeric(neg)
#     )
#   }) %>%
#     dplyr::bind_rows()
#   permalistmuscle[[i]] = pcgenetable
#   
#   #PERMANOVA
#   ad = adonis2(matrix ~ tissue,
#                data = meta,
#                strata = meta$individual,
#                method = "euclidean",
#                permutations = 9999)
#   
#   permamusc = rbind(permamusc, data.frame(celltype = i, R2 = ad$R2[1], F = ad$F[1], pval = ad$`Pr(>F)`[1]))
# 
# }
# 
# permamusc = permamusc %>%
#   arrange(desc(R2))
#   
# write.xlsx(permamusc, file.path(plottingDir, "PERMANOVA_HOXGENES_MUSCLE.xlsx"), overwrite = TRUE)
# names(permalistmuscle) = make.unique(
#   substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", names(permalistmuscle)), 1, 31)
# )
# write.xlsx(permalistmuscle, file.path(hoxpcamuscleDir, "celltype_HOX_muscle_PCA_genes_per_PC.xlsx"), overwrite = TRUE)
# 
# 
# ########### PERMANOVA FOR ALL VARIABLE GENES PER CELL TYPE ##############
# plan(sequential)
# profpseudo3 = AggregateExpression(UMct, assays = "RNA", return.seurat = TRUE, group.by = c("sample", "tissue", "celltype"))
# plan(multicore, workers = 12)
# 
# #Define variable gene set, NEEDS TO BE DONE GLOBALLY SHOULD NOT BE DONE PER SET 
# hvgenes = VariableFeatures(UMct)
# 
# #Add informative metadata
# profparts3 = do.call(rbind, strsplit(Cells(profpseudo3), "_"))
# profpseudo3$sample = profparts3[,1]
# profpseudo3$tissue = profparts3[,2]
# profpseudo3$celltype = profparts3[,3]
# Idents(profpseudo3) = profpseudo3$celltype
# meta_map3 = unique(UMct@meta.data[, c("sample", "individual", "protocol")])
# rownames(meta_map3) = meta_map3$sample
# profpseudo3$individual = meta_map3[profpseudo3$sample, "individual"]
# profpseudo3$protocol = meta_map3[profpseudo3$sample, "protocol"]
# profcounts3 = GetAssayData(profpseudo3, assay = "RNA", layer = "counts")
# profmeta3 = profpseudo3@meta.data[, c("sample", "tissue", "celltype", "individual", "protocol")]
# rownames(profmeta3) = colnames(profcounts3)
# 
# #For skin ################################
# varyskin = data.frame()
# varylistskin = list()
# for(i in unique(profmeta3$celltype)) {
#   
#   keep = profmeta3$celltype == i &
#     profmeta3$tissue %in% c("Neck", "LowBack")
#   message("RUNNING PERMANOVA FOR ", i)
#   
#   count = profcounts3[, keep]
#   meta = profmeta3[keep,]
#   
#   pairtab = table(meta$individual, meta$tissue)
#   
#   if(!all(c("Neck", "LowBack") %in% colnames(pairtab))) {
#     message("Skipping ", i, ": missing Neck or LowBack")
#     next
#   }
#   
#   if(length(unique(meta$individual)) < 2 ||
#      any(pairtab[ , c("Neck", "LowBack"), drop = FALSE] == 0)) {
#     message("Skipping ", i, ": not enough paired individuals")
#     next
#   }
#   
#   meta = droplevels(meta)
#   meta$tissue = factor(meta$tissue, levels = c("Neck", "LowBack"))
#   meta$individual = factor(meta$individual)
#   
#   ###
#   #This section normalizes the count data to avoid highly expressing HOX genes dominating the statistics calculation
#   #Better because adonis2 works by calculating distances
#   dds = DESeqDataSetFromMatrix(
#     countData = round(count),
#     colData = meta,
#     design = ~ tissue
#   )
#   vsd = varianceStabilizingTransformation(dds, blind = TRUE) #VST - Variance stabilizing transformation
#   matrix = t(assay(vsd)[hvgenes, ])
#   message(i, ":", ncol(matrix), " variable genes")
#   ###
#   
#   if(ncol(matrix) < 2 || nrow(matrix) < 3) {
#     message("Skipping: ", i,  " too few variable genes or samples")
#     next
#   }
#   
#   #Include PCA plots 
#   pca = prcomp(matrix, center = TRUE, scale. = FALSE)
#   scores = as.data.frame(pca$x)
#   scores$tissue = meta$tissue
#   scores$individual = meta$individual
#   percentVar = 100 * pca$sdev^2 / sum(pca$sdev^2)
#   maxpc = ncol(pca$x)
#   
#   for(j in seq(1, maxpc - 1, by = 2)){
#     
#     p = ggplot(scores, aes_string(x = paste0("PC", j),
#                                   y = paste0("PC", j + 1),
#                                   colour = "tissue",
#                                   shape = "individual")) +
#       geom_point(size = 3) +
#       theme_bw() +
#       labs(
#         x = paste0("PC", j, " (", round(percentVar[j], 1), "%)"),
#         y = paste0("PC", j + 1, " (", round(percentVar[j + 1], 1), "%)"),
#         title = paste0("PC", j, " vs PC", j + 1),
#         subtitle = i
#       )
#     
#     ggsave(file.path(hvpcaskinDir, 
#                      paste0("PCA_hvg_skin_PC", j, "_PC", j + 1, "_", make.names(i),".pdf")), 
#            p, width = 10, height = 8)
#   }
#   
#   #Store PC genes
#   loadings = pca$rotation
#   pcgenetable = lapply(seq_len(ncol(loadings)), function(pc) {
#     
#     pos = sort(loadings[, pc], decreasing = TRUE)[1:25]
#     neg = sort(loadings[, pc], decreasing = FALSE)[1:25]
#     data.frame(
#       PC = paste0("PC_", pc),
#       Positive_gene = names(pos),
#       Positive_loading = as.numeric(pos),
#       Negative_gene = names(neg),
#       Negative_loading = as.numeric(neg)
#     )
#   }) %>%
#     dplyr::bind_rows()
#   varylistskin[[i]] = pcgenetable
#   
#   #PERMANOVA
#   ad = adonis2(matrix ~ tissue,
#                data = meta,
#                strata = meta$individual,
#                method = "euclidean",
#                permutations = 9999)
#   
#   varyskin = rbind(varyskin, data.frame(celltype = i, R2 = ad$R2[1], F = ad$F[1], pval = ad$`Pr(>F)`[1]))
# 
# }
# 
# varyskin = varyskin %>%
#   arrange(desc(R2))
# 
# write.xlsx(varyskin, file.path(plottingDir, "PERMANOVA_VARIABLEGENES_SKIN.xlsx"), overwrite = TRUE)
# names(varylistskin) = make.unique(
#   substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", names(varylistskin)), 1, 31)
# )
# write.xlsx(varylistskin, file.path(hvpcaskinDir, "celltype_HVG_skin_PCA_genes_per_PC.xlsx"), overwrite = TRUE)
# 
# 
# 
# #For muscle
# varymusc = data.frame()
# varylistmuscle = list()
# for(i in unique(profmeta3$celltype)) {
#   
#   keep = profmeta3$celltype == i &
#     profmeta3$tissue %in% c("Bicep", "Thigh")
#   message("RUNNING PERMANOVA FOR ", i)
#   
#   count = profcounts3[, keep]
#   meta = profmeta3[keep,]
#   
#   pairtab = table(meta$individual, meta$tissue)
#   
#   if(!all(c("Bicep", "Thigh") %in% colnames(pairtab))) {
#     message("Skipping ", i, ": missing Bicep or Thigh")
#     next
#   }
#   
#   if(length(unique(meta$individual)) < 2 ||
#      any(pairtab[ , c("Bicep", "Thigh"), drop = FALSE] == 0)) {
#     message("Skipping ", i, ": not enough paired individuals")
#     next
#   }
#   
#   meta = droplevels(meta)
#   meta$tissue = factor(meta$tissue, levels = c("Bicep", "Thigh"))
#   meta$individual = factor(meta$individual)
#   
#   ###
#   #This section normalizes the count data to avoid highly expressing HOX genes dominating the statistics calculation
#   #Better because adonis2 works by calculating distances
#   dds = DESeqDataSetFromMatrix(
#     countData = round(count),
#     colData = meta,
#     design = ~ tissue
#   )
#   vsd = varianceStabilizingTransformation(dds, blind = TRUE) #VST - Variance stabilizing transformation
#   matrix = t(assay(vsd)[hvgenes, ])
#   message(i, ":", ncol(matrix), " variable genes genes")
#   ###
#   
#   if(ncol(matrix) < 2 || nrow(matrix) < 3) {
#     message("Skipping: ", i,  " too few variable genes or samples")
#     next
#   }
#   
#   #Include PCA plots 
#   pca = prcomp(matrix, center = TRUE, scale. = FALSE)
#   scores = as.data.frame(pca$x)
#   scores$tissue = meta$tissue
#   scores$individual = meta$individual
#   percentVar = 100 * pca$sdev^2 / sum(pca$sdev^2)
#   maxpc = ncol(pca$x)
#   
#   for(j in seq(1, maxpc - 1, by = 2)){
#     
#     p = ggplot(scores, aes_string(x = paste0("PC", j),
#                                   y = paste0("PC", j + 1),
#                                   colour = "tissue",
#                                   shape = "individual")) +
#       geom_point(size = 3) +
#       theme_bw() +
#       labs(
#         x = paste0("PC", j, " (", round(percentVar[j], 1), "%)"),
#         y = paste0("PC", j + 1, " (", round(percentVar[j + 1], 1), "%)"),
#         title = paste0("PC", j, " vs PC", j + 1),
#         subtitle = i
#       )
#     
#     ggsave(file.path(hvpcamuscleDir, 
#                      paste0("PCA_hvg_muscle_PC", j, "_PC", j + 1, "_", make.names(i),".pdf")), 
#            p, width = 10, height = 8)
#   }
#   
#   #Store PC genes
#   loadings = pca$rotation
#   pcgenetable = lapply(seq_len(ncol(loadings)), function(pc) {
#     
#     pos = sort(loadings[, pc], decreasing = TRUE)[1:25]
#     neg = sort(loadings[, pc], decreasing = FALSE)[1:25]
#     data.frame(
#       PC = paste0("PC_", pc),
#       Positive_gene = names(pos),
#       Positive_loading = as.numeric(pos),
#       Negative_gene = names(neg),
#       Negative_loading = as.numeric(neg)
#     )
#   }) %>%
#     dplyr::bind_rows()
#   varylistmuscle[[i]] = pcgenetable
#   
#   #PERMANOVA
#   ad = adonis2(matrix ~ tissue,
#                data = meta,
#                strata = meta$individual,
#                method = "euclidean",
#                permutations = 9999)
#   
#   varymusc = rbind(varymusc, data.frame(celltype = i, R2 = ad$R2[1], F = ad$F[1], pval = ad$`Pr(>F)`[1]))
#   
# }
# 
# varymusc = varymusc %>%
#   arrange(desc(R2))
# 
# write.xlsx(varymusc, file.path(plottingDir, "PERMANOVA_VARIABLEGENES_MUSCLE.xlsx"), overwrite = TRUE)
# names(varylistmuscle) = make.unique(
#   substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", names(varylistmuscle)), 1, 31)
# )
# write.xlsx(varylistmuscle, file.path(hvpcamuscleDir, "celltype_HVG_muscle_PCA_genes_per_PC.xlsx"), overwrite = TRUE)
# 
# 
# ########### PERMANOVA FOR ALL GENES PER CELL TYPE ##############
# plan(sequential)
# profpseudo4 = AggregateExpression(UMct, assays = "RNA", return.seurat = TRUE, group.by = c("sample", "tissue", "celltype"))
# plan(multicore, workers = 12)
# 
# #Define variable gene set, NEEDS TO BE DONE GLOBALLY SHOULD NOT BE DONE PER SET 
# allgenes = rownames(UMct)
# 
# #Add informative metadata
# profparts4 = do.call(rbind, strsplit(Cells(profpseudo4), "_"))
# profpseudo4$sample = profparts4[,1]
# profpseudo4$tissue = profparts4[,2]
# profpseudo4$celltype = profparts4[,3]
# Idents(profpseudo4) = profpseudo4$celltype
# meta_map4 = unique(UMct@meta.data[, c("sample", "individual", "protocol")])
# rownames(meta_map4) = meta_map4$sample
# profpseudo4$individual = meta_map4[profpseudo4$sample, "individual"]
# profpseudo4$protocol = meta_map4[profpseudo4$sample, "protocol"]
# profcounts4 = GetAssayData(profpseudo4, assay = "RNA", layer = "counts")
# profmeta4 = profpseudo4@meta.data[, c("sample", "tissue", "celltype", "individual", "protocol")]
# rownames(profmeta4) = colnames(profcounts4)
# 
# #For skin ################################
# allskin = data.frame()
# for(i in unique(profmeta4$celltype)) {
#   
#   keep = profmeta4$celltype == i &
#     profmeta4$tissue %in% c("Neck", "LowBack")
#   message("RUNNING PERMANOVA FOR ", i)
#   
#   count = profcounts4[, keep]
#   meta = profmeta4[keep,]
#   
#   pairtab = table(meta$individual, meta$tissue)
#   
#   if(!all(c("Neck", "LowBack") %in% colnames(pairtab))) {
#     message("Skipping ", i, ": missing Neck or LowBack")
#     next
#   }
#   
#   if(length(unique(meta$individual)) < 2 ||
#      any(pairtab[ , c("Neck", "LowBack"), drop = FALSE] == 0)) {
#     message("Skipping ", i, ": not enough paired individuals")
#     next
#   }
#   
#   meta = droplevels(meta)
#   meta$tissue = factor(meta$tissue, levels = c("Neck", "LowBack"))
#   meta$individual = factor(meta$individual)
#   
#   ###
#   #This section normalizes the count data to avoid highly expressing HOX genes dominating the statistics calculation
#   #Better because adonis2 works by calculating distances
#   dds = DESeqDataSetFromMatrix(
#     countData = round(count),
#     colData = meta,
#     design = ~ tissue
#   )
#   vsd = varianceStabilizingTransformation(dds, blind = TRUE) #VST - Variance stabilizing transformation
#   
#   matrix = t(assay(vsd))
#   message(i, ":", ncol(matrix), " variable genes")
#   ###
#   
#   if(ncol(matrix) < 2 || nrow(matrix) < 3) {
#     message("Skipping: ", i,  " too few variable genes or samples")
#     next
#   }
#   
#   ad = adonis2(matrix ~ tissue,
#                data = meta,
#                strata = meta$individual,
#                method = "euclidean",
#                permutations = 9999)
#   
#   allskin = rbind(allskin, data.frame(celltype = i, R2 = ad$R2[1], F = ad$F[1], pval = ad$`Pr(>F)`[1]))
#    
# }
# 
# allskin = allskin %>%
#   arrange(desc(R2))
# 
# write.xlsx(allskin, file.path(plottingDir, "PERMANOVA_ALLGENES_SKIN.xlsx"), overwrite = TRUE)
# 
# 
# 
# #For muscle
# allmusc = data.frame()
# for(i in unique(profmeta4$celltype)) {
#   
#   keep = profmeta4$celltype == i &
#     profmeta4$tissue %in% c("Bicep", "Thigh")
#   message("RUNNING PERMANOVA FOR ", i)
#   
#   count = profcounts4[, keep]
#   meta = profmeta4[keep,]
#   
#   pairtab = table(meta$individual, meta$tissue)
#   
#   if(!all(c("Bicep", "Thigh") %in% colnames(pairtab))) {
#     message("Skipping ", i, ": missing Bicep or Thigh")
#     next
#   }
#   
#   if(length(unique(meta$individual)) < 2 ||
#      any(pairtab[ , c("Bicep", "Thigh"), drop = FALSE] == 0)) {
#     message("Skipping ", i, ": not enough paired individuals")
#     next
#   }
#   
#   meta = droplevels(meta)
#   meta$tissue = factor(meta$tissue, levels = c("Bicep", "Thigh"))
#   meta$individual = factor(meta$individual)
#   
#   ###
#   #This section normalizes the count data to avoid highly expressing HOX genes dominating the statistics calculation
#   #Better because adonis2 works by calculating distances
#   dds = DESeqDataSetFromMatrix(
#     countData = round(count),
#     colData = meta,
#     design = ~ tissue
#   )
#   vsd = varianceStabilizingTransformation(dds, blind = TRUE) #VST - Variance stabilizing transformation
#   
#   matrix = t(assay(vsd))
#   message(i, ":", ncol(matrix), " variable genes")
#   ###
#   
#   if(ncol(matrix) < 2 || nrow(matrix) < 3) {
#     message("Skipping: ", i,  " too few variable genes or samples")
#     next
#   }
#   
#   ad = adonis2(matrix ~ tissue,
#                data = meta,
#                strata = meta$individual,
#                method = "euclidean",
#                permutations = 9999)
#   
#   allmusc = rbind(allmusc, data.frame(celltype = i, R2 = ad$R2[1], F = ad$F[1], pval = ad$`Pr(>F)`[1]))
#   
# }
# 
# allmusc = allmusc %>%
#   arrange(desc(R2))
# 
# write.xlsx(allmusc, file.path(plottingDir, "PERMANOVA_ALLGENES_MUSCLE.xlsx"), overwrite = TRUE)
# 
# 
# 
# 
# ########################### Explicit stats to confirm support of collinearity principle
# 
# #### SKIN
# colskin = data.frame()
# heatlist.skin = list()
# #Will just use profpseudo and profmeta again since they have the structure I want
# for( i in unique(profmeta$celltype)) {
#   
#   keep = profmeta$celltype == i &
#     profmeta$tissue %in% c("Neck", "LowBack")
#   
#   count = profcounts[, keep]
#   meta = profmeta[keep,]
#   
#   pairtab = table(meta$individual, meta$tissue)
#   
#   if(!all(c("Neck", "LowBack") %in% colnames(pairtab))) next
#   
#   if(length(unique(meta$individual)) < 2 ||
#      any(pairtab[ , c("Neck", "LowBack"), drop = FALSE] == 0)) next
#   
#   meta = droplevels(meta)
#   meta$tissue = factor(meta$tissue, levels = c("Neck", "LowBack"))
#   meta$individual = factor(meta$individual)
#   
#   dds = DESeqDataSetFromMatrix(
#     countData = round(count),
#     colData = meta,
#     design = ~ tissue
#   )
#   vsd = varianceStabilizingTransformation(dds, blind = TRUE) #VST - Variance stabilizing transformation
#   hoxkeep = grepl("^HOX", rownames(assay(vsd)))
#   matrix = assay(vsd)[hoxkeep, ]
#   hoxmat = matrix[grepl("\\d+$", rownames(matrix)), , drop = FALSE] #Keep only HOX genes with numbers at the end of their names
#   
#   if(nrow(hoxmat) < 3) next
#   
#   #Make long
#   long = as.data.frame(as.matrix(hoxmat))
#   long$gene = rownames(long)
#   long = long %>%
#     pivot_longer(cols = -gene, names_to = "sample_id", values_to = "expression")
#   
#   meta2 = meta %>%
#     mutate(sample_id = rownames(meta)) %>%
#     dplyr::select(sample_id, sample, individual, tissue)
#   
#   long = left_join(long, meta2, by = "sample_id")
#   
#   #Calculate anterior-posterior delta
#   wide = long %>%
#     dplyr::select(gene, individual, tissue, expression) %>%
#     pivot_wider(names_from = tissue, values_from = expression) %>%
#     filter(!is.na(Neck) & !is.na(LowBack)) %>%
#     mutate(delta = LowBack - Neck,
#            paralog = as.numeric(str_extract(gene, "\\d+$")),
#            celltype = i) %>%
#     filter(!is.na(paralog))
#   
#   #Store info for heatmap
#   heatlist.skin[[i]] = wide %>%
#     group_by(gene) %>%
#     summarize(mean_delta = mean(delta), .groups = "drop") %>%
#     mutate(celltype = i)
#   
#   #Test collinearity principle
#   fit = lmer(delta ~ paralog + (1|individual), data = wide)
#   print(summary(fit))
#   
#   #Spearman on mean delta per gene
#   mean_delta = wide %>%
#     group_by(gene, paralog) %>%
#     summarize(mean_delta = mean(delta, na.rm = TRUE), .groups = "drop")
#   
#   cor_res = cor.test(mean_delta$paralog, mean_delta$mean_delta, method = "spearman", exact = FALSE)
#   
#   colskin = rbind(
#     colskin,
#     data.frame(
#       celltype = i,
#       spearman_rho = unname(cor_res$estimate),
#       spearman_pval = cor_res$p.value,
#       lmer_slope = summary(fit)$coefficients["paralog", "Estimate"],
#       lmer_pval = summary(fit)$coefficients["paralog", "Pr(>|t|)"],
#       n_hox_genes = length(unique(wide$gene)),
#       n_individuals = length(unique(wide$individual))
#     ))
#   
# }
# 
# #Save numbers
# colskin = colskin %>% arrange(desc(lmer_slope))
# write.xlsx(colskin, file.path(plottingDir, "COLLINEARITY_HOXGENES_SKIN.xlsx"), overwrite = TRUE)
# 
# ### Make heatmap
# #Setting order
# skinheatdf = bind_rows(heatlist.skin) %>%
#   mutate(cluster = str_extract(gene, "HOX[ABCD]"),
#          paralog = as.numeric(str_extract(gene, "\\d+$"))) %>%
#   filter(!is.na(cluster), !is.na(paralog))
# 
# skingeneorder = skinheatdf %>%
#   distinct(gene, cluster, paralog) %>%
#   arrange(paralog, cluster) %>%
#   pull(gene)
# skinheatdf$gene = factor(skinheatdf$gene, levels = skingeneorder)
# 
# #Making matrix
# skinheatmat = skinheatdf %>%
#   dplyr::select(gene, celltype, mean_delta) %>%
#   pivot_wider(names_from = celltype, values_from = mean_delta) %>%
#   arrange(gene) %>%
#   as.data.frame()
# rownames(skinheatmat) = skinheatmat$gene
# skinheatmat$gene = NULL
# skinheatmat = as.matrix(skinheatmat)
# 
# #Making the color gradient symmetrical around zero for easier interpretation
# skinmax = max(abs(skinheatmat), na.rm = TRUE)
# skinbreaks = seq(-skinmax, skinmax, length.out = 101)
# colors = colorRampPalette(c("#00468B", "white", "#AD000A"))(100)
# 
# #Plot
# pheatmap.skin = pheatmap(skinheatmat,
#                          cluster_rows = FALSE,
#                          cluster_cols = FALSE,
#                          color = colors,
#                          breaks = skinbreaks,
#                          main = "Mean HOX delta (LowBack - Neck) per cell type",
#                          fontsize_row = 6,
#                          fontsize_col = 8,
#                          scale = "none",
#                            angle_col = 315)
# ggsave(file.path(plottingDir, "skin_celltype_collinearity_check.pdf"), 
#        pheatmap.skin, width = 5, height = 7)
# 
# #### MUSCLE
# colmuscle = data.frame()
# heatlist.muscle = list()
# #Will just use profpseudo and profmeta again since they have the structure I want
# for( i in unique(profmeta$celltype)) {
#   
#   keep = profmeta$celltype == i &
#     profmeta$tissue %in% c("Bicep", "Thigh")
#   
#   count = profcounts[, keep]
#   meta = profmeta[keep,]
#   
#   pairtab = table(meta$individual, meta$tissue)
#   
#   if(!all(c("Bicep", "Thigh") %in% colnames(pairtab))) next
#   
#   if(length(unique(meta$individual)) < 2 ||
#      any(pairtab[ , c("Bicep", "Thigh"), drop = FALSE] == 0)) next
#   
#   meta = droplevels(meta)
#   meta$tissue = factor(meta$tissue, levels = c("Bicep", "Thigh"))
#   meta$individual = factor(meta$individual)
#   
#   dds = DESeqDataSetFromMatrix(
#     countData = round(count),
#     colData = meta,
#     design = ~ tissue
#   )
#   vsd = varianceStabilizingTransformation(dds, blind = TRUE) #VST - Variance stabilizing transformation
#   hoxkeep = grepl("^HOX", rownames(assay(vsd)))
#   matrix = assay(vsd)[hoxkeep, ]
#   hoxmat = matrix[grepl("\\d+$", rownames(matrix)), , drop = FALSE] #Keep only HOX genes with numbers at the end of their names
#   
#   if(nrow(hoxmat) < 3) next
#   
#   #Make long
#   long = as.data.frame(as.matrix(hoxmat))
#   long$gene = rownames(long)
#   long = long %>%
#     pivot_longer(cols = -gene, names_to = "sample_id", values_to = "expression")
#   
#   meta2 = meta %>%
#     mutate(sample_id = rownames(meta)) %>%
#     dplyr::select(sample_id, sample, individual, tissue)
#   
#   long = left_join(long, meta2, by = "sample_id")
#   
#   #Calculate anterior-posterior delta
#   wide = long %>%
#     dplyr::select(gene, individual, tissue, expression) %>%
#     pivot_wider(names_from = tissue, values_from = expression) %>%
#     filter(!is.na(Bicep) & !is.na(Thigh)) %>%
#     mutate(delta = Thigh - Bicep,
#            paralog = as.numeric(str_extract(gene, "\\d+$")),
#            celltype = i) %>%
#     filter(!is.na(paralog))
#   
#   #Store info for heatmap
#   heatlist.muscle[[i]] = wide %>%
#     group_by(gene) %>%
#     summarize(mean_delta = mean(delta), .groups = "drop") %>%
#     mutate(celltype = i)
#   
#   #Test collinearity principle
#   fit = lmer(delta ~ paralog + (1|individual), data = wide)
#   print(summary(fit))
#   
#   #Spearman on mean delta per gene
#   mean_delta = wide %>%
#     group_by(gene, paralog) %>%
#     summarize(mean_delta = mean(delta, na.rm = TRUE), .groups = "drop")
#   
#   cor_res = cor.test(mean_delta$paralog, mean_delta$mean_delta, method = "spearman", exact = FALSE)
#   
#   colmuscle = rbind(
#     colmuscle,
#     data.frame(
#       celltype = i,
#       spearman_rho = unname(cor_res$estimate),
#       spearman_pval = cor_res$p.value,
#       lmer_slope = summary(fit)$coefficients["paralog", "Estimate"],
#       lmer_pval = summary(fit)$coefficients["paralog", "Pr(>|t|)"],
#       n_hox_genes = length(unique(wide$gene)),
#       n_individuals = length(unique(wide$individual))
#     ))
#   
# }
# 
# #Save numbers
# colmuscle = colmuscle %>% arrange(desc(lmer_slope))
# write.xlsx(colmuscle, file.path(plottingDir, "COLLINEARITY_HOXGENES_MUSCLE.xlsx"), overwrite = TRUE)
# 
# ### Make heatmap
# #Setting order
# muscleheatdf = bind_rows(heatlist.muscle) %>%
#   mutate(cluster = str_extract(gene, "HOX[ABCD]"),
#          paralog = as.numeric(str_extract(gene, "\\d+$"))) %>%
#   filter(!is.na(cluster), !is.na(paralog))
# 
# musclegeneorder = muscleheatdf %>%
#   
#   distinct(gene, cluster, paralog) %>%
#   arrange(paralog, cluster) %>%
#   pull(gene)
# muscleheatdf$gene = factor(muscleheatdf$gene, levels = musclegeneorder)
# 
# #Making matrix
# muscleheatmat = muscleheatdf %>%
#   dplyr::select(gene, celltype, mean_delta) %>%
#   pivot_wider(names_from = celltype, values_from = mean_delta) %>%
#   arrange(gene) %>%
#   as.data.frame()
# rownames(muscleheatmat) = muscleheatmat$gene
# muscleheatmat$gene = NULL
# muscleheatmat = as.matrix(muscleheatmat)
# 
# #Making the color gradient symmetrical around zero for easier interpretation
# musclemax = max(abs(muscleheatmat), na.rm = TRUE)
# musclebreaks = seq(-musclemax, musclemax, length.out = 101)
# colors = colorRampPalette(c("#00468B", "white", "#AD000A"))(100)
# 
# #Plot
# pheatmap.muscle = pheatmap(muscleheatmat,
#                            cluster_rows = FALSE,
#                            cluster_cols = FALSE,
#                            color = colors,
#                            breaks = musclebreaks,
#                            main = "Mean HOX delta (Thigh - Bicep) per cell type",
#                            fontsize_row = 6,
#                            fontsize_col = 8,
#                            scale = "none",
#                            angle_col = 315)
# ggsave(file.path(plottingDir, "muscle_celltype_collinearity_check.pdf"), 
#        pheatmap.muscle, width = 5, height = 7)
# 
# 
# 
# 
# 
# ###Finished
# message("<><><><><><><><><><> All done <><><><><><><><><><>")
