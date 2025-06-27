library(ggplot2)
library(gprofiler2)
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(dplyr)
library(glmGamPoi)
library(qs)

fap.markers <- c("Pdgfra", "Ly6a", "Col1a1", "Col3a1", "Pi16", "Postn") 
tenocyte.markers <- c("Fmod", "Tnmd", "Scx")
endothelial.markers <- c("Pecam1", "Kdr", "Egfl7", "Cd36", "Cd34")
satellite.markers <- c("Pax7", "Ncam1", "Myod1") 
pericyte.markers <- c("Acta2", "Abcc9", "Kcnj8", "Rgs5", "Pdgfrb")
momp.markers <- c("Adgre1", "Fcgr1", "Cd68", "Ly6c2", "Chil3", "Lyz2")
mo.markers <- c("Plac8", "Ifitm6", "Ly6c2", "Chil3", "Ifitm3")
inflammatorymp.markers <- c("Arg1", "Cxcl1", "Thbs1", "Ccl6", "Ccl9") 
residentmp.markers <- c("Pltp", "Selenop", "C1qc", "Ms4a7", "C1qa", "H2-Aa")
fibroticmp.markers <- c("Gpnmb", "Fabp5", "Spp1", "Syngr1", "Ctsd")
ifnmp.markers <- c("Cxcl10", "Isg15", "Rsad2", "Ifit3", "Ms4a4c")
proliferativemp.markers <- c("Stmn1", "Tuba1b", "Pclaf", "Tubb5", "Top2a")
dendritic.markers <- c("Itgax", "H2-Aa", "H2-Eb1", "Cd74", "Cd209a", "Xcr1")
neutrophil.markers <- c("S100a8", "S100a9", "Retnlg", "Cxcr2", "Cxcl2") 
t.markers <- c("Cd3d", "Cd3e", "Cd3g", "Trdc", "Tcrg-C1") 
b.markers <- c("Iglc1", "Ighm", "Igkc", "Cd79a", "Ly6d") 
nk.markers <- c("Nkg7", "Ccl5", "Trbc2", "Cd3g")
nkt.markers <- c("Gata3", "Nkg7", "Trbc2", "Tcrg-C1", "Cd69") 
schwann.markers <- c("Sox10", "S100b") 
mesothelial.markers <- c("Upk1b", "Upk3b", "Lrrn4")

satellite.quiescence <- c("Tgfb1", "Jag1", "Jag2")
satellite.activation <- c("Hgf", "Fgf2", "Igf1", "Pdgfb", "Il1a", "Il1b", "Tnfsf12", "Il6", "Cxcl10", "Cxcl12", "Adamts1", "Nampt", "Glul", "Fn1", "Col6a1", "Glud1", "Tgfb1", "Cxcl14")
satellite.differentiation <- c("Igf1", "Il10", "Cxcl10")
myocyte.fusion <- c("Igf1", "Gdf3", "Il1a", "Il1b", "Spp1")

M1 <- c("Cd68", "Marco", "Tnf", "Il1b", "Il6", "Il12a", "Il23a", "Il12b", "Nos2", "Nfkbiz", "Stat1", "Socs1", "Irf5")
M2 <- c("Mrc1", "Retnla", "Cd163", "Il4ra", "Il10", "Tgfb1", "Igf1", "Arg1", "Pparg", "Socs3", "Stat6", "Irf4")
Additional <- c("Il1a", "Ccl2", "Ccl7", "Cxcl2", "Pdgfa", "Fn1", "Ccn2", "Spp1", "Mmp2", "Mmp9", "Mmp12", "Mmp14", "Timp1", "Timp2", "Nampt", "Glul", "Glud1", "Adamts1", "Cxcl10", "Tnfsf12", "Ly6c2", "Lyve1", "Timd4", "Csf1r", "Adgre1", "Fcgr1")

individual <- c("Tnf", "Il1b", "Ifna1", "Ifna2", "Ifna4", "Ifna5", "Ifna6", "Ifna7", "Ifna8", "Ifna10", "Ifna13", "Ifna14", "Ifna16", "Ifna17", "Ifna21", "Ifnb1", "Ifng", "Ltbp4", "Mmp14")
updateMetadata <- function(seuratObject) {
  x <- seuratObject
  x$log10GenesPerUMI <- log10(x$nFeature_RNA) / log10(x$nCount_RNA)
  x$mitoRatio <- PercentageFeatureSet(object = x, pattern = "^mt-")
  x$mitoRatio <- x@meta.data$mitoRatio / 100
  metadata <- x@meta.data
  
  # Add cell IDs to metadata
  metadata$cells <- rownames(metadata)
  
  # Rename columns
  metadata <- metadata %>%
    dplyr::rename(nUMI = nCount_RNA,
                  nGene = nFeature_RNA)
  
  x@meta.data <- metadata
  
  return(x)
}

plotQC <- function(seuratObject) {
  metadata <- seuratObject@meta.data
  
  vln <- VlnPlot(seuratObject, features = c("nGene", "nUMI", "mitoRatio"), ncol = 3, pt.size = 0.000001)
  print(vln)
  
  umi.density <- metadata %>% 
    ggplot(aes(x=nUMI)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density")
  
  print(umi.density)
  
  # Visualize the distribution of genes detected per cell via histogram
  gene.density <- metadata %>% 
    ggplot(aes(x=nGene)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10()
  
  print(gene.density)
  
  # Visualize the distribution of genes detected per cell via boxplot
  cells.v.genes <- metadata %>% 
    ggplot(aes(y=log10(nGene))) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells vs NGenes")
  
  print(cells.v.genes)
  
  quality <- metadata %>% 
    ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250)
  
  print(quality)
  
  
  # Visualize the distribution of mitochondrial gene expression detected per cell
  mito.density <- metadata %>% 
    ggplot(aes(x=mitoRatio)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic()
  
  print(mito.density)
}

remove_doublets <- function(seuratObject) {
  sce <- as.SingleCellExperiment(seuratObject, assay = "RNA")
  sce <- scDblFinder(sce, samples = seuratObject@meta.data$orig.ident)
  
  seuratObject@meta.data$doublet.class <- sce$scDblFinder.class
  seuratObject <- subset(seuratObject, subset=doublet.class == "singlet")
  return(seuratObject)
}

check_cycle <- function(x, returnObject = TRUE) {
  # Normalize the counts
  seurat_phase <- NormalizeData(x)
  seurat_phase <- FindVariableFeatures(seurat_phase)
  seurat_phase <- ScaleData(seurat_phase, features = rownames(seurat_phase))
  
  s.genes <- gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
  g2m.genes <- gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
  
  # Score cells for cell cycle
  seurat_phase <- CellCycleScoring(seurat_phase, 
                                   g2m.features = g2m.genes, 
                                   s.features = s.genes)
  seurat_phase$CC.Difference <- seurat_phase$S.Score - seurat_phase$G2M.Score
  if (returnObject){
    return(seurat_phase)
  }
  else {
    x@meta.data$Phase <- seurat_phase@meta.data$Phase
    x@meta.data$CC.Difference <- seurat_phase@meta.data$CC.Difference
    return(x)
  }
}

dim_cluster <- function(seuratObject, res=0.5, cluster.name=NULL) {
  features <- VariableFeatures(seuratObject)
  seuratObject <- RunPCA(seuratObject, features = features)
  seuratObject <- FindNeighbors(seuratObject, dims = 1:10)
  seuratObject <- FindClusters(seuratObject, resolution = res, cluster.name = cluster.name)
  seuratObject <- RunUMAP(seuratObject, dims = 1:10)
  return(seuratObject)
}

output_markers <- function(seuratObject, cluster_on = "seurat_clusters", n_markers = 30, prefix) {
  Idents(object = seuratObject) <- cluster_on
  markers <-FindAllMarkers(seuratObject, 
                           only.pos = FALSE,
                           min.pct =  0.25, 
                           min.diff.pct = 0.1)
  
  markers <- markers %>% arrange(cluster, desc(avg_log2FC))
  write.csv(markers, paste0(prefix, "_all", ".csv"))
  
  pos.markers <- filter(markers, avg_log2FC > 0)
  marker.name <- paste(prefix, "markers", n_markers, sep="_")
  pos.markers <- pos.markers %>%
    group_by(cluster) %>%
    slice_max(n = n_markers, order_by = avg_log2FC)
  write.csv(pos.markers, paste0(marker.name, ".csv"))
  return(markers)
}

output_FC <- function(seuratObject, cluster_on, ident, prefix) {
  Idents(object = seuratObject) <- cluster_on
  genes <- FoldChange(seuratObject, ident.1=ident)
  genes$symbol <- row.names(genes)
  genes <- genes %>% arrange(desc(avg_log2FC))
  genes <- dplyr::select(genes, symbol, avg_log2FC)
  write.table(genes, paste0(prefix, ".rnk"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
}

plotMarkers <- function(seuratObject, cluster_on = "seurat_clusters") {
  Idents(object = seuratObject) <- cluster_on
  fap <- DotPlot(seuratObject, features = fap.markers) + ggtitle("FAP Markers")
  print(fap)
  tenocyte <- DotPlot(seuratObject, features = tenocyte.markers) + ggtitle("Tenocyte Markers")
  print(tenocyte)
  endothelial <- DotPlot(seuratObject, features = endothelial.markers) + ggtitle("Endothelial Markers")
  print(endothelial)
  satellite <- DotPlot(seuratObject, features = satellite.markers) + ggtitle("Satellite Markers")
  print(satellite)
  pericyte <- DotPlot(seuratObject, features = pericyte.markers) + ggtitle("Pericyte Markers")
  print(pericyte)
  momp <- DotPlot(seuratObject, features = momp.markers) + ggtitle("MO/MP Markers")
  print(momp)
  dendritic <- DotPlot(seuratObject, features = dendritic.markers) + ggtitle("Dendritic Markers")
  print(dendritic)
  neutrophil <- DotPlot(seuratObject, features = neutrophil.markers) + ggtitle("Neutrophil Markers")
  print(neutrophil)
  t <- DotPlot(seuratObject, features = t.markers) + ggtitle("T Markers")
  print(t)
  b <- DotPlot(seuratObject, features = b.markers) + ggtitle("B Markers")
  print(b)
  nk <- DotPlot(seuratObject, features = nk.markers) + ggtitle("NK Markers")
  print(nk)
  nkt <- DotPlot(seuratObject, features = nkt.markers) + ggtitle("NKT Markers")
  print(nkt)
  schwann <- DotPlot(seuratObject, features = schwann.markers) + ggtitle("Schwann Markers")
  print(schwann)
  mesothelial <- DotPlot(seuratObject, features = mesothelial.markers) + ggtitle("Mesothelial Markers")
  print(mesothelial)
}

plotMompMarkers <- function(seuratObject, cluster_on = "seurat_clusters", prefix = '', save = FALSE) {
  Idents(object = seuratObject) <- cluster_on
  mo <- DotPlot(seuratObject, features = mo.markers) + ggtitle("MO Markers")
  print(mo)
  inflammatory <- DotPlot(seuratObject, features = inflammatorymp.markers) + ggtitle("Pro-Inflammatory MP Markers")
  print(inflammatory)
  resident <- DotPlot(seuratObject, features = residentmp.markers) + ggtitle("Resident-Like MP Markers")
  print(resident)
  fibrotic <- DotPlot(seuratObject, features = fibroticmp.markers) + ggtitle("Pro-Fibrotic MP Markers")
  print(fibrotic)
  ifn <- DotPlot(seuratObject, features = ifnmp.markers) + ggtitle("IFN-Activated MP Markers")
  print(ifn)
  proliferative <- DotPlot(seuratObject, features = proliferativemp.markers) + ggtitle("Proliferative MP Markers")
  print(proliferative)
  
  if (save) {
    ggsave(paste0(prefix, "MO.png"), mo, width = 7, height = 4.32, units = "in")
    ggsave(paste0(prefix, "Inflammatory.png"), inflammatory, width = 7, height = 4.32, units = "in")
    ggsave(paste0(prefix, "ResLike.png"), resident, width = 7, height = 4.32, units = "in")
    ggsave(paste0(prefix, "Fibrotic.png"), fibrotic, width = 7, height = 4.32, units = "in")
    ggsave(paste0(prefix, "IFN.png"), ifn, width = 7, height = 4.32, units = "in")
    ggsave(paste0(prefix, "Proliferating.png"), proliferative, width = 7, height = 4.32, units = "in")
  }
}