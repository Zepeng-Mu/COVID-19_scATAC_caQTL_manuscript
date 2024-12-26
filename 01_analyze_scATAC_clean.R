####################################################################################################
## Script purpose: Basic analysis on combined COVID-19, Benaglio, and You clean data
## Author: Zepeng Mu
####################################################################################################
library(tidyverse)
library(data.table)
library(ArchR)
library(batchelor)
library(patchwork)
library(Matrix)
library(matrixStats)
library(sparseMatrixStats)
library(ComplexHeatmap)
library(ArchRWrappers)
library(RColorBrewer)
library(Seurat)
library(ini)
library(pbmcapply)
library(ggtree)
library(BSgenome.Hsapiens.UCSC.hg19)
"%&%" <- function (a, b) paste0(a, b)
source("r_functions.R")
source("scatac/azimuth_markers.R")

addArchRThreads(threads = 3)
addArchRGenome("hg19")

genomeAnnot <- getGenomeAnnotation()

gencode_anno_list <- readRDS("archR_geneAnno_gencodev19New.rds")
gencode_anno <- createGeneAnnotation(
  genome = "hg19",
  TSS = gencode_anno_list$TSS,
  exons = gencode_anno_list$exons,
  genes = gencode_anno_list$genes,
  annoStyle = "ENSEMBL"
)

# Metadata ----
metaDt <- fread("summary_patients_w_admixture_and_captureInfo_edRAG.txt")
metaDtNew <- fread("metadata_simple.csv")
smpsCovid <- metaDt$ID_VIW_PBMC[metaDt$name_sequencing_scATACseq != ""]
smpsBenaglio <- list.dirs("/project/yangili1/zpmu/benaglio2020/data/cellranger_out21", recursive = F, full.names = F)
smpsYou <- list.dirs("/project/yangili1/zpmu/you2021/data/cellranger_out21", recursive = F, full.names = F)

# Load ArchR project ----
projPBMC2 <- loadArchRProject("archR_3study_cmb_amulet_clean_v1")

projPBMC2$Disease <- ifelse(projPBMC2$Sample %in% metaDtNew$sampleID[metaDtNew$infection_status == "case"], "COVID-19",
                            ifelse(projPBMC2$Sample %in% c("M1", "M2", "M4", "M5", "M6", "S1", "S2", "S3"),
                                   "COVID-19 convalescent", "Healthy"))

saveArchRProject(projPBMC2,
                 "archR_3study_cmb_amulet_clean_v1")

projPBMC2 <- loadArchRProject("archR_3study_cmb_amulet_clean_v1")

# Iterative LSI ----
## reducedMNN ----
# See previous script

# Cluster ----
projPBMC2 <- addClusters(
  input = projPBMC2,
  reducedDims = "reducedMNN",
  method = "Seurat",
  name = "ClustersMNN",
  resolution = 0.8,
  corCutOff = 0.75,
  dimsToUse = 1:30,
  scaleDims = T,
  maxClusters = 50,
  force = T
)

# UMAP ----
projPBMC2 <- addUMAP(
  ArchRProj = projPBMC2,
  reducedDims = "reducedMNN",
  name = "UMAPMNN",
  corCutOff = 0.75,
  dimsToUse = 1:30,
  scaleDims = T,
  nNeighbors = 30,
  minDist = 0.6,
  spread = 1.2,
  metric = "cosine",
  threads = 20,
  force = T
)

plotEmbedding(
  ArchRProj = projPBMC2,
  colorBy = "cellColData",
  name = "predictedGroupL1_Co",
  embedding = "UMAPMNN",
  baseSize = 5,
  size = 0.1,
  pal = L1_colors,
  labelMeans = F,
  rastr = T,
  plotAs = "points",
  title = ""
) +
  theme_void() +
  labs(title = "") +
  theme(legend.position = "none")

ggsave(width = 6, height = 6, dpi = 960,
       filename = "/project/yangili1/zpmu/covid19/figs/MS_fig/umapmnn_L1cell.png")

plotEmbedding(
  ArchRProj = projPBMC2,
  colorBy = "cellColData",
  name = "predictedGroupL2_Co",
  embedding = "UMAPMNN",
  baseSize = 5,
  size = 0.1,
  pal = L2_colors,
  labelMeans = F,
  rastr = T,
  plotAs = "points",
  title = ""
) +
  theme_void() +
  labs(title = "") +
  theme(legend.position = "none")

ggsave(width = 6, height = 6, dpi = 960,
       filename = "/project/yangili1/zpmu/covid19/figs/MS_fig/umapmnn_L2cell.png")

g <- plotCellCol(getCellColData(projPBMC2))
g <- g + plot_layout(ncol = 1)

ggsave("/project/yangili1/zpmu/covid19/figs/MS_fig/cellCollData.pdf",
       width = 8, height = 6, plot = g)

# MAGIC imputation ----
projPBMC2 <- addImputeWeights(
  ArchRProj = projPBMC2,
  dimsToUse = 1:30,
  scaleDims = T,
  corCutOff = 0.75,
  reducedDims = "reducedMNN",
  threads = 4
)

saveArchRProject(
  ArchRProj = projPBMC2,
  outputDirectory = "archR_3study_cmb_amulet_clean_v1",
  load = F
)

## L1 markers ----
lapply(c("MS4A1", "CD3E", "CD8A", "NCR1", "S100A8", "FLT3"), function(x) {
  g <- plotEmbedding(
    ArchRProj = projPBMC2,
    colorBy = "GeneScoreMatrix",
    name = x,
    embedding = "UMAPMNN",
    rastr = T,
    baseSize = 5,
    size = 0.1,
    labelMeans = F,
    plotAs = "points"
  ) +
    theme_void() +
    theme(legend.position = "none", plot.title = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "in"))

  ggsave(plot = g, width = 3, height = 3,
         filename = str_glue("/project/yangili1/zpmu/covid19/figs/MS_fig/umapmnn_GAscore_{x}.png"))
})

# Differential genes ----
## L2 cells ----
L2Markers <- getMarkerFeatures(
  ArchRProj = projPBMC2,
  useMatrix = "GeneScoreMatrix",
  groupBy = "predictedGroupL2_Co",
  bias = c("log10(nFrags)", "TSSEnrichment"),
  useSeqnames = c("chr"%&%1:22),
  threads = 3
)

write_rds(L2Markers, "archR_3study_cmb_amulet_clean_v1/cellL2_diffGene.rds")
L2Markers <- read_rds("archR_3study_cmb_amulet_clean_v1/cellL2_diffGene.rds")

plotMarkerHeatmap(L2Markers, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")

heatmapGS <- plotMarkerHeatmap(
  seMarker = L2Markers,
  cutOff = "FDR < 0.005 & Log2FC > 1",
  transpose = F,
  returnMatrix = T
)

rownames(heatmapGS) <- str_remove(rownames(heatmapGS), "^chr[0-9+]:")
myCol <- circlize::colorRamp2(breaks = seq(-2, 2, length.out = 9), colors = ArchRPalettes$blueYellow)

h <- Heatmap(
  t(heatmapGS),
  col = myCol, border = "black",
  border_gp = gpar(lwd = 0.5),
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 3, fontface = "italic"),
  show_column_dend = F, show_row_dend = F,
  column_labels = ifelse(rownames(heatmapGS) %in% unlist(azimuth_pbmc_L2), rownames(heatmapGS), ""),
  show_heatmap_legend = F,
  raster_quality = 10
)

pdf("/project/yangili1/zpmu/covid19/figs/MS_fig/marker_heatmap_L21.pdf", width = 4.3, height = 1.5)
draw(h)
dev.off()

## Sinto input ----
for (smp in unique(projPBMC2$Sample)) {
  getCellColData(projPBMC2, select = c("Sample", "predictedGroupL1_Co"), drop = F) %>%
    as.data.frame() %>%
    filter(Sample == smp) %>%
    mutate(
      bc = str_remove(rownames(.), fixed(smp %&% "#")),
      predictedGroupL1_Co = str_replace_all(predictedGroupL1_Co, " ", "."),
      tmp = Sample %&% "-" %&% predictedGroupL1_Co
    ) %>%
    dplyr::select(bc, tmp) %>%
    fwrite(
      file = str_glue("group_bigwig/L1/L1_{smp}_bc.txt"),
      quote = F, row.names = F, col.names = F, sep = "\t"
    )
}

for (smp in unique(projPBMC2$Sample)) {
  getCellColData(projPBMC2, select = c("Sample", "predictedGroupL2_Co"), drop = F) %>%
    as.data.frame() %>%
    filter(Sample == smp) %>%
    mutate(
      bc = str_remove(rownames(.), fixed(smp %&% "#")),
      predictedGroupL2_Co = str_replace_all(predictedGroupL2_Co, " ", "."),
      tmp = Sample %&% "-" %&% predictedGroupL2_Co
    ) %>%
    dplyr::select(bc, tmp) %>%
    fwrite(
      file = str_glue("group_bigwig/L2/L2_{smp}_bc.txt"),
      quote = F, row.names = F, col.names = F, sep = "\t"
    )
}

# Plot tracks with pgt ----
cells <- str_replace(unique(projPBMC2$predictedGroupL2_Co), " ", ".")

genes1 <- c(
  "LYZ", "MME", "CEBPB", "CD3D", "CD4", "CD8A", "CD8B", "PAX5", "MS4A1", "EBF1",
  "ZNF860", "PLCG2", "NKG7", "NCR1", "IL1B", "CCL3", "IL6R", "IL6", "CCL20")

markersGS <- read_rds("archR_3study_cmb_amulet_clean_v1/cellL2_diffGene.rds")
diffGene <- getMarkers(markersGS, cutOff = "FDR < 0.01 & Log2FC > 1")
tmpGene <- lapply(diffGene, function(x) x$name)
genes <- intersect(unlist(azimuth_pbmc_L2), unlist(tmpGene))

iniFile <- "scatac/fig1AllL2.ini"
for (i in unique(unlist(azimuth_pbmc_L2_small))) {
  if (file.exists(str_glue("L2_{i}.pdf"))) {
    next
  }

  message(i)
  tmp <- gencode_anno$genes[gencode_anno$genes$symbol == i]
  ext <- max(2e3, round(width(tmp) * 0.1))
  start(tmp) <- start(tmp) - round(ext)
  end(tmp) <- end(tmp) + round(ext)

  if (length(tmp) > 1) {
    next
  }

  maxBW <- pbmclapply(cells, function (cell) {
    tmpScore <- rtracklayer::import(
      str_glue("group_bigwig/L2/{cell}.bw"),
      which = tmp
    )
    return(max(tmpScore$score))
  }, mc.cores = 3) %>% unlist() %>% max()

  ini1 <- read.ini(iniFile)
  for (cell in cells) {
    ini1[[cell]]$file <- str_glue("group_bigwig/L2/{cell}.bw")
    ini1[[cell]]$max_value <- ceiling(maxBW * 1)
    ini1[[cell]]$number_of_bins <- width(tmp) %/% 25
  }

  write.ini(ini1, filepath = str_glue("group_bigwig/L2/L2_{i}_{seqnames(tmp)}_{start(tmp)}_{end(tmp)}.ini"))

  system(
    str_glue("grep -w {i} /project/yangili1/zpmu/gencode/v19/gencode.v19.annotation.gtf | grep -v 'AS1' > ",
             "group_bigwig/L2/tmp.gtf")
  )

  cmd <- str_glue(
    "pgt --tracks group_bigwig/L2/",
    "L2_{i}_{seqnames(tmp)}_{start(tmp)}_{end(tmp)}.ini ",
    "--width 5 --region {seqnames(tmp)}:{start(tmp)}-{end(tmp)} ",
    "--outFileName L2_{i}.pdf"
  )

  system(cmd)
}

# Cell type proportion ----
## L1 cell types ----
cellL1 <- c("B", "CD4 T", "CD8 T", "other T", "NK", "Mono", "DC", "other")
id <- factor("D"%&%1:length(unique(projPBMC2$Sample)), levels = "D"%&%1:length(unique(projPBMC2$Sample)))
names(id) <- unique(projPBMC2$Sample)
projPBMC2$ID <- id[projPBMC2$Sample]

propDf <- getCellColData(projPBMC2, select = c("ID", "predictedGroupL1_Co")) %>%
  as.data.frame() %>%
  dplyr::group_by(ID, predictedGroupL1_Co) %>%
  dplyr::tally() %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup() %>%
  pivot_wider(id_cols = ID, names_from = predictedGroupL1_Co, values_from = prop)

propMtrx <- propDf %>% dplyr::select(-ID) %>% as.matrix()
rownames(propMtrx) <- propDf$ID
propMtrx[which(is.na(propMtrx))] <- 0
propClu <- hclust(dist(propMtrx, method = "canberra"), method = "ward.D2")

annoDf <- getCellColData(projPBMC2) %>%
  as.data.frame() %>%
  dplyr::distinct(ID, .keep_all = T) %>%
  dplyr::mutate(node = 1:59) %>%
  dplyr::select(node, ID, Disease)

p <- ggtree(propClu, ladderize = F, lwd = 0.2) %<+% annoDf +
  cowplot::theme_cowplot() +
  layout_dendrogram() +
  geom_tippoint(aes(color = Disease), shape = 16, size = 2, position = position_nudge(x = -0.05)) +
  scale_color_manual(limits = c("Healthy", "COVID-19", "COVID-19 convalescent"), values = c("grey70", "magenta2", "green2")) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "in"))

plotDf <- getCellColData(projPBMC2, select = c("ID", "predictedGroupL1_Co", "Disease")) %>%
  as.data.frame() %>%
  dplyr::group_by(ID, predictedGroupL1_Co) %>%
  dplyr::tally() %>%
  dplyr::ungroup()

xLab <- plotDf %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(lab = str_glue("{ID} ({sum(n)})")) %>%
  dplyr::distinct(ID, .keep_all = T)

g <- ggplot(plotDf, aes(x = factor(ID, levels = propClu$labels[propClu$order]), y = n,
                        fill = factor(predictedGroupL1_Co, levels = cellL1))) +
  cowplot::theme_cowplot(font_size = 10, font_family = "Helvetica", line_size = 0.4) +
  geom_col(color = NA, position = "fill", width = 0.7, size = 0.3) +
  scale_fill_manual(
    limits = cellL1,
    values = L1_colors[cellL1],
    name = ""
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(labels = xLab$lab[match(propClu$labels[propClu$order], xLab$ID)]) +
  labs(x = "59 Samples", y = "L1 cell type proportions") +
  theme(axis.text.x = element_blank())

pOut <- p / g + plot_layout(heights = c(0.7, 2))
ggsave("/project/yangili1/zpmu/covid19/figs/MS_fig/all_barplot_3study_fixed.pdf",
       width = 8, height = 2, plot = pOut & theme(legend.position = "none"))

lgd <- cowplot::get_legend(g)
ggsave("/project/yangili1/zpmu/covid19/figs/MS_fig/all_barplot_3study_lgd_fixed.pdf",
       width = 2, height = 3, plot = lgd)

lgd1 <- cowplot::get_legend(p)
ggsave("/project/yangili1/zpmu/covid19/figs/MS_fig/all_barplot_3study_tree_lgd_fixed.pdf",
       width = 2, height = 1, plot = lgd1)

## L2 cell types ----
plotDf <- getCellColData(projPBMC2, select = c("Sample", "predictedGroupL2_Co", "Disease")) %>%
  as.data.frame() %>%
  dplyr::group_by(Sample, predictedGroupL2_Co) %>%
  dplyr::tally() %>%
  dplyr::ungroup()

g <- ggplot(plotDf, aes(x = Sample, y = n, fill = predictedGroupL2_Co)) +
  geom_col(color = NA, position = "fill", width = 0.75, size = 0.3) +
  scale_fill_manual(
    limits = unique(plotDf$predictedGroupL2_Co),
    values = L2_colors[unique(plotDf$predictedGroupL2_Co)],
    name = ""
  )

lgd <- cowplot::get_legend(g)
ggsave("/project/yangili1/zpmu/covid19/figs/MS_fig/cellL2_lgd.pdf",
       width = 5, height = 5, plot = lgd)

mclust::adjustedRandIndex(projPBMC2$ClustersMNN, projPBMC2$predictedGroupL1_Co)
mclust::adjustedRandIndex(projPBMC2$ClustersMNN, projPBMC2$predictedGroupL2_Co)

cM <- as.matrix(confusionMatrix(projPBMC2$ClustersMNN, projPBMC2$predictedGroupL2_Co))
cM <- cM / rowSums(cM)

pdf("/project/yangili1/zpmu/covid19/figs/MS_fig/azimuth_confusion.pdf", width = 6, height = 6)
Heatmap(as.matrix(cM), name = "Row\nproportion", border = T,
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D",
        col = colorRampPalette(brewer.pal(9, "YlGnBu"))(100))
dev.off()

for (tmpL1 in names(cellL2)) {
  message(tmpL1)
  propDf <- getCellColData(projPBMC2, select = c("ID", "predictedGroupL1_Co", "predictedGroupL2_Co")) %>%
    as.data.frame() %>%
    dplyr::filter(predictedGroupL1_Co == tmpL1) %>%
    dplyr::group_by(ID, predictedGroupL2_Co) %>%
    dplyr::tally() %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup() %>%
    pivot_wider(id_cols = ID, names_from = predictedGroupL2_Co, values_from = prop)

  propMtrx <- propDf %>% dplyr::select(-ID) %>% as.matrix()
  rownames(propMtrx) <- propDf$ID
  propMtrx[which(is.na(propMtrx))] <- 0
  propClu <- hclust(dist(propMtrx, method = "canberra"), method = "ward.D2")

  annoDf <- getCellColData(projPBMC2) %>%
    as.data.frame() %>%
    dplyr::distinct(ID, .keep_all = T) %>%
    dplyr::mutate(node = 1:59) %>%
    dplyr::select(node, ID, Disease)

  p <- ggtree(propClu, ladderize = F, lwd = 0.2) %<+% annoDf +
    cowplot::theme_cowplot() +
    layout_dendrogram() +
    geom_tippoint(aes(color = Disease), shape = 16, size = 2, position = position_nudge(x = -0.05)) +
    scale_color_manual(limits = c("Healthy", "COVID-19", "COVID-19 convalescent"), values = c("grey70", "magenta2", "green2")) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "in"), legend.position = "none")

  plotDf <- getCellColData(projPBMC2, select = c("ID", "predictedGroupL1_Co", "predictedGroupL2_Co", "Disease")) %>%
    as.data.frame() %>%
    dplyr::filter(predictedGroupL1_Co == tmpL1) %>%
    dplyr::group_by(ID, predictedGroupL2_Co) %>%
    dplyr::tally() %>%
    dplyr::ungroup()

  xLab <- plotDf %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(lab = str_glue("{ID} ({sum(n)})")) %>%
    dplyr::distinct(ID, .keep_all = T)

  g <- ggplot(plotDf, aes(x = factor(ID, levels = propClu$labels[propClu$order]), y = n,
                          fill = factor(predictedGroupL2_Co, levels = cellL2[[tmpL1]]))) +
    cowplot::theme_cowplot(font_size = 10, font_family = "Helvetica", line_size = 0.4) +
    geom_col(color = NA, position = "fill", width = 0.7, size = 0.3) +
    scale_fill_manual(limits = levels(cellL2[[tmpL1]]), values = L2_colors[levels(cellL2[[tmpL1]])],
                      name = "") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(labels = xLab$lab[match(propClu$labels[propClu$order], xLab$ID)]) +
    labs(x = "59 Samples", y = str_glue("L2 cell type proportions: {tmpL1}")) +
    theme(axis.text.x = element_blank())

  pOut <- p / g + plot_layout(heights = c(0.7, 2))
  ggsave(str_glue("/project/yangili1/zpmu/covid19/figs/MS_fig/all_barplot_3study_{tmpL1}_fixed.pdf"),
         width = 12, height = 3, plot = pOut)
}

# Markers ----
markersGS <- getMarkerFeatures(
  ArchRProj = projPBMC2,
  useMatrix = "GeneScoreMatrix",
  groupBy = "ClustersMNN",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  threads = 5
)

write_rds(markersGS, "archR_3study_cmb_amulet_clean_v1/markerGS_clu.rds")

markersGS1 <- getMarkerFeatures(
  ArchRProj = projPBMC2,
  useMatrix = "GeneScoreMatrix",
  groupBy = "predictedGroupL2_Co",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  threads = 6
)

write_rds(markersGS1, "archR_3study_cmb_amulet_clean_v1/markerGS_L2.rds")

markersGS <- read_rds("archR_3study_cmb_amulet_clean_v1/markerGS_clu.rds")
diffGene <- getMarkers(markersGS, cutOff = "FDR < 0.01 & Log2FC > 0.5")

# Group coverages ----
projPBMC2 <- addGroupCoverages(ArchRProj = projPBMC2, groupBy = "ClustersMNN", sampleLabels = "Study",
                               excludeChr = c("chrY", "chrM"), threads = 10, force = T)

projPBMC2 <- addReproduciblePeakSet(
  ArchRProj = projPBMC2,
  reproducibility = "2",
  excludeChr = c("chrY", "chrM"),
  pathToMacs2 = "/scratch/midway3/zepengmu/miniconda3/envs/R42/bin/macs3",
  groupBy = "ClustersMNN",
  minCells = 100,
  threads = 8,
  force = T
)

peakSet <- getPeakSet(projPBMC2)
peakSet$name <- as.character(str_glue("{seqnames(peakSet)}_{start(peakSet)}-{end(peakSet)}"))
names(peakSet) <- peakSet$name
peakSet <- ArchR:::.fastAnnoPeaks(peaks = peakSet, BSgenome = "hg19", geneAnnotation = gencode_anno, promoterRegion = c(3000, 500))
projPBMC2 <- addPeakSet(projPBMC2, peakSet, force = T)
projPBMC2 <- addPeakMatrix(projPBMC2, ceiling = Inf, threads = 10, force = T)

saveArchRProject(
  ArchRProj = projPBMC2,
  outputDirectory = "archR_3study_cmb_amulet_clean_v1",
  load = F
)

# Motif annotation ----
projPBMC2 <- addMotifAnnotations(
  ArchRProj = projPBMC2,
  motifSet = "cisbp",
  annoName = "cisbp",
  force = T
)

projPBMC2 <- addBgdPeaks(projPBMC2)

projPBMC2 <- addDeviationsMatrix(
  ArchRProj = projPBMC2,
  peakAnnotation = "cisbp",
  matrixName = "cisbpMatrix",
  threads = 6,
  force = T
)

saveArchRProject(
  ArchRProj = projPBMC2,
  outputDirectory = "archR_3study_cmb_amulet_clean_v1",
  load = F,
  threads = 3
)

# Co-accessible peaks ----
peakSet <- getPeakSet(projPBMC2)
names(peakSet) <- peakSet$name

projPBMC2 <- addCoAccessibility(
  ArchRProj = projPBMC2,
  reducedDims = "reducedMNN",
  dimsToUse = 1:30,
  maxDist = 5e5,
  scaleDims = T,
  corCutOff = 0.75,
  threads = 5
)

saveArchRProject(
  ArchRProj = projPBMC2,
  outputDirectory = "archR_3study_cmb_amulet_clean_v1",
  load = F,
  threads = 3
)

cA <- getCoAccessibility(
  ArchRProj = projPBMC2,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = F
)

write_rds(cA, file = "archR_3study_cmb_amulet_clean_v1/coAcc.rds")
