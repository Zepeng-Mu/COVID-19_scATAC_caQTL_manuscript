####################################################################################################
## Script purpose: Topic model analysis on all cells from 3 datasets
## Author: Zepeng Mu
####################################################################################################
library(tidyverse)
library(data.table)
library(ArchR)
library(fastTopics)
library(Matrix)
library(matrixStats)
library(sparseMatrixStats)
library(cowplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(lme4)
library(RColorBrewer)
library(ggrepel)
library(colorspace)
library(pbmcapply)
library(ini)
"%&%" <- function(a, b) paste0(a, b)
source("r_functions.R")
source("azimuth_markers.R")
addArchRThreads(4)
setDTthreads(4)

getKnnGroup <- function(ArchRProj, reducedDims, knnIteration = 500, overlapCutoff = 0.8, k = 100, ...) {
  rD <- getReducedDims(ArchRProj = ArchRProj, reducedDims = reducedDims, ...)
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)
  knnObj <- ArchR:::.computeKNN(data = rD, query = rD[idx, ], k = k)

  #Determin Overlap
  keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn == 0, ]

  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList

  return(knnObj)
}

getKnnGroupMtrx <- function(mtrx, knnIteration = 500, overlapCutoff = 0.8, k = 100, ...) {
  idx <- sample(seq_len(nrow(mtrx)), knnIteration, replace = !nrow(mtrx) >= knnIteration)
  knnObj <- ArchR:::.computeKNN(data = mtrx, query = mtrx[idx, ], k = k)

  #Determin Overlap
  keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn == 0, ]

  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(mtrx)[knnObj[x, ]]
  }) %>% SimpleList

  return(knnObj)
}

# Load ArchR project ----
gencode_anno_list <- readRDS("archR_geneAnno_gencodev19New.rds")
gencode_anno <- createGeneAnnotation(
  genome = "hg19",
  TSS = gencode_anno_list$TSS,
  exons = gencode_anno_list$exons,
  genes = gencode_anno_list$genes,
  annoStyle = "ENSEMBL"
)

projPBMC <- loadArchRProject("archR_3study_cmb_amulet_clean_v1")
peakSet <- getPeakSet(projPBMC)

cellGrp <- getKnnGroup(projPBMC, reducedDims = "reducedMNN", scale = T, knnIteration = 4000, k = 50)

## Gene to peak link ----
cA <- read_rds("archR_3study_cmb_amulet_clean_v1/coAcc.rds")
cAMeta <- metadata(cA)$peakSet
names(cAMeta) <- str_glue("{seqnames(cAMeta)}_{start(cAMeta)}-{end(cAMeta)}")
cAMeta <- peakSet[names(cAMeta)]
cADf <- as.data.frame(cA)
cADf <- cADf %>%
  mutate(queryName = names(cAMeta)[queryHits],
         queryType = cAMeta$peakType[queryHits],
         queryNearestGene = cAMeta$nearestGene[queryHits],
         queryNearestTSS = cAMeta$nearestTSS[queryHits],
         subjectName = names(cAMeta)[subjectHits],
         subjectType = cAMeta$peakType[subjectHits],
         subjectNearestGene = cAMeta$nearestGene[subjectHits],
         subjectNearestTSS = cAMeta$nearestTSS[subjectHits],
         dist = distance(resize(cAMeta[queryHits], 1), resize(cAMeta[subjectHits], 1)),
         linkType = queryType%&%"-"%&%subjectType) %>%
  mutate(linkType = case_when(
    linkType == "Exonic-Distal" ~ "Distal-Exonic",
    linkType == "Intronic-Distal" ~ "Distal-Intronic",
    linkType == "Promoter-Distal" ~ "Distal-Promoter",
    linkType == "Intronic-Exonic" ~ "Exonic-Intronic",
    linkType == "Promoter-Exonic" ~ "Exonic-Promoter",
    linkType == "Promoter-Intronic" ~ "Intronic-Promoter",
    T ~ linkType
  ))

cAPromoter <- cADf[cADf$queryType == "Promoter", ]

# Get peak matrix and peak set ----
peakCntMtrx <- getMatrixFromProject(
  ArchRProj = projPBMC,
  useMatrix = "PeakMatrix",
  useSeqnames = c("chr"%&%1:22, "chrX"),
  threads = 3
)

peakCntValue <- assay(peakCntMtrx)
rownames(peakCntValue) <- rowData(peakCntMtrx)$name
peakCntValue <- t(peakCntValue)

metaDt <- colData(peakCntMtrx)
rrPeak <- rowRanges(peakCntMtrx)
names(rrPeak) <- rrPeak$name
rm(peakCntMtrx)
gc()

# Model fitting ----
fit <- read_rds("archR_3study_cmb_amulet_clean_v1/fastTopics/peak_fastTopics_fitCond_K20.rds")

peakCntValue <- peakCntValue[rownames(fit$L), ]
peakCntValue <- peakCntValue[, rownames(fit$F)]
fitF <- getF(fit, peakCntValue)

cellGrpL <- getKnnGroupMtrx(scale(fit$L), scale = T, knnIteration = 4000, k = 50)

metaDt <- getCellColData(projPBMC)
metaDt$cellNames <- rownames(metaDt)
metaDt <- metaDt[rownames(fit$L), ]

twentyCol <- c("pink1", "violet", "magenta", "mediumpurple1", "purple", "turquoise2", "skyblue",
               "steelblue", "blue2", "orange", "palevioletred", "red3", "springgreen2",
               "yellowgreen", "palegreen4", "wheat2", "tan", "tan3", "brown", "grey70")

names(twentyCol) <- "k"%&%1:20

## Plot UMAP ----
UMAPMNN <- getEmbedding(projPBMC, "UMAPMNN")
UMAPMNN <- UMAPMNN[rownames(fit$L), ]
p <- pbmclapply("k"%&%1:20, function(k) {
  plotTopicsEmbedding(UMAPMNN, fit, k, normalize = T, colorTitle = k, randomize = T) +
    theme_cowplot(font_size = 5, font_family = "Helvetica", line_size = 0.4) +
    guides(color = "none", fill = "none") +
    ggtitle(label = "", subtitle = k) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "in"),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0, vjust = 1))
}, mc.cores = 2)

pdf("all_umap_topic20.pdf", width = 4, height = 6)
cowplot::plot_grid(plotlist = p, ncol = 4)
dev.off()

# K means clustering for figures ----
fit <- readRDS("archR_3study_cmb_amulet_clean_v1/fastTopics/peak_fastTopics_fitCond_K20.rds")
Lscale <- scale(fit$L, center = T, scale = F)
distL <- dist(Lscale, method = "euclidian")
distL <- distances::distances(Lscale)

## Figure 2 ----
my_structure_call <- function(dat, colors, ticks = NULL, font.size = 6) {
  ggplot(dat, aes_string(x = "sample", y = "prop", color = "topic", fill = "topic")) +
    theme_cowplot(font_size = 10, font_family = "Helvetica", line_size = 0.4) +
    ggrastr::rasterise(geom_col(), dpi = 1600) +
    scale_x_continuous(
      limits = c(0, max(dat$sample) + 1), breaks = ticks, labels = names(ticks),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 1), labels = c("0.0", "1.0")) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = "", y = "") +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = ggtext::element_markdown(vjust = c(0, 1)),
          legend.position = "none", panel.border = element_rect(colour = "black"))
}

pca <- prcomp(fit$L)$x
tmpClu <- kmeans(pca, centers = 30, iter.max = 100)$cluster
cell2plot <- lapply(unique(tmpClu), function(x) {
  tmp <- tmpClu[tmpClu == x]
  tmpCell <- sample(names(tmp), length(tmp) / 50)
  return(tmpCell)
}) %>% unlist() %>% unique()

metaDt <- getCellColData(projPBMC)
metaDt$cellNames <- rownames(metaDt)
metaDt <- metaDt[rownames(fit$L), ]
metaDt2plot <- metaDt[cell2plot, ]
cellOrder <- fastTopics::structure_plot_default_embed_method(select_loadings(fit, cell2plot))

g <- structure_plot(
  select_loadings(fit, cell2plot),
  colors = twentyCol,
  topics = 1:20,
  loadings_order = order(cellOrder),
  gap = 70,
  grouping = factor(metaDt2plot$predictedGroupL1_Co,
                    levels = c("B", "CD4 T", "CD8 T", "other T", "NK", "Mono", "DC", "other")),
  ggplot_call = my_structure_call
)

ggsave(filename = "topics20.pdf", height = 2, width = 6, plot = g, dpi = 1960)

for (tmpK in seq(6, 18, 2)) {
  message(tmpK)
  tmpFit <- readRDS(str_glue("archR_3study_cmb_amulet_clean_v1/fastTopics/",
                             "peak_fastTopics_fitCond_K{tmpK}.rds"))
  g <- structure_plot(
    select_loadings(tmpFit, cell2plot),
    colors = twentyCol,
    topics = 1:tmpK,
    loadings_order = order(cellOrder),
    gap = 70,
    grouping = factor(metaDt2plot$predictedGroupL1_Co,
                      levels = c("B", "CD4 T", "CD8 T", "other T", "NK", "Mono", "DC", "other")),
    ggplot_call = my_structure_call
  )

  ggsave(filename = str_glue("topics{tmpK}.pdf"), height = 2, width = 6, plot = g, dpi = 1600)
}

### Structure plot legend ----
g <- ggplot(mapping = aes_(x = factor(1:20), y = 1, fill = twentyCol)) +
  theme_void(base_size = 4) +
  geom_tile(color = "grey30") +
  scale_x_discrete(labels = "k"%&%1:20, position = "bottom") +
  scale_fill_identity() +
  theme(legend.position = "none", axis.text.x = element_text()) +
  coord_fixed()

ggsave(filename = "topicsLegend.pdf", height = 0.5, width = 2.5, plot = g)

### Heatmap by K ----
metaDt <- metaDt[rownames(fit$L), ]

ht <- plotKbyGroup(
  fit, metaDt$predictedGroupL2_Co, norm = T,
  column_title = NA,
  show_row_dend = F, show_column_dend = F,
  width = unit(0.06 * 20, "in"),
  height = unit(0.06 * 21, "in"),
  column_title_gp = gpar(fontsize = 5),
  row_title_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  row_names_gp = gpar(fontsize = 5),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  heatmap_legend_param = list(
    title = "Mean loadings", legend_direction = "horizontal", border = NULL,
    legend_width = unit(0.6, "in"), labels_gp = gpar(fontsize = 5),
    title_gp = gpar(fontsize = 5), grid_height = unit(0.07, "in")
  ),
  border_gp = gpar(lwd = 0.5),
  use_raster = T,
  raster_quality = 10
)

pdf("topic20_by_cellL2.pdf")
draw(ht, heatmap_legend_side = "top")
dev.off()

ht <- plotKbyGroup(
  fit, metaDt$Sample, norm = F,
  row_title = "Samples",
  column_title = NA,
  show_row_dend = F, show_column_dend = F,
  width = unit(0.06 * 20, "in"),
  height = unit(0.06 * length(unique(metaDt$Sample)), "in"),
  column_title_gp = gpar(fontsize = 5),
  row_title_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  row_names_gp = gpar(fontsize = 5),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  heatmap_legend_param = list(
    title = "Mean Loadings", legend_direction = "horizontal", border = NULL,
    legend_width = unit(0.6, "in"), labels_gp = gpar(fontsize = 5),
    title_gp = gpar(fontsize = 5), grid_height = unit(0.07, "in")
  ),
  border_gp = gpar(lwd = 0.5),
  use_raster = T,
  raster_quality = 10
)

pdf("topic20_by_sample.pdf")
draw(ht, heatmap_legend_side = "top")
dev.off()

ht <- plotKbyGroup(
  fit, metaDt$Disease,
  row_title = "Disease",
  column_title = NA,
  show_row_dend = F, show_column_dend = F,
  width = unit(0.06 * 20, "in"),
  height = unit(0.06 * 3, "in"),
  column_title_gp = gpar(fontsize = 5),
  row_title_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  row_names_gp = gpar(fontsize = 5),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  heatmap_legend_param = list(
    title = "Mean Loadings", legend_direction = "horizontal", border = NULL,
    legend_width = unit(0.6, "in"), labels_gp = gpar(fontsize = 5),
    title_gp = gpar(fontsize = 5), grid_height = unit(0.07, "in")
  ),
  border_gp = gpar(lwd = 0.5),
  use_raster = T,
  raster_quality = 10
)

pdf("topic20_by_disease_fixed.pdf")
draw(ht, heatmap_legend_side = "top")
dev.off()

p <- lapply("k"%&% 1:20, function(k) {
  g <- ggplot() +
    theme_cowplot(font_size = 10, font_family = "Helvetica", line_size = 0.4) +
    geom_hline(yintercept = c(0, 1), lty = 2, color = "grey", lwd = 0.4) +
    stat_ecdf(aes(fit$L[metaDt$Disease == "Healthy", k], color = "Healthy"), linewidth = 0.4) +
    stat_ecdf(aes(fit$L[metaDt$Disease == "COVID-19 convalescent", k], color = "COVID-19 convalescent"), linewidth = 0.4) +
    stat_ecdf(aes(fit$L[metaDt$Disease == "COVID-19", k], color = "COVID-19"), linewidth = 0.4) +
    scale_color_manual(limits = c("Healthy", "COVID-19 convalescent", "COVID-19"), values = c("grey70", "green2", "magenta2")) +
    scale_x_continuous(breaks = c(0, 1), labels = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 1), labels = c(0, 1)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ylab(k) +
    theme(aspect.ratio = 1, legend.position = "none", plot.margin = unit(c(0.1, 0.1, 0, 0), "in"),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0, vjust = 1))

  return(g)
})

pdf("ecdf_disease_topic20_fixed.pdf", width = 6, height = 6)
cowplot::plot_grid(plotlist = p, ncol = 4)
dev.off()

# B cell topics ----
structure_plot(
  select_loadings(fit, which(metaDt$predictedGroupL1_Co == "B")),
  colors = twentyCol,
  topics = 1:20,
  gap = 40,
  grouping = factor(metaDt$predictedGroupL2_Co[metaDt$predictedGroupL1_Co == "B"],
                    levels = c("B naive", "B intermediate", "B memory", "Plasmablast"))
)

fitLB <- fit$L[metaDt$cellNames[metaDt$predictedGroupL1_Co == "B"], ]

ggplot() +
  theme_cowplot(font_size = 10, line_size = 0.4, font_family = "Helvetica") +
  geom_hline(yintercept = c(0, 1), lty = 2, col = "grey60", lwd = 0.3) +
  stat_ecdf(aes(fitLB[metaDt$cellNames[metaDt$predictedGroupL2_Co == "B naive"], "k1"],
                color = "B naive"), lwd = 0.4) +
  stat_ecdf(aes(fitLB[metaDt$cellNames[metaDt$predictedGroupL2_Co == "B intermediate"], "k1"],
                color = "B intermediate"), lwd = 0.4) +
  stat_ecdf(aes(fitLB[metaDt$cellNames[metaDt$predictedGroupL2_Co == "B memory"], "k1"],
                color = "B memory"), lwd = 0.4) +
  stat_ecdf(aes(fitLB[metaDt$cellNames[metaDt$predictedGroupL2_Co == "Plasmablast"], "k1"],
                color = "Plasmablast"), lwd = 0.4) +
  scale_color_manual(limits = c("B naive", "B intermediate", "B memory", "Plasmablast"),
                     values = c("#fff143", "#ffc773", "#ffb61e", "#ffa400")) +
  labs(x = "k2: B naive", y = "") +
  theme(legend.position = c(0.7, 0.4))

ggsave(filename = "k1_Bcell_ecdf.pdf", width = 3, height = 3)

# Correlate with metadata ----
## Disease ----
metaDt <- getCellColData(projPBMC)
metaDt$cellNames <- rownames(metaDt)
metaDt <- metaDt[rownames(fit$L), ]
metaDt <- as.data.frame(metaDt)

outDf <- data.frame()
for (i in 1:20) {
  tmpL <- fit$L[, i]
  tmpDt <- metaDt %>% mutate(tmpL = tmpL[rownames(.)])

  tmpDt <- tmpDt %>%
    filter(Disease != "COVID-19 convalescent" & tmpL > 0.05) %>%
    mutate(Disease = factor(Disease, levels = c("Healthy", "COVID-19")))

  fitLmer <- lmer(scale(tmpL) ~ Disease + (1|Donor) + scale(log10(nFrags)) + scale(TSSEnrichment),
                  data = tmpDt, REML = F)

  fitLmer1 <- lmer(scale(tmpL) ~ (1|Donor) + scale(log10(nFrags)) + scale(TSSEnrichment),
                   data = tmpDt, REML = F)

  summ <- summary(fitLmer)

  lrt <- anova(fitLmer1, fitLmer)
  outDf <- rbind(outDf, data.frame(k = "k"%&%i, Estimate = summ$coefficients["DiseaseCOVID-19", "Estimate"],
                                   se = summ$coefficients["DiseaseCOVID-19", "Std. Error"],
                                   tval = summ$coefficients["DiseaseCOVID-19", "t value"],
                                   pval = lrt$`Pr(>Chisq)`[2]))
}

ggplot(outDf, aes(x = k, y = Estimate, color = twentyCol[1:20])) +
  theme_cowplot(font_size = 10, line_size = 0.4) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.4) +
  geom_point(aes(size = -log10(pval))) +
  scale_size(name = "", range = c(0.1, 2)) +
  geom_errorbar(aes(ymin = Estimate - se / 2, ymax = Estimate + se / 2),
                linewidth = 0.4, width = 0.2) +
  scale_x_discrete(limits = "k"%&%1:20) +
  scale_color_identity() +
  labs(x = "") +
  theme(panel.grid.major.x = element_line(color = "lightgrey", linetype = 3, linewidth = 0.2),
        legend.position = c(0.7, 0.3), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

ggsave("topic20_enrichment_fixed.pdf", width = 4, height = 2.4)

## TSSEnrichment ----
ggPoint(fit$L[, 2], metaDt$TSSEnrichment, colorDensity = T, ylabel = "TSSEnrichment",
        xlabel = "k2 loadings", rastr = T)
ggsave("k2_TSSEnrichment.pdf", width = 3, height = 3)

## TF Enrichment ----
fitFZ <- t(scale(t(fitF)))
fOrder <- apply(t(scale(t(fitFZ))), 2, function(x) length(x) - rank(x) + 1)
fOrder <- fOrder[rownames(fitFZ), ]
fQuantile <- fOrder / nrow(fOrder) * 100

topicSe <- SummarizedExperiment(
  assays = list(fOrder = fOrder, fZ = fitFZ, fQuantile = fQuantile),
  rowData = as.data.frame(peakSet[rownames(fitFZ)]),
  metadata = list(Params = list(useMatrix = "PeakMatrix"))
)

enrichMotifs1 <- peakAnnoEnrichment(
  seMarker = topicSe,
  ArchRProj = projPBMC,
  peakAnnotation = "cisbp",
  cutOff = "fOrder < 3000"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs1, n = 3, transpose = F, rastr = F)
draw(heatmapEM, heatmap_legend_side = "bottom", column_title_gp = gpar(fontsize = 6),
     annotation_legend_side = "top")

enrichP <- as.matrix(assay(enrichMotifs1, "mlog10Padj"))
# enrichP <- enrichP[rowMaxs(enrichP) > 3,  ]
tf2show <- unique(unlist(lapply(1:ncol(enrichP), function(x) names(sort(enrichP[, x], decreasing = T)[1:5]))))
tmpLab <- unique(unlist(lapply(1:ncol(enrichP), function(x) names(sort(enrichP[, x], decreasing = T)[1:4]))))
tmpLab <- c("PAX5_709", "SP6_275", "SP5_279", "KLF6_165", "SP2_232", "RUNX1_733", "TBX10_791", "EOMES_788",
            "RUNX2_732", "BCL11A_194", "BCL11B_825", "IRF4_632", "FOSL1_142", "JUND_124", "FOSL2_105",
            "EGR1_195", "TCF7L2_762", "LEF1_760", "TCF7L1_763", "TCF7_750", "NFKB1_719", "RELA_722",
            "REL_721", "CEBPA_155", "CEBPD_152", "CEBPB_140", "POU2F2_609", "POU5F1B_622", "POU2F3_613",
            "POU3F4_619", "POU2F1_614", "GATA2_388", "GATA1_383", "GATA3_384", "GATA5_385", "WT1_266",
            "ZBTB7A_258", "SMAD5_866", "DNMT1_301", "JUN_143", "JUNB_139", "TFAP2C_3", "ZNF263_159",
            "FOS_137", "RORA_658", "RORB_865", "TFAP2D_2", "BATF_129")
enrichP <- enrichP / rowMaxs(enrichP)
enrichP <- enrichP[tf2show, ]

mlog10Padj <- as.data.frame(assay(enrichMotifs1, "mlog10Padj")) %>%
  rownames_to_column("TF") %>%
  pivot_longer(cols = k1:k20, names_to = "topic", values_to = "minusLog10Padj")

Enrichment <- as.data.frame(assay(enrichMotifs1, "Enrichment")) %>%
  rownames_to_column("TF") %>%
  pivot_longer(cols = k1:k20, names_to = "topic", values_to = "Enrichment")

outDf <- mlog10Padj %>%
  inner_join(Enrichment, by = c("topic", "TF")) %>%
  dplyr::group_by(topic) %>%
  dplyr::slice_max(minusLog10Padj, n = 30) %>%
  dplyr::ungroup() %>%
  filter(minusLog10Padj >= 10) %>%
  dplyr::arrange(topic, desc(minusLog10Padj)) %>%
  mutate(TF = str_remove(TF, "_[0-9]+$"))

tmpLab <- c("AC0016|KLF|C2H2_ZF", "AC0039|TFAP2A/TFAP2C|AP-2", "AC0507|PAX|Paired_box",
            "AC0046|EBF|['Rel_homology_region_(RHR)_factors']", "AC0380|PAX|Homeodomain,Paired_box",
            "AC0017|RUNX|Runt", "AC0036|BCL11A/BCL11B|C2H2_ZF", "AC0631|T/TBX|T-box",
            "AC0430|TCF7L/LEF|Sox", "AC0208|NFATC/ZNF|Rel", "AC0239|BACH/NFE|bZIP",
            "AC0242|FOSL/JUND|bZIP", "AC0292|ZBTB|C2H2_ZF",
            "AC0083|ZBTB|C2H2_ZF", "AC0274|PAX|Homeodomain,Paired_box",
            "AC0408|RORC/RORA|Nuclear_receptor", "AC0628|SMAD|SMAD", "AC0635|TBX/EOMES|T-box",
            "AC0118|ZBTB7A|C2H2_ZF", "AC0410|RARA/RARG|Nuclear_receptor", "AC0414|NR4A/NR2C|Nuclear_receptor",
            "AC0196|NFKB/RELA|Rel", "AC0190|IRF/ZBED|IRF", "AC0227|SPI/BCL11A|Ets",
            "AC0005|POU3F/POU1F|Homeodomain,POU", "AC0009|PAX/VSX|Homeodomain",
            "AC0252|GATA/GATA1||TAL|GATA", "AC0123|FOXO|Forkhead", "AC0213|STAT/STAT5A|STAT", "AC0236|BATF/IRF|bZIP",
            "AC0210|BCL/BCL6B|C2H2_ZF", "AC0311|CTCF|C2H2_ZF", "AC0144|IKZF/NFKB|C2H2_ZF", "AC0170|FLI/ZNF|Ets")

tmpLab1 <- unique(unlist(lapply(1:ncol(enrichP), function(x) names(sort(enrichP[, x], decreasing = T))[1])))

tmpLab <- unique(c(tmpLab, tmpLab1))

annoTF <- rowAnnotation(
  gene = anno_mark(
    at = which(rownames(enrichP) %in% tmpLab),
    labels = str_remove(rownames(enrichP)[rownames(enrichP) %in% tmpLab], "_[0-9]+"),
    labels_gp = gpar(fontsize = unit(4, "pt")), link_gp = gpar(lwd = 0.5),
    link_width = unit(0.15, "in")
  )
)

myCol <- circlize::colorRamp2(breaks = c(0, 1), colors = c("white", "purple3"))
ht <- Heatmap(
  enrichP,
  show_row_names = F, right_annotation = annoTF, col = myCol,
  show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = T,
  clustering_method_rows = "ward.D2",
  clustering_distance_rows = "canberra",
  clustering_method_columns = "ward.D2",
  clustering_distance_columns = "canberra",
  column_names_gp = gpar(fontsize = 5),
  column_title_side = "top",
  border = "black",
  border_gp = gpar(lwd = 0.5),
  column_title = "-log10(adj. P) of TF enrichment",
  column_title_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(
    title = "", legend_direction = "horizontal", border = NULL,
    legend_width = unit(1, "in"), labels_gp = gpar(fontsize = 5),
    grid_height = unit(0.2, "in")
  ),
  height = unit(4, "in"), width = unit(0.07 * 20, "in"),
  use_raster = T,
  raster_quality = 5
)

pdf("topicTFEnrichment_cisbp.pdf", width = 3, height = 5)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

enrichP <- as.matrix(assay(enrichMotifs1, "mlog10Padj"))
as.data.frame(enrichP) %>%
  rownames_to_column("TF") %>%
  mutate(TF = str_remove(TF, "_[0-9]+")) %>%
  dplyr::arrange(desc(k20)) %>%
  mutate(x = row_number()) %>%
  ggplot(aes(x = x, y = k20)) +
  cowplot::theme_cowplot(font_size = 10, line_size = 0.4, font_family = "Helvetica") +
  geom_point(size = 1) +
  geom_text_repel(aes(label = TF), data = as.data.frame(enrichP) %>%
                    rownames_to_column("TF") %>%
                    mutate(TF = str_remove(TF, "_[0-9]+")) %>%
                    dplyr::arrange(desc(k20)) %>%
                    mutate(x = row_number()) %>% filter(k20 > 5),
                  xlim = c(10, NA), size = 3) +
  labs(x = "Rank order", y = "-log10(P) in k20")

ggsave("topicTFEnrichment_k20.pdf", width = 3, height = 3)

## Correlation with TF Z score ----
devMtrx <- getMatrixFromProject(projPBMC, "vierstraMatrix", threads = 5)
chromVarZ <- t(assay(devMtrx, "z"))
chromVarZ <- chromVarZ[rownames(fit$L), ]
ggPoint(fit$L[, 4], chromVarZ[, "BACH2"], colorDensity = T) +
  geom_smooth(method = "gam", color = "black", lty = 2) +
  geom_text(aes(x = Inf, y = Inf, label = "BACH2", vjust = 2, hjust = 2), color = "black")

g <- ggPoint(fit$L[, 2], chromVarZ[, "LEF1"], baseSize = 7,
             colorDensity = T, rastr = T, title = "LEF1",
             xlabel = "Loadings k2", ylabel = "chromVar Z score", legendSize = 0.1) +
  geom_smooth(method = "gam", color = "black", lty = 2, lwd = 0.3) +
  cowplot::theme_cowplot(font_size = 5, line_size = 0.3) +
  theme(legend.position = "none")

ggsave("loading_ChromVar_LEF1.pdf",
       width = unit(1.5, "in"), height = unit(1.5, "in"), plot = g)

g <- ggPoint(fit$L[, 4], chromVarZ[, "TBX21"], baseSize = 7,
             colorDensity = T, rastr = T, title = "TBX21",
             xlabel = "Loadings k4", ylabel = "chromVar Z score", legendSize = 0.1) +
  geom_smooth(method = "gam", color = "black", lty = 2, lwd = 0.3) +
  cowplot::theme_cowplot(font_size = 5, line_size = 0.3) +
  theme(legend.position = "none")

ggsave("loading_ChromVar_TBX21.pdf",
       width = unit(1.5, "in"), height = unit(1.5, "in"), plot = g)

corMtrx <- cor(as.matrix(chromVarZ), fit$L, method = "s")
tf <- c("TBX21", "TBR1", "EOMES", "RUNX1", "RUNX3", "ETS1", "ERG", "FOXP1", "ETV4", "IRF7", "GATA6", "GATA5", "NRF1")
ggplot(mapping = aes_(x = 1:nrow(corMtrx), y = sort(corMtrx[, 4], decreasing = T))) +
  theme_bw() +
  geom_col(aes(color = names(sort(corMtrx[, 4], decreasing = T)) %in% tf), fill = NA) +
  scale_color_manual(values = c("grey", "red")) +
  ggrepel::geom_text_repel(aes(
    x = which(names(sort(corMtrx[, 4], decreasing = T)) %in% tf),
    y = sort(corMtrx[tf, 4], decreasing = T),
    label = names(sort(corMtrx[tf, 4], decreasing = T))
  ),
  force = 10, max.overlaps = 50, force_pull = 0.1, nudge_x = 100
  ) +
  theme(legend.position = "none")

corMtrx <- cor(as.matrix(chromVarZ)[fit$L[, 2] > 0.1, ], fit$L[fit$L[, 2] > 0.1, 2], method = "s")
tf <- c("TCF7L2", "TCF7", "LEF1", "TCF7L1", "ZNF444", "ZSCAN1", "VEZF1", "ZBTB7B", "SNAI1")
ggplot(mapping = aes_(x = 1:nrow(corMtrx), y = sort(corMtrx[, 1], decreasing = T))) +
  theme_bw() +
  geom_col(aes(color = names(sort(corMtrx[, 1], decreasing = T)) %in% tf), fill = NA) +
  scale_color_manual(values = c("grey", "red")) +
  ggrepel::geom_text_repel(aes(
    x = which(names(sort(corMtrx[, 1], decreasing = T)) %in% tf),
    y = sort(corMtrx[tf, 1], decreasing = T),
    label = names(sort(corMtrx[tf, 1], decreasing = T))
  ),
  force = 30, max.overlaps = 50, force_pull = 0, xlim = c(10, 300),
  ) +
  theme(legend.position = "none")

corMtrx <- cor(as.matrix(chromVarZ)[fit$L[, 10] > 0.1, ], fit$L[fit$L[, 10] > 0.1, 10], method = "s")
tf <- c("TCF7L2", "TCF7", "LEF1", "TCF7L1", "ZEB1", "TCF4", "ZBTB14", "EGR4", "EGR1", "ZBTB7B",
        "EGR3", "EGR2", "HOXA9", "E2F5", "ZBTB7A", "KLF11")
ggplot(mapping = aes_(x = 1:nrow(corMtrx), y = sort(corMtrx[, 10], decreasing = T))) +
  theme_bw() +
  geom_col(aes(color = names(sort(corMtrx[, 10], decreasing = T)) %in% tf), fill = NA) +
  scale_color_manual(values = c("grey", "red")) +
  ggrepel::geom_text_repel(aes(
    x = which(names(sort(corMtrx[, 10], decreasing = T)) %in% tf),
    y = sort(corMtrx[tf, 10], decreasing = T),
    label = names(sort(corMtrx[tf, 10], decreasing = T))
  ),
  force = 30, max.overlaps = 50, force_pull = 0, xlim = c(10, 300),
  ) +
  theme(legend.position = "none")

Heatmap(corMtrx)

grpVierstraMtrxValueL_list <- pbmcapply::pbmclapply(cellGrpL, function(x) {
  return(sparseMatrixStats::colMeans2(chromVarZ, rows = x))
}, mc.cores = 15)

grpVierstraMtrxValueL <- matrix(unlist(grpVierstraMtrxValueL_list), byrow = T, nrow = length(grpVierstraMtrxValueL_list))
colnames(grpVierstraMtrxValueL) <- colnames(chromVarZ)
rownames(grpVierstraMtrxValueL) <- "G"%&%1:length(cellGrpL)

write_rds(grpVierstraMtrxValueL, "archR_3study_cmb_amulet_clean_v1/grpVierstraMtrxValueL.rds")

## Correlation with GA score ----
### Low-overlapping aggregates ----
grpGsMtrxValue <- ArchR:::.safelapply("chr"%&%1:22, function(x) {
  gsMtrxTmp <- getMatrixFromProject(projPBMC, "GeneScoreMatrix", threads = 4, useSeqnames = x)
  rownames(gsMtrxTmp) <- rowData(gsMtrxTmp)$name
  gsMtrxValueTmp <- assay(gsMtrxTmp, "GeneScoreMatrix")
  rm(gsMtrxTmp)
  gc()

  gsMtrxValueTmp <- gsMtrxValueTmp[sparseMatrixStats::rowSums2(gsMtrxValueTmp) > 0, ]
  gsMtrxValueTmp <- log1p(gsMtrxValueTmp)
  gsMtrxValueTmp <- t(gsMtrxValueTmp)
  gsMtrxValueTmp <- gsMtrxValueTmp[rownames(fit$L), ]

  tmpGrpGsMtrxValue <- ArchR:::.safelapply(cellGrp, function(x) {
    return(sparseMatrixStats::colMeans2(gsMtrxValueTmp, rows = x))
  }, threads = 3) %>% Reduce(rbind, .)

  colnames(tmpGrpGsMtrxValue) <- colnames(gsMtrxValueTmp)
  rownames(tmpGrpGsMtrxValue) <- "G"%&%1:length(cellGrp)

  rm(gsMtrxValueTmp)
  gc()

  return(tmpGrpGsMtrxValue)
}, threads = 3) %>% Reduce(cbind, .)

grpGsMtrxValueL <- ArchR:::.safelapply("chr"%&%1:22, function(x) {
  gsMtrxTmp <- getMatrixFromProject(projPBMC, "GeneScoreMatrix", threads = 4, useSeqnames = x)
  rownames(gsMtrxTmp) <- rowData(gsMtrxTmp)$name
  gsMtrxValueTmp <- assay(gsMtrxTmp, "GeneScoreMatrix")
  rm(gsMtrxTmp)
  gc()

  gsMtrxValueTmp <- gsMtrxValueTmp[sparseMatrixStats::rowSums2(gsMtrxValueTmp) > 0, ]
  gsMtrxValueTmp <- log1p(gsMtrxValueTmp)
  gsMtrxValueTmp <- t(gsMtrxValueTmp)
  gsMtrxValueTmp <- gsMtrxValueTmp[rownames(fit$L), ]

  tmpGrpGsMtrxValue <- ArchR:::.safelapply(cellGrpL, function(x) {
    return(sparseMatrixStats::colMeans2(gsMtrxValueTmp, rows = x))
  }, threads = 3) %>% Reduce(rbind, .)

  colnames(tmpGrpGsMtrxValue) <- colnames(gsMtrxValueTmp)
  rownames(tmpGrpGsMtrxValue) <- "G"%&%1:length(cellGrpL)

  rm(gsMtrxValueTmp)
  gc()

  return(tmpGrpGsMtrxValue)
}, threads = 3) %>% Reduce(cbind, .)

grpFitL <- pbmclapply(cellGrp, function(x) {
  return(colMeans2(fit$L[x, ], useNames = T))
}, mc.cores = 4) %>% Reduce(rbind, .)

colnames(grpFitL) <- colnames(fit$L)
rownames(grpFitL) <- "G"%&%1:length(cellGrp)

grpFitLL <- pbmclapply(cellGrpL, function(x) {
  return(colMeans2(fit$L[x, ], useNames = T))
}, mc.cores = 4) %>% Reduce(rbind, .)

colnames(grpFitLL) <- colnames(fit$L)
rownames(grpFitLL) <- "G"%&%1:length(cellGrpL)

grpGsMtrxValue1 <- grpGsMtrxValue[, colMaxs(grpGsMtrxValue, useNames = T) > 0.3]
grpGsMtrxValueL1 <- grpGsMtrxValueL[, colMaxs(grpGsMtrxValueL, useNames = T) > 0.3]

corMtrxNewL <- pbmclapply(1:ncol(grpFitLL), function(x) {
  message(str_glue("Calculating k{x}...\n"))
  tmpL <- grpFitLL[, x] > 0.01
  tmp <- cor(grpGsMtrxValueL1[tmpL, ], grpFitLL[tmpL, x], method = "s")
}, mc.cores = 2) %>% Reduce(cbind, .)

colnames(corMtrxNewL) <- colnames(grpFitLL)

corGeneListL <- lapply(1:ncol(corMtrxNewL), function(x) {
  tmp <- names(sort(corMtrxNewL[, x], decreasing = T))[1:50]
  return(tmp)
})

# Trajectory ----
## k1 ----
min.max <- function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
projPBMC$k1 <- fit$L[projPBMC$cellNames, 1]
projPBMC$k1[projPBMC$predictedGroupL1_Co != "B" | projPBMC$predictedGroupL2_Co == "Plasmablast"] <- NA
projPBMC$k1 <- 1 - min.max(projPBMC$k1)
projPBMC$k1rank <- rank(projPBMC$k1, na.last = "keep") / max(rank(projPBMC$k1, na.last = "keep"), na.rm = T) * 100

tmpDf <- getCellColData(projPBMC, c("Donor", "k1rank", "predictedGroupL2_Co")) %>%
  as.data.frame() %>%
  filter(!is.na(k1rank)) %>%
  mutate(bin = cut_width(k1rank, width = 1))

plotDf <- pbmclapply(levels(tmpDf$bin), function(x) {
  tmp <- tmpDf %>% filter(bin == x)
  tmpProp <- sum(tmp$predictedGroupL2_Co != "B naive") / nrow(tmp)
  data.frame(x, tmpProp)
}, mc.cores = 10) %>% Reduce(rbind, .) %>%
  mutate(rmean = zoo::rollmean(tmpProp, k = 5, fill = NA))

ggplot(plotDf, aes(x = x, y = rmean * 100, group = 1)) +
  theme_cowplot(font_size = 10, line_size = 0.4) +
  geom_line(col = "pink1", linewidth = 1) +
  labs(x = "k1 cutoff", y = "Non-naive B cells\nin bin (%)") +
  scale_x_discrete(limits = levels(tmpDf$bin)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(filename = "Bmemory_prop_k1.pdf", height = 1.2, width = 2.5)

k1Bin <- cut_interval(projPBMC$k1, n = 5)
projPBMC$k1Bin <- k1Bin

projPBMC <- addTrajectory(
  ArchRProj = projPBMC,
  name = "k1Traj",
  groupBy = "k1Bin",
  trajectory = levels(k1Bin),
  reducedDims = "reducedMNN",
  embedding = NULL,
  force = T
)

p <- plotTrajectory(projPBMC, trajectory = "k1Traj", embedding = "UMAPMNN",
                    colorBy = "cellColData", name = "k1Traj", plotAs = "points")

trajMM  <- getTrajectory(ArchRProj = projPBMC, name = "k1rank", threads = 4,
                         useMatrix = "GeneScoreMatrix", log2Norm = T)

rownames(trajMM) <- str_remove(rownames(trajMM), "^chr[0-9X]+:")
trajMM <- trajMM[rownames(corMtrxNewL), ]

test <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0.9,
                              returnMatrix = T, grepExclude = "RP[1-9]+|AC[.]")

geneOrder <- rownames(test)
rownames(test) <- ifelse(rownames(test) %in% geneB, rownames(test), "")
heatCol <- circlize::colorRamp2(
  breaks = seq(-2, 2, length.out = 9),
  colors = c("#3361A5", "#248AF3", "#14B3FF", "#88CEEF", "#C1D5DC", "#EAD397", "#FDB31A", "#E42A2A", "#A31D1D")
)

h1 <- Heatmap(test, col = heatCol, cluster_rows = F, cluster_columns = F, show_row_names = T, show_column_names = F,
              row_names_gp = gpar(fontsize = unit(5, "pt"), fontface = "italic"), use_raster = T, raster_quality = 5,
              width = unit(0.8, "in"), show_heatmap_legend = F)

h2 <- Heatmap(-corMtrxNewL[geneOrder, "k1", drop = F], cluster_rows = F, cluster_columns = F,
              show_row_names = F, show_column_names = F, use_raster = T, raster_quality = 5,
              width = unit(0.1, "in"), show_heatmap_legend = F,
              col = circlize::colorRamp2(breaks = c(-0.5, 0, 0.5), colors = c("blue", "white", "red")))

pdf("k1_Bcell_traj.pdf", width = 1.5, height = 2.1)
draw(h2 + h1)
dev.off()

## k17 (COVID-19) ----
projPBMC$Disease1 <- ifelse(projPBMC$Disease == "COVID-19 convalescent", NA, projPBMC$Disease)

min.max <- function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
projPBMC$k17 <- fit$L[projPBMC$cellNames, 17]

getCellColData(projPBMC, select = c("predictedGroupL1_Co", "k17")) %>%
  as.data.frame() %>%
  ggplot(., aes(y = k17, x = predictedGroupL1_Co)) +
  geom_boxplot()

projPBMC$k17[projPBMC$k17 <= 0.01 | !projPBMC$predictedGroupL1_Co %in% c("CD4 T", "CD8 T", "NK", "other T")] <- NA
projPBMC$k17 <- min.max(projPBMC$k17)
projPBMC$k17rank <- rank(projPBMC$k17, na.last = "keep") / max(rank(projPBMC$k17, na.last = "keep"), na.rm = T) * 100

tmpDf <- getCellColData(projPBMC, c("Donor", "k17rank", "Disease1")) %>%
  as.data.frame() %>%
  filter(!is.na(k17rank) & !is.na(Disease1)) %>%
  mutate(bin = cut_width(k17rank, width = 1))

plotDf <- pbmclapply(levels(tmpDf$bin), function(x) {
  tmp <- tmpDf %>% filter(bin == x)
  tmpProp <- sum(tmp$Disease1 == "COVID-19") / nrow(tmp)
  data.frame(x, tmpProp)
}, mc.cores = 3) %>% Reduce(rbind, .) %>%
  mutate(rmean = zoo::rollmean(tmpProp, k = 5, fill = NA))

ggplot(plotDf, aes(x = x, y = rmean * 100, group = 1)) +
  theme_cowplot(font_size = 10, line_size = 0.4) +
  geom_line(col = "tan", linewidth = 1) +
  labs(x = "k17 cutoff", y = "COVID-19 cells\nin bin (%)") +
  scale_x_discrete(limits = levels(tmpDf$bin)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(filename = "covid19_prop_k17_fixed.pdf",
       height = 1.2, width = 2.5)

projPBMC$k17Bin <- cut_interval(projPBMC$k17rank, n = 5, labels = 1:5)
projPBMC$k17Bin <- ifelse(is.na(projPBMC$k17Bin), NA, projPBMC$k17Bin%&%projPBMC$Disease1)

### Differential genes ----
projPBMC$new <- ifelse(is.na(projPBMC$k17) | is.na(projPBMC$Disease1), NA, projPBMC$Disease)

getCellColData(projPBMC, select = c("new", "k17")) %>%
  as.data.frame() %>%
  ggplot(., aes(y = k17, color = new)) +
  stat_ecdf()

allTest <- getMarkerFeatures(
  ArchRProj = projPBMC,
  useMatrix = "GeneScoreMatrix",
  groupBy = "new",
  bias = c("log10(nFrags)", "TSSEnrichment"),
  maxCells = 1e4,
  useGroups = "COVID-19",
  bgdGroups = "Healthy",
  useSeqnames = c("chr"%&%1:22),
  threads = 4
)

binList2 <- pbmclapply(2:5, function(i) {
  tmp <- getMarkerFeatures(
    ArchRProj = projPBMC,
    useMatrix = "GeneScoreMatrix",
    groupBy = "k17Bin",
    bias = c("log10(nFrags)", "TSSEnrichment"),
    maxCells = 1e4,
    useGroups = i%&%"COVID-19",
    bgdGroups = c("1Healthy", "1COVID-19"),
    useSeqnames = c("chr"%&%1:22),
    threads = 2
  )
}, mc.cores = 1)

allTestRes <- as.data.frame(getMarkers(allTest, cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.5")[[1]])
allTestRes$group <- "6"
binRes <- lapply(1:4, function(x) {
  tmp <- as.data.frame(getMarkers(binList2[[x]], cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.5")[[1]])
  tmp$group <- as.character(x)
  return(tmp)
}) %>% Reduce(rbind, .)

cmbDf <- rbind(allTestRes, binRes)
cmbDf %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(type = Log2FC > 0) %>%
  dplyr::count(type) %>%
  mutate(n = case_when(type == FALSE ~ -n, .default = n)) %>%
  ggplot(aes(as.factor(group), n, fill = type)) +
  cowplot::theme_cowplot(font_size = 10, line_size = 0.4) +
  geom_col() +
  scale_fill_brewer(palette = "Set1", limits = c(TRUE, FALSE)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = c("6", "1", "2", "3", "4"),
                   labels = c("All k17 cells", "(20,40]", "(40,60]", "(60,80]", "(80,100]")) +
  labs(x = "k17 quintile", y = "Number of DEGs") +
  coord_flip() +
  theme(legend.position = "none")

ggsave("k17_raDiffGene_barplot_fixed.pdf", height = 2, width = 3)

diffGeneList <- lapply(1:4, function(x) {
  tmp <- as.data.frame(getMarkers(binList2[[x]], cutOff = "FDR <= 0.05 & Log2FC >= 0.5")[[1]])
  tmp$name
  return(tmp$name)
})

names(diffGeneList) <- as.character(1:4)
diffGeneList[["6"]] <- allTestRes$name[allTestRes$Log2FC >= 0.5 & allTestRes$FDR <= 0.05]

goResTraj <- compareCluster(
  diffGeneList,
  fun = "enrichGO", OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL",
  pvalueCutoff = 0.1, qvalueCutoff = 0.1, readable = T,
  universe = allTest@elementMetadata$name
)

goResTrajDf <- as.data.frame(goResTraj)
goResTrajDf <- goResTrajDf %>%
  rowwise() %>%
  mutate(gRatio = eval(parse(text = GeneRatio)),
         BgRatioNew = eval(parse(text = BgRatio))) %>%
  mutate(OR = gRatio / BgRatioNew, Description = as.character(str_to_sentence(Description))) %>%
  filter(Count >= 5 & qvalue < 0.05)

go2plot <- c(
  "Mononuclear cell migration", "Chemokine-mediated signaling pathway", "Lymphocyte chemotaxis",
  "Leukocyte mediated immunity", "T cell activation",
  "Acute-phase response", "Regulation of signaling receptor activity"
)

goResTrajDf %>%
  filter(Description %in% go2plot & qvalue < 0.05) %>%
  ggplot(aes(Cluster, Description, col = -log10(qvalue), size = OR)) +
  cowplot::theme_cowplot(font_size = 10, line_size = 0.4, font_family = "Helvetica") +
  geom_point() +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35),
                   limits = go2plot) +
  scale_color_gradient2(low = "blue2", mid = "white", high = "red3", midpoint = 0) +
  scale_size_continuous(limits = c(1, 30), breaks = c(10, 20, 30), range = c(1, 9)) +
  labs(x = "", y = "GO BP") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("k17Traj_up-gene_GO_BP_fixed.pdf", width = 4.5, height = 3.5)

binRes <- binRes %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(rho = cor(Log2FC, as.integer(group), method = "s"))

geneHigh <- binRes %>%
  filter(group %in% 4 & Log2FC > 0 & rho > 0.6) %>%
  pull(name) %>% unique()

geneLow <- binRes %>%
  filter(group %in% 4 & Log2FC < 0 & rho < -0.6) %>%
  pull(name) %>% unique()

rownames(trajMM) <- str_remove(rownames(trajMM), "^chr[0-9X]+:")
trajMM <- trajMM[rownames(corMtrxNewL), ]
test <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0.8,
                              returnMatrix = T)

geneOrder <- rownames(test)
covidGene <- unique(binRes$name)
geneHigh <- geneHigh[!str_detect(geneHigh, "^RP11|^AP00|^AC[0-9]|^AC[0-9]|^CT")]
geneLow <- geneLow[!str_detect(geneLow, "^RP11|^AP00|^AC[0-9]|^AC[0-9]|^CT")]
enrichedGene <- unique(str_split(c("C4A|ORM1|ORM2|SAA1|VCAM1|PRG4|ORM1|ORM2|SAA1|CR2|KLF7|IRS2|DACT1|ZNF667|ZNF667-AS1|C4A|ORM1|ORM2|VCAM1|EFNA5|CCL4|NOG|MSC|TSHZ2|C4A|ORM1|VCAM1|CD22|CR2|KCNJ15|FCRL2|C4A|TMEM45B|TSHZ2|FBLN7|PROX1-AS1|KCNJ15|CCL5|VCAM1|ARRDC4|MB21D2|C4A|CCR5|IL4|IL13|NKG7|ORM1|ORM2|PTGIR|SAA1|CCL3|CCL4|CCL5|VCAM1|ZEB2|CXCL13|IL36A|CCR2|CCR5|IL4|IL13|SAA1"), pattern = fixed("|"))[[1]])
rownames(test) <- ifelse(rownames(test) %in% c(enrichedGene, "PRG4", "CCL3", "CCR2"), rownames(test), "")

heatCol <- circlize::colorRamp2(
  breaks = seq(-2, 2, length.out = 9),
  colors = c("#3361A5", "#248AF3", "#14B3FF", "#88CEEF", "#C1D5DC", "#EAD397", "#FDB31A", "#E42A2A", "#A31D1D")
)

h1 <- Heatmap(test, col = heatCol, cluster_rows = F, cluster_columns = F, show_row_names = T, show_column_names = F,
              row_names_gp = gpar(fontsize = unit(5, "pt"), fontface = "italic"), use_raster = T, raster_quality = 5,
              width = unit(0.75, "in"), show_heatmap_legend = F)

h2 <- Heatmap(corMtrxNewL[geneOrder, "k17", drop = F], cluster_rows = F, cluster_columns = F,
              show_row_names = F, show_column_names = F, use_raster = T, raster_quality = 5,
              width = unit(0.08, "in"), show_heatmap_legend = F,
              col = circlize::colorRamp2(breaks = c(-0.5, 0, 0.5), colors = c("blue", "white", "red")))

pdf("k17_traj_fixed.pdf", width = 1.5, height = 2.1)
draw(h2 + h1)
dev.off()

### Plot with pgt ----
for (smp in unique(projPBMC$Sample)) {
  getCellColData(projPBMC, select = c("Sample", "k17Bin"), drop = F) %>%
    as.data.frame() %>%
    filter(Sample == smp & !is.na(k17Bin)) %>%
    mutate(bc = str_remove(rownames(.), fixed(smp%&%"#")),
           grp = Sample%&%"_"%&%k17Bin) %>%
    dplyr::select(bc, grp) %>%
    fwrite(file = str_glue("group_bigwig/{smp}_k17Bin_bc_fixed.txt"),
           quote = F, row.names = F, col.names = F, sep = "\t")
}

smps <- c(1:5%&%"Healthy", 1:5%&%"COVID-19")
iniFile <- "fig2.ini"
for (i in enrichedGene) {
  if (file.exists(str_glue("k17_{i}.pdf"))) {
    next
  }

  message(i)
  tmp <- ArchR:::extendGR(gencode_anno$genes[gencode_anno$genes$symbol == i],
                          upstream = 5e3, downstream = 5e3)

  maxBW <- pbmclapply(smps, function (smp) {
    tmpScore <- rtracklayer::import(
      str_glue("group_bigwig/{smp}.bw"),
      which = tmp
    )
    return(max(tmpScore$score))
  }, mc.cores = 2) %>% unlist() %>% max()

  ini1 <- read.ini(iniFile)
  for (smp in smps) {
    ini1[[smp]]$max_value <- ceiling(maxBW * 1)
    ini1[[smp]]$number_of_bins <- width(tmp) %/% 25
  }

  write.ini(ini1, filepath = str_glue("group_bigwig/k17_{i}_{seqnames(tmp)}_{start(tmp)}_{end(tmp)}.ini"))

  system(
    str_glue("grep -w {i} gencode.v19.annotation.gtf | grep -v 'AS1' > group_bigwig/tmp.gtf")
  )

  cmd <- str_glue(
    "pgt ",
    "--tracks group_bigwig/",
    "k17_{i}_{seqnames(tmp)}_{start(tmp)}_{end(tmp)}.ini ",
    "--width 7 --region {seqnames(tmp)}:{start(tmp)}-{end(tmp)} ",
    "--outFileName k17_{i}_fixed.pdf"
  )

  system(cmd)
}

smps <- c(1:5%&%"Healthy", 1:5%&%"RA")
iniFile <- "fig2.ini"
for (i in peaks) {
  if (file.exists(str_glue("k17_{i}.pdf"))) {
    next
  }

  message(i)
  tmp <- ArchR:::extendGR(peakSet[i], upstream = 2.5e3, downstream = 2.5e3)

  maxBW <- pbmclapply(smps, function (smp) {
    tmpScore <- rtracklayer::import(
      str_glue("group_bigwig/{smp}.bw"),
      which = tmp
    )
    return(max(tmpScore$score))
  }, mc.cores = 2) %>% unlist() %>% max()

  ini1 <- read.ini(iniFile)
  for (smp in smps) {
    ini1[[smp]]$max_value <- ceiling(maxBW * 1.01)
    ini1[[smp]]$number_of_bins <- width(tmp) %/% 25
  }

  write.ini(ini1, filepath = str_glue("group_bigwig/k17_{i}_{seqnames(tmp)}_{start(tmp)}_{end(tmp)}.ini"))

  cmd <- str_glue(
    "pgt --tracks group_bigwig/",
    "k17_{i}_{seqnames(tmp)}_{start(tmp)}_{end(tmp)}.ini ",
    "--width 7 --region {seqnames(tmp)}:{start(tmp)}-{end(tmp)} ",
    "--outFileName k17_{i}.pdf"
  )

  system(cmd)
}

# Enrichment analysis using F matrix ----
TSS <- gencode_anno_list$TSS
TSS <- TSS[seqnames(TSS) != "chrY"]
tssOvPeak <- findOverlaps(TSS, peakSet)

fitFZ <- t(scale(t(fitF)))

tssList <- lapply(1:ncol(fitF), function(x) {
  message(str_glue("Running {x}...\n"))
  tssOvPeak <- findOverlaps(TSS, peakSet)
  tmpDf <- data.frame(
    tss = TSS$symbol[tssOvPeak@from],
    peak = names(peakSet)[tssOvPeak@to]
  ) %>%
    mutate(score = fitF[peak, x]) %>%
    arrange(desc(score))

  gl <- tmpDf$score
  names(gl) <- tmpDf$tss
  return(na.exclude(gl))
})

names(tssList) <- "k"%&%seq_along(tssList)

nearestList <- lapply(1:ncol(fitF), function(x) {
  message(str_glue("Running {x}...\n"))
  nearestIdx <- nearest(TSS, peakSet)
  tmpDf <- data.frame(
    SYMBOL = TSS$symbol,
    peak = names(peakSet)[nearestIdx]
  ) %>%
    dplyr::mutate(score = fitF[peak, x]) %>%
    dplyr::arrange(desc(score))

  gl <- tmpDf$score
  names(gl) <- tmpDf$SYMBOL
  return(na.exclude(gl))
})

names(nearestList) <- "k"%&%seq_along(nearestList)

coAccList <- lapply(1:ncol(fitF), function(x, y) {
  message(str_glue("Running {x}...\n"))
  tssOvPeak <- findOverlaps(TSS, peakSet)
  tmpDf <- data.frame(
    SYMBOL = TSS$symbol[tssOvPeak@from],
    queryName = names(peakSet)[tssOvPeak@to]
  ) %>%
    left_join(dplyr::select(cADf, queryName, subjectName, correlation), by = "queryName") %>%
    filter(!is.na(subjectName)) %>%
    mutate(subjectScore = y[subjectName, x]) %>%
    group_by(SYMBOL) %>%
    summarise(sumScore = sum(subjectScore, na.rm = T))

  tmpDf1 <- data.frame(
    SYMBOL = TSS$symbol[tssOvPeak@from],
    queryName = names(peakSet)[tssOvPeak@to]
  ) %>%
    mutate(queryScore = y[queryName, x]) %>%
    filter(!is.na(queryScore))

  gl <- tmpDf$sumScore
  names(gl) <- tmpDf$SYMBOL

  gl1 <- tmpDf1$queryScore
  names(gl1) <- tmpDf1$SYMBOL
  gl1[names(gl)] <- gl1[names(gl)] + gl

  return(sort(gl1, decreasing = T))
}, fitF)

names(coAccList) <- "k"%&%seq_along(coAccList)

geneScorePM <- scoreGeneByPeak(
  genes = gencode_anno_list$genes, peaks = peakSet, scoreMat = fitF,
  blacklist = getBlacklist(projPBMC), method = "sum"
)

geneScorePM <- as.matrix(geneScorePM)
geneScoreTop <- lapply(1:ncol(geneScorePM), function(x) {
  return(names(sort(geneScorePM[, x], decreasing = T))[1:500])
})

names(geneScoreTop) <- "k"%&%seq_along(geneScoreTop)

geneScorePMZ <- t(scale(t(geneScorePM)))
geneScoreList <- lapply(1:ncol(geneScorePM), function(x) {
  tmp <- sort(geneScorePM[, x], decreasing = T)[1:1500]
  tmp <- tmp[geneScorePMZ[names(tmp), x] >= 1.5]
  return(tmp)
})

names(geneScoreList) <- "k"%&%seq_along(geneScoreList)

lapply(geneScoreList, function(x) names(x)[1:15])
top15 <- unique(unlist(lapply(geneScoreList, function(x) names(x)[1:15])))
top15 <- top15[!str_detect(top15, "RP[0-9]+|AC[0-9+]")]
gene2anno <- c(
  "CXCR4", "FOXP1", "CD74", "TNFRSF13C", "RAD51B", "JUND", "ING1", "IGFLR1", "ZEB2", "FAM53B",
  "PRF1", "NFATC2", "ADRB2", "HIVEP3", "S1PR5", "GZMB", "CCDC88", "INPP4A", "SH2B3", "ZMIZ1",
  "RUNX1", "FOS", "ZFP36", "JUNB", "DUSP2", "DUSP1", "NR4A2", "TNFAIP3", "KDM6B", "HOOK2",
  "BACH2", "BCL11B", "CCR7", "LEF1", "ITPKB", "TCF7", "CD247", "TRGV9", "TRGV10", "ZBTB16",
  "CCL2", "ETV6", "RAPGEF1", "IL8", "IL1B", "RARA", "IFI30", "TNFSF9", "EBF1", "MAFB", "NEAT1",
  "LYN", "RXRA", "CLEC11A", "SOX4", "INPP5A", "CISH", "ZFP36L2", "ZFP36L1", "TRBC2", "CXCL3",
  "S100A10", "CXCL2", "EGR1", "JARID2", "IRF7", "KDM4B", "FYN", "ITGAE", "RORA", "BHLHE40",
  "TNFRSF25", "IL7R", "RUNX2")

score2plot <- geneScorePMZ[gene2anno, ]
myCol <- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("steelblue", "white", "coral"))
ht <- Heatmap(
  score2plot,
  col = myCol, show_row_names = T,
  clustering_method_rows = "ward.D2", show_row_dend = F,
  clustering_method_columns = "ward.D2", show_column_dend = F,
  column_names_gp = gpar(fontsize = 5),
  row_names_gp = gpar(fontsize = 5, fontface = "italic"),
  heatmap_legend_param = list(
    title = "", border = NULL,
    legend_height = unit(0.4, "in"), labels_gp = gpar(fontsize = 5), grid_height = unit(0.07, "in")
  ),
  height = unit(2.5, "in"),
  width = unit(0.06 * 20, "in"),
  use_raster = T,
  raster_quality = 5
)

pdf("topicGeneScorev1.pdf", height = 3, width = 3)
draw(ht, heatmap_legend_side = "right")
dev.off()

peakTypeDf <- as.data.frame(table(peakSet$peakType))
peakTypeDf$set <- "All"
for (i in 1:ncol(fitFZ)) {
  tmp <- as.data.frame(table(peakSet[names(sort(fitFZ[, i], decreasing = T))[1:1500]]$peakType))
  tmp$set <- "k"%&%i
  peakTypeDf <- rbind(peakTypeDf, tmp)
}

ggplot(peakTypeDf, aes(set, Freq, fill = Var1)) +
  theme_cowplot(font_size = 5, font_family = "Helvetica", line_size = 0.4) +
  geom_col(position = "fill") +
  scale_x_discrete(limits = c("All", "k"%&%1:20)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set1", name = "") +
  labs(x = "", y = "")

ggsave(width = 3.5, height = 1, filename = "topic20_peak_anno.pdf")

1 - phyper(q = 1259, m = 45568, n = 89127+21748+171303, k = 3000)

## Compare 4 peak2gene strategies ----
immunexUT <- fread("immunexUT_specific_genes.csv", sep = ",")

### k7 ----
utGene <- immunexUT %>%
  dplyr::select(`Cell type`, Gene) %>%
  filter(`Cell type` %in% c("Naïve CD4", "Naïve CD8")) %>%
  pull(Gene)

k7gene <- unique(c(azimuth_pbmc_L2$CD4_naive, azimuth_pbmc_L2$CD8_naive, utGene))
k7gene <- unique(c(azimuth_pbmc_L2$CD4_naive, azimuth_pbmc_L2$CD8_naive))
k7gene <- k7gene[k7gene %in% rownames(geneScorePM)]
k7Mtrx <- matrix(data = tssList$k7[k7gene], nrow = 1, dimnames = list("one", k7gene))
k7Mtrx <- rbind(k7Mtrx, matrix(data = nearestList$k7[k7gene], nrow = 1, dimnames = list("two", k7gene)))
k7Mtrx <- rbind(k7Mtrx, matrix(data = coAccList$k7[k7gene], nrow = 1, dimnames = list("three", k7gene)))
k7Mtrx <- rbind(k7Mtrx, matrix(data = geneScorePM[k7gene, "k7"], nrow = 1, dimnames = list("four", k7gene)))
k7Mtrx <- k7Mtrx[, colSums2(is.na(k7Mtrx), useNames = T) < 4]

k7Mtrx[k7Mtrx > 50] <- 50
gene2mark <- unique(c("LEF1", "TCF7", "CCR7", "IL7R", "CD4", "CD8A", "MAL", "LDHB", "NOSIP"))
ha <- columnAnnotation(gene = anno_mark(
  at = na.exclude(match(gene2mark, k7gene)), labels = gene2mark,
  labels_gp = gpar(fontsize = 5, fontface = "italic"), link_gp = gpar(lwd = 0.4),
  link_width = unit(0.1, "in"), side = "bottom"))

ht <- Heatmap(
  k7Mtrx / rowMaxs(k7Mtrx, na.rm = T, useNames = T), column_order = order(k7Mtrx[4, ], decreasing = T),
  bottom_annotation = ha,
  col = circlize::colorRamp2(c(0, 1), c("white", "red3")),
  cluster_rows = F, cluster_columns = F, clustering_method_columns = "ward.D2", show_column_dend = F, show_column_names = F,
  column_title_side = "top", column_title = "Normalized gene score (k7)", column_title_gp = gpar(fontsize = 5),
  show_row_names = F,
  height = unit(0.6, "in"), width = unit(1.5, "in"),
  border = "black", border_gp = gpar(lwd = 0.5),
  heatmap_legend_param = list(
    title = "", legend_direction = "horizontal",
    legend_width = unit(0.7, "in"), labels_gp = gpar(fontsize = 5),
    grid_height = unit(0.1, "in"), at = c(0, 1), labels = c("Min", "Max")
  ),
  use_raster = T
)

pdf("k7_score_heatmap.pdf")
draw(ht, heatmap_legend_side = "bottom")
dev.off()
