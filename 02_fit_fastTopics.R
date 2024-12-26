####################################################################################################
## Script purpose: Topic model analysis on three datasets using various K
## Author: Zepeng Mu
####################################################################################################
library(tidyverse)
library(data.table)
library(ArchR)
library(fastTopics)
library(Matrix)
library(matrixStats)
library(sparseMatrixStats)
library(GenomicRanges)
library(ArchRWrappers)
"%&%" <- function(a, b) paste0(a, b)
source("r_functions.R")

addArchRThreads(20)

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

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

# Get count matrix ----
peakCntMtrx <- getMatrixFromProject(
  ArchRProj = projPBMC,
  useMatrix = "PeakMatrix",
  useSeqnames = c("chr"%&%1:22, "chrX"),
  threads = 5
)

peakCntValue <- assay(peakCntMtrx)
rownames(peakCntValue) <- rowData(peakCntMtrx)$name
peakCntValue <- peakCntValue[sparseMatrixStats::rowSums2(peakCntValue) > 0, ]
peakCntValue <- t(peakCntValue)

rm(peakCntMtrx)

# fitPred ----
nMain <- 100
nRefine <- 200

if (file.exists("archR_3study_cmb_amulet_clean_v1/fastTopics/rndCell.rds")) {
  message("Using existing cells!")
  rndCellNew <- read_rds("archR_3study_cmb_amulet_clean_v1/fastTopics/rndCell.rds")
} else {
  rndCell <- sample(1:nrow(peakCntValue), 1e4, replace = F)
  subMtrx0 <- peakCntValue[rndCell, ]
  tmpPeak <- which(sparseMatrixStats::colSums2(subMtrx0) == 0)

  addCell <- c()
  for (i in tmpPeak) {
    tmpCell <- which(peakCntValue[, i] > 0)
    tmpCellNew <- sample(tmpCell, round(length(tmpCell) / 5))
    addCell <- append(addCell, tmpCellNew)
  }

  rndCellNew <- unique(c(rndCell, addCell))
  write_rds(rndCellNew, "archR_3study_cmb_amulet_clean_v1/fastTopics/rndCell.rds")
}

subMtrx <- peakCntValue[rndCellNew, ]

## Fit 6 topics ----
K <- 6
if (!file.exists("archR_3study_cmb_amulet_clean_v1/fastTopics/peak_fastTopics_fitCond0_K"%&%K%&%".rds")) {
  fit0 <- fit_topic_model(
    subMtrx,
    k = K,
    numiter.main = nMain,
    numiter.refine = nRefine,
    control.main = list(nc = 12),
    control.refine = list(nc = 12),
    verbose = "progressbar"
  )

  write_rds(fit0, "archR_3study_cmb_amulet_clean_v1/fastTopics/peak_fastTopics_fitCond0_K"%&%K%&%".rds")
}

if (!file.exists("archR_3study_cmb_amulet_clean_v1/fastTopics/peak_fastTopics_fitCond_K"%&%K%&%".rds")) {
  fit0 <- read_rds("archR_3study_cmb_amulet_clean_v1/fastTopics/peak_fastTopics_fitCond0_K"%&%K%&%".rds")
  predL <- predict(fit0, peakCntValue, control = list(nc = 12))
  fitPred <- init_poisson_nmf(X = peakCntValue, L = predL, init.method = "random", control = list(nc = 10))
  write_rds(fitPred, "archR_3study_cmb_amulet_clean_v1/fastTopics/peak_fastTopics_fitCond_K"%&%K%&%".rds")
}

## Fit more topics conditioning on first 6 ----
for (tmpK in c(8, 10, 12, 14, 16, 18, 20)) {
  message(str_glue("Running with {tmpK} K...\n"))
  if (!file.exists("archR_3study_cmb_amulet_clean_v1/fastTopics/peak_fastTopics_fitCond0_K"%&%tmpK%&%".rds")) {
    fitInit0 <- read_rds("archR_3study_cmb_amulet_clean_v1/fastTopics/peak_fastTopics_fitCond0_K"%&%(tmpK - 2)%&%".rds")
    uniMtrxL <- matrix(1 / tmpK, nrow = nrow(fitInit0$L), ncol = 2)
    uniMtrxF <- matrix(1 / nrow(fitInit0$F), nrow = nrow(fitInit0$F), ncol = 2)

    fitCond0 <- fitInit0
    fitCond0$L <- cbind(fitInit0$L, uniMtrxL)
    fitCond0$F <- cbind(fitInit0$F, uniMtrxF)
    colnames(fitCond0$L) <- colnames(fitCond0$F) <- "k"%&%1:tmpK

    fitNew0 <- init_poisson_nmf(subMtrx, F = fitCond0$F, L = fitCond0$L)

    fitOut0 <- fit_poisson_nmf(
      subMtrx,
      fit0 = fitNew0,
      numiter = nMain,
      method = "em",
      control = list(nc = 12),
    )

    fitOut0 <- fit_poisson_nmf(
      subMtrx,
      fit0 = fitOut0,
      method = "scd",
      numiter = nRefine,
      control = list(nc = 12),
    )

    fitOut0 <- poisson2multinom(fitOut0)
    write_rds(fitOut0, "archR_3study_cmb_amulet_clean_v1/fastTopics/peak_fastTopics_fitCond0_K"%&%tmpK%&%".rds")
  }

  if (!file.exists("archR_3study_cmb_amulet_clean_v1/fastTopics/peak_fastTopics_fitCond_K"%&%tmpK%&%".rds")) {
    fitOut0 <- read_rds("archR_3study_cmb_amulet_clean_v1/fastTopics/peak_fastTopics_fitCond0_K"%&%tmpK%&%".rds")
    predL <- predict(fitOut0, peakCntValue, control = list(nc = 5))
    fitPred <- init_poisson_nmf(X = peakCntValue, L = predL, init.method = "random")
    write_rds(fitPred, "archR_3study_cmb_amulet_clean_v1/fastTopics/peak_fastTopics_fitCond_K"%&%tmpK%&%".rds")
  }
}

message("Finished from R!")
