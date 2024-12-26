####################################################################################################
## Script purpose: Interaction with topics for lead caQTLs
## Author: Zepeng Mu
####################################################################################################
library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(lme4)
library(pbmcapply)
library(qvalue)
library(patchwork)
"%&%" <- function(a, b) paste0(a, b)
source("r_functions.R")

# Set up parameters ----
myArgs <- commandArgs(trailingOnly = T)
chrm <- myArgs[1]
mtrx <- myArgs[2]
vcf <- myArgs[3]
cell <- myArgs[4]
outFile <- myArgs[5]
smp <- fread("harmonised_smp_3study.txt", header = F)$V1
sigPeak <- fread("sigPeak.txt", header = F)$V1 # A list of significant peaks

fit <- read_rds("peak_fastTopics_fitCond_K20.rds")
fitL <- fit$L

qtl <- fread(str_glue("poi_L1_{cell}_3study_noGC_formatted_250kb_lead.txt.gz"),
             sep = "\t", data.table = F) %>%
  dplyr::filter(qval < 0.1)

outDf <- data.frame()

## Prepare matrices ----
message(str_glue("Running chr{chrm}...\n"))
peakMtrx <- read_rds(mtrx)
rownames(peakMtrx) <- rowData(peakMtrx)$name
peakMtrx <- peakMtrx[intersect(rownames(peakMtrx), qtl$peak), ]
peakCnt <- assay(peakMtrx, "PeakMatrix")
peakCnt <- t(peakCnt)
metaDt <- as.data.frame(colData(peakMtrx))

peakGr <- rowRanges(peakMtrx)

rm(peakMtrx)
gc()

lsi <- read_rds("IterativeLSI_3study.rds")
lsi <- lsi[rownames(metaDt), ]

# Load genotype PC ----
genoPC <- fread("1kg_3study_chr1-22.merged.1kgfmt.pruned.eigenvec", data.table = F) %>%
  dplyr::select("Sample" = V2, "PC1" = V3, "PC2" = V4, "PC3" = V5, "PC4" = V6, "PC5" = V7) %>%
  filter(Sample %in% metaDt$Sample)

rownames(genoPC) <- genoPC$Sample
genoPC <- scale(genoPC[, "PC" %&% 1:5])
genoPC <- genoPC[metaDt$Sample, ]

## Define topic to test for each cell type
kList <- list("B" = c(1, 11), "CD4.T" = c(6, 7, 17), "CD8.T" = c(3, 6, 7, 14, 17, 18, 19),
              "NK" = c(17), "Mono" = c(10, 12, 15), "DC" = c(4, 10, 12, 15),
              "other.T" = c(3, 6, 8, 14, 17, 18, 19))

for (k in kList[[cell]]) {
  message(str_glue("Running k{k}...\n"))
  tmpDt <- metaDt %>%
    cbind(genoPC) %>%
    cbind(tmpL = fitL[rownames(.), k]) %>%
    dplyr::filter(tmpL > 0.05) %>%
    cbind(scale(lsi[rownames(.), "LSI" %&% 1:5])) %>%
    dplyr::mutate(TSSEnrichment = scale(TSSEnrichment)[, 1], MTRatio = scale(MTRatio)[, 1],
                  nFragsNew = scale(log10(nFrags))[, 1], scaleL = scale(tmpL)[, 1])

  tmpOut <- pbmclapply(names(peakGr), function(tmpPeak) {
    snp <- qtl %>% filter(peak == tmpPeak) %>% pull(pos)
    tmpGr <- peakGr[tmpPeak]
    chrm <- str_remove(as.character(seqnames(tmpGr)), "chr")

    tmpDS <- fread(cmd = str_glue("tabix -h {vcf} chr{chrm}:{snp-1}-{snp+1}"),
                   data.table = F, header = F,
                   col.names = c("chr", "pos", "ref", "alt", smp))

    ref <- tmpDS$ref[tmpDS$pos == snp]
    alt <- tmpDS$alt[tmpDS$pos == snp]

    tmpDt <- tmpDt %>% mutate(pheno = peakCnt[rownames(.), tmpPeak])

    tmpDt$DS <- as.double(tmpDS[tmpDS$pos == snp, tmpDt$Sample])
    tmpDt <- tmpDt %>%
      dplyr::select(pheno, tmpL, scaleL, DS, PC1:PC5, LSI1:LSI5, MTRatio, TSSEnrichment,
                    nFrags, nFragsNew, Sample, Donor, Study)

    tmpRes <- tryCatch(
      expr = {
        res <- glmer(
          pheno ~ DS + PC1 + PC2 + PC3 + PC4 + PC5 + LSI1 + LSI2 + LSI3 + LSI4 + LSI5 + MTRatio +
            TSSEnrichment + nFragsNew + (1 | Sample) + (1 | Donor) + (1 | Study) + scaleL + DS:scaleL,
          data = tmpDt, family = poisson, nAGQ = 0, control = glmerControl(optimizer = "bobyqa", calc.derivs = F)
        )

        res0 <- glmer(
          pheno ~ DS + PC1 + PC2 + PC3 + PC4 + PC5 + LSI1 + LSI2 + LSI3 + LSI4 + LSI5 + MTRatio +
            TSSEnrichment + nFragsNew + (1 | Sample) + (1 | Donor) + (1 | Study) + scaleL,
          data = tmpDt, family = poisson, nAGQ = 0, control = glmerControl(optimizer = "bobyqa", calc.derivs = F)
        )

        resSumm <- summary(res)$coefficients
        cmpRes <- anova(res0, res)

        data.frame(
          cell = cell, topic = k, peak = tmpPeak, chr = as.integer(chrm), snp.pos = snp, ref = ref,
          alt = alt, N = nrow(tmpDt),
          effect = resSumm["DS", "Estimate"],
          se = resSumm["DS", "Std. Error"],
          z = resSumm["DS", "z value"],
          p.chisq = cmpRes$`Pr(>Chisq)`[2],
          pval = resSumm["DS", "Pr(>|z|)"],
          effectL = resSumm["DS:scaleL", 1],
          seL = resSumm["DS:scaleL", 2],
          zL = resSumm["DS:scaleL", "z value"],
          pvalL = resSumm["DS:scaleL", 4]
        )
      },
      error = function(e) {}
    )
  }, mc.cores = 2) %>% Reduce(rbind, .)

  outDf <- rbind(outDf, tmpOut)
}

fwrite(outDf, sep = "\t", quote = F, col.names = T, row.names = F, file = outFile)
