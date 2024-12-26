####################################################################################################
## Script purpose: scQTL in L1 cells with Poisson regression
## Author: Zepeng Mu
####################################################################################################
library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(lme4)
library(pbmcapply)
"%&%" <- function(a, b) paste0(a, b)

RhpcBLASctl::blas_set_num_threads(1)

# Set up parameters ----
windowSize <- 250e3 / 2 # A 250Kb window around peak center
chrm <- commandArgs(trailingOnly = T)[1]
inChunk <- commandArgs(trailingOnly = T)[2]
totalChunk <- as.integer(commandArgs(trailingOnly = T)[3])
mtrx <- commandArgs(trailingOnly = T)[4]
vcf <- commandArgs(trailingOnly = T)[5]
cell <- commandArgs(trailingOnly = T)[6]
outName <- commandArgs(trailingOnly = T)[7]

# List of candidate peaks ----
allRes <- fread("10kb_chr1-22_3study_onetime_allCell_hornet_6PC_lead_pval.txt.gz") %>%
  filter(emp.q.val < 0.1) %>%
  pull(peak)

pcL1 <- fread("nPC.txt") %>% filter(cell != "other")
cellRes <- lapply(pcL1$cell, function(x) {
  pc <- pcL1$nPC[pcL1$cell == x]
  tmp <- fread(str_glue("10kb_chr1-22_3study_onetime_L1_{x}_hornet_{pc}PC_lead_pval.txt.gz")) %>%
    filter(emp.q.val < 0.1) %>%
    pull(peak)

  return(tmp)
}) %>% unlist() %>% unique()

cmbPeak <- unique(c(allRes, cellRes))
rm(allRes, cellRes)
gc()

## Prepare matrices ----
smp <- fread("harmonised_smp_3study.txt", header = F)$V1
smpMeta <- fread("summary_patients_w_admixture_and_captureInfo_edRAGNew.txt") %>%
  filter(name_sequencing_scATACseq != "" & ID_VIW_PBMC != "09334-20")

smpCovid <- smpMeta$sample_ID
ID_Covid <- smpMeta$ID_VIW_PBMC

peakMtrx <- read_rds(mtrx)
peakMtrx <- peakMtrx[intersect(rownames(peakMtrx), cmbPeak), peakMtrx$Sample %in% ID_Covid]
metaDt <- as.data.frame(colData(peakMtrx))
lsi <- read_rds("IterativeLSI_3study.rds")
lsi <- lsi[rownames(metaDt), ]
gc()

### Scale and normalize counts ----
peakGr <- rowRanges(peakMtrx)
peakGr <- peakGr[!(seqnames(peakGr) == "chr6" & start(peakGr) > 25e6 & end(peakGr) < 35e6)]

peakGr$chunk <- cut(1:length(peakGr), breaks = totalChunk, labels = "chunk"%&%1:totalChunk)
peakGr <- peakGr[peakGr$chunk == inChunk]
peakMtrx <- peakMtrx[names(peakGr)]
peakCnt <- assay(peakMtrx)
peakCnt <- t(peakCnt)

# Load genotype PC ----
genoPC <- fread("1kg_3study_chr1-22.merged.1kgfmt.pruned.eigenvec", data.table = F) %>%
  dplyr::select("Sample" = V2, "PC1" = V3, "PC2" = V4, "PC3" = V5, "PC4" = V6, "PC5" = V7) %>%
  filter(Sample %in% metaDt$Sample)

rownames(genoPC) <- genoPC$Sample
genoPC <- scale(genoPC[, "PC"%&%1:5])
genoPC <- genoPC[metaDt$Sample, ]

metaDt <- metaDt %>%
  cbind(genoPC) %>%
  cbind(scale(lsi[, "LSI"%&%1:5])) %>%
  mutate(TSSEnrichment = scale(TSSEnrichment)[, 1], MTRatio = scale(MTRatio)[, 1],
         nFrags = scale(log10(nFrags))[, 1])

# Mapping sc-caQTL ----
form <- cnt ~ DS + PC1 + PC2 + PC3 + PC4 + PC5 + LSI1 + LSI2 + LSI3 + LSI4 + LSI5 + MTRatio + TSSEnrichment + nFrags + (1|Sample) + (1|Donor)
outDf <- lapply(1:length(peakGr), function(i) {
  tmpPeak <- names(peakGr)[i]
  tmpGr <- resize(peakGr[tmpPeak], windowSize * 2, fix = "center")
  start(tmpGr) <- max(start(tmpGr), 0)

  tmpDt <- metaDt %>% cbind(cnt = peakCnt[, tmpPeak])

  tmpDS <- fread(
    cmd = str_glue("tabix -h {vcf} chr{chrm}:{start(tmpGr)}-{end(tmpGr)}"),
    data.table = F, header = T, col.names = c("chr", "pos", "ref", "alt", smp)
  )

  tmpDS <- tmpDS %>%
    distinct(pos, .keep_all = T) %>%
    filter(!(chr == "chr6" & pos > 25e6 & pos < 35e6)) %>%
    dplyr::select(chr, pos, ref, alt, !!ID_Covid)

  tmpRes <- pbmclapply(tmpDS$pos, function(snp) {
    tmpDt$DS <- as.double(tmpDS[tmpDS$pos == snp, tmpDt$Sample])
    maf <- min(sum(round(tmpDt$DS)) / nrow(tmpDt) * 0.5, 1 - sum(round(tmpDt$DS)) / nrow(tmpDt) * 0.5)
    ref <- tmpDS$ref[tmpDS$pos == snp]
    alt <- tmpDS$alt[tmpDS$pos == snp]

    tmpTbl <- table(round(tmpDt$DS), tmpDt$cnt)
    if(nrow(tmpTbl) > 1 & ncol(tmpTbl) > 1 & min(table(round(as.double(tmpDS[tmpDS$pos == snp, ID_Covid])))) >= 2) {

      tmpOut <- tryCatch(
        {
          res <- glmer(
            formula = form, data = tmpDt, family = poisson, nAGQ = 0,
            control = glmerControl(optimizer = "bobyqa", calc.derivs = F)
          )

          resSumm <- summary(res)
          data.frame(
            peak = tmpPeak, chr = as.integer(chrm), snp.pos = snp, ref = ref,
            alt = alt, maf = maf, sdY = sd(tmpDt$cnt), N = nrow(tmpDt),
            effect = resSumm$coefficients["DS", "Estimate"],
            se = resSumm$coefficients["DS", "Std. Error"],
            z = resSumm$coefficients["DS", "z value"],
            pval = resSumm$coefficients["DS", "Pr(>|z|)"]
          )
        },
        error = function(e) { # Return NA if model not converge
          data.frame(
            peak = tmpPeak, chr = as.integer(chrm), snp.pos = snp, ref = ref,
            alt = alt, maf = maf, sdY = sd(tmpDt$cnt), N = nrow(tmpDt),
            effect = NA, se = NA, z = NA, pval = NA
          )
        }
      )

      return(tmpOut)
    }
  }, mc.cores = 2) %>% Reduce(rbind, .)

  return(tmpRes)

}) %>% Reduce(rbind, .)

fwrite(outDf, outName, quote = F, sep = "\t", col.names = F, row.names = F)
