####################################################################################################
## Script purpose: COLOC GWAS with metasoft QTL in L1 cells
## Author: Zepeng Mu
####################################################################################################
library(coloc)
library(tidyverse)
library(data.table)
library(qvalue)
"%&%" <- function(a, b) paste0(a, b)

setDTthreads(4)

## Process GWAS
## Replace <dataDir> with your data dir
gwas <- fread("Trans_all_auto-10-2021_formatted.txt.gz",
              data.table = F)

loci <- fread("Trans_all_auto-loci-1e-07.txt")

## Load cell type to run COLOC
cellL1 <- commandArgs(trailingOnly = T)[1]

outDf <- data.frame(peak = character(), caQtlCell = character(), colocSnp = character(),
                    lociName = character(), nsnps = integer(), PP.H0.abf = double(), PP.H1.abf = double(),
                    PP.H2.abf = double(), PP.H3.abf = double(), PP.H4.abf = double(), testLead = character())

for (cell in cellL1) {
  ## Sample size
  N <- length(fread(str_glue("smp_to_use_L1_{cell}.txt"), header = F)$V1)

  ## QTL summ stats
  qtl <- fread(str_glue("poi_L1_{cell}_3study_noGC_formatted_250kb.txt.gz"), sep = "\t", data.table = F)

  qtlSig <- qtl %>%
    dplyr::filter(NUM_OF_STUDY > 1 & abs(BETA_RE) < 30) %>%
    dplyr::group_by(peak) %>%
    dplyr::mutate(padj = p.adjust(PVALUE_RE2, method = "bonferroni")) %>%
    dplyr::filter(PVALUE_RE2 == min(PVALUE_RE2, na.rm = T)) %>%
    dplyr::distinct(peak, .keep_all = T) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(qval = qvalue(padj)$qvalues) %>%
    dplyr::filter(qval < 0.1)

  for (i in unique(loci$CHR)) {
    chrLoci <- loci %>% dplyr::filter(CHR == i)
    chrGwas <- gwas %>%
      dplyr::filter(CHR == i) %>%
      dplyr::mutate(sid = CHR%&%":"%&%BP)

    ## MAF calculated from QTL samples
    chrMaf <- fread(
      str_glue("chr{i}.{cell}.filtered.AF.txt.gz"),
      header = F, col.names = c("chr", "pos", "ref", "alt", "af")
    ) %>%
      dplyr::mutate(sid = str_remove(chr, "chr")%&%":"%&%pos,
                    maf = case_when(af > 0.5 ~ 1 - af, T ~ af)) %>%
      dplyr::select(sid, maf)

    chrNomQtl <- qtl %>%
      dplyr::filter(NUM_OF_STUDY > 1 & chr == i) %>%
      dplyr::inner_join(chrMaf, by = "sid")

    for (j in 1:nrow(chrLoci)) {
      lociName <- chrLoci$SNPID[j]
      lociPos <- chrLoci$BP[j]
      upStream <- chrLoci$start[j]
      dnStream <- chrLoci$end[j]
      tmpGene <- qtlSig %>%
        dplyr::filter(chr == i & pos < dnStream & pos > upStream) %>%
        pull(peak) %>%
        unique()

      if (length(tmpGene) > 0) {
        for (g in tmpGene) {
          tmpQtl <- chrNomQtl %>% filter(peak == g & pos < dnStream & pos > upStream)
          tmpLeadPos <- qtlSig$pos[qtlSig$peak == g]

          tmpGwas <- chrGwas %>%
            dplyr::filter(BP %in% tmpQtl$pos & BP < dnStream & BP > upStream) %>%
            dplyr::distinct(sid, .keep_all = T)

          tmpQtl <- tmpQtl %>% filter(pos %in% tmpGwas$BP)

          ## Check if lead SNPs for both QTL and GWAS are in an overlapping region
          testLead <- ifelse(lociPos > min(tmpGwas$BP) & lociPos < max(tmpGwas$BP) & tmpLeadPos > min(tmpGwas$BP) & tmpLeadPos < max(tmpGwas$BP), "Yes", "No")

          if (nrow(tmpGwas) > 1) {
            message(g)

            my.res <- coloc.abf(
              dataset1 = list(
                N = N,
                pvalues = tmpQtl$PVALUE_RE2,
                MAF = tmpQtl$maf,
                type = "quant",
                snp = tmpQtl$sid
              ),
              dataset2 = list(
                beta = tmpGwas$Beta,
                varbeta = tmpGwas$SE ^ 2,
                type = "cc",
                s = 35871 / (240149 + 35871),
                snp = tmpGwas$sid
              )
            )

            colocSnp <- my.res$results$snp[which.max(my.res$results$SNP.PP.H4)]
            outDf <- rbind(outDf, data.frame(gene = g, caQtlCell = cell, colocSnp, lociName,
                                             nsnps = my.res$summary[1],
                                             PP.H0.abf = my.res$summary[2], PP.H1.abf = my.res$summary[3],
                                             PP.H2.abf = my.res$summary[4], PP.H3.abf = my.res$summary[5],
                                             PP.H4.abf = my.res$summary[6], testLead = testLead,
                                             row.names = NULL))
          }
        }
      }
    }
  }
}

fwrite(outDf, quote = F, col.names = T, row.names = F, sep = "\t",
       file = str_glue("poi_pval_Trans_all_5e5_{cellL1}_L1_metasoft.txt"))
