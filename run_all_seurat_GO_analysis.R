# 細胞クラスターで特異的に発現しているヒト特異的遺伝子とマウス特異的遺伝子のGO解析
renv::restore()
library(rbioapi)
library(tidyverse)
library(SingleCellExperiment)
library(gridExtra)
library(cowplot)
setwd("/Users/kariyayama/Labolatory/JST_Mouse_Human/Tabula/Integration")
source("run_seurat_GO_analysis.R")

input <- commandArgs(trailingOnly = TRUE)[1]
outdir <- paste(curdir, input, sep = "/")
test_type <- "FISHER"
cutoff <- 0.05
hsapiens <- 9606
mmuculus <- 10090
# BP <- "GO:0008150"
# BP: GO:0008150
# MF: GO:0003674
# CC: GO:0005575
BP <- "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP"

# 全てのIdentに対して実行
run_GOanalysis_for_all_ident("Hs", hsapiens)
run_GOanalysis_for_all_ident("Mm", mmuculus)

# GO解析の結果をプロットして可視化
# 定数
threshold <- -log10(0.05)
plot_para <- element_text(angle = 90, vjust = 0.5, hjust = 0)
height = 16
width <- function(x)  x * 0.05 + 5

# ヒトのオーソログがある遺伝子と種特異的遺伝子間の比較
plot_SS_ortholog_compare("Hs")
plot_SS_ortholog_compare("Mm")

