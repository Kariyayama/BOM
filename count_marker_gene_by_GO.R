# ライブラリ
library(rbioapi)
library(tidyverse)
library(SingleCellExperiment)
library(gridExtra)
library(cowplot)
library(here)
library(dplyr)
library(HGNChelper)
source("run_seurat_GO_analysis.R")
source("apply_sctype.R")
source("~/Labolatory/scRNA-seq/ScType/sc-type/R/gene_sets_prepare.R")
source("~/Labolatory/scRNA-seq/ScType/sc-type/R/sctype_score_.R")
db_ = "~/Labolatory/scRNA-seq/ScType/sc-type/master/ScTypeDB_full.xlsx";

# sctypeでラベル
run_sctype <- function(seurat, organ){
    db_ = "~/Labolatory/scRNA-seq/ScType/sc-type/ScTypeDB_full.xlsx"
    # get cell-type-specific gene sets from our in-built database (DB)
    gs_list = gene_sets_prepare(db_, organ)

    # get cell-type by cell matrix
    es.max = sctype_score(scRNAseqData = seurat@assays$SCT@scale.data, scaled = TRUE,
                          gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
    seurat@meta.data -> seurat.meta.data

    # merge by cluster
    clusters <- seurat$integrated_seurat_cluster
    cL_resutls = do.call("rbind",
                         lapply(unique(clusters),
                                function(cl){
                                    seurat.meta.data %>%
                                        filter(integrated_seurat_cluster == cl) %>%
                                        rownames -> es.max.list
                                    es.max[, es.max.list] %>%
                                        rowSums %>%
                                        sort(decreasing = !0) -> es.max.cl
                                    data.frame(cluster = cl,
                                               type = names(es.max.cl),
                                               scores = es.max.cl,
                                               ncells = sum(clusters==cl)) %>%
                                    head(10)
                                }
                         )
    )

    sctype_scores <- cL_resutls %>%
      group_by(cluster) %>%
      top_n(n = 1, wt = scores) %>%
      mutate(customclassif = paste(cluster, type, sep = "_"))

    # set low-confident (low ScType score) clusters to "unknown"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
    print(sctype_scores[,1:3])

    seurat@meta.data$customclassif = ""
    for(j in unique(sctype_scores$cluster)){
      cl_type = sctype_scores[sctype_scores$cluster==j,]
      seurat@meta.data$customclassif[seurat@meta.data$integrated_seurat_cluster == j] <- as.character(cl_type$customclassif[1])
    }
    return(seurat)
}

# ヒトのマーカー遺伝子数を対応関係ごとに合計
# マーカー遺伝子の読み込み
extract_Hs_count_gene <- function(organ){
    paste(curdir, "/", organ,"/Marker_gene/Hs_", organ,
            ".seurat_integrate.celltype_pick",
            target, "markergene.tsv", sep = "") %>%
        read_tsv(show_col_types = FALSE) -> Hs_all
    # 全マーカー遺伝子数の導出
    Hs_all %>%
      dplyr::count(cluster) -> Hs_all_count

    # マーカー遺伝子と遺伝子対抗関係を結合、遺伝子対抗関係ごとに合計数を計算
    Hs_all %>%
      inner_join(Hs_homology, by = c("gene" = "Gene_name")) %>%
      group_by(cluster) %>%
      dplyr::count(Mouse_homology_type.x) %>%
      mutate(SP = "Hs") %>%
      dplyr::rename(homology_type = Mouse_homology_type.x) %>%
      inner_join(Hs_all_count, by = "cluster") %>%
      mutate(ratio = n.x/n.y)
}

# マウスのマーカー遺伝子数を対応関係ごとに合計
# マーカー遺伝子の読み込み
extract_Mm_count_gene <- function(organ){
    paste(curdir, "/", organ,"/Marker_gene/Mm_", organ,
            ".seurat_integrate.celltype_pick",
            target, "markergene.tsv", sep = "") %>%
        read_tsv(show_col_types = FALSE) -> Mm_all
    # 全マーカー遺伝子数の導出
    Mm_all %>%
      dplyr::count(cluster) -> Mm_all_count

    # マーカー遺伝子と遺伝子対抗関係を結合、遺伝子対抗関係ごとに合計数を計算
    Mm_all %>%
      inner_join(Mm_homology, by = c("gene" = "Gene_name")) %>%
      group_by(cluster) %>%
      dplyr::count(Human_homology_type.x) %>%
      mutate(SP = "Mm") %>%
      dplyr::rename(homology_type = Human_homology_type.x) %>%
      inner_join(Mm_all_count, by = "cluster") %>%
      mutate(ratio = n.x/n.y)
}


